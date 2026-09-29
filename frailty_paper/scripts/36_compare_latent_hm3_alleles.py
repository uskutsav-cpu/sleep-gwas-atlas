#!/usr/bin/env python3
"""Audit latent-frailty alleles against the pinned EUR HapMap3 GRCh37 map.

This is read-only QC. It compares rsID, chromosome, position and allele pairs
for the intersection with the versioned HM3 map. It does not harmonize effects,
drop variants, infer ambiguous palindromic orientation, or establish build for
variants outside that map.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys


ACCESSIONS = [f"GCST906240{i:02d}" for i in range(46, 53)]
ALLELES = set("ACGT")
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
ORIENTATION_FIELDS = [
    "same_orientation",
    "swapped_orientation",
    "complemented_same",
    "complemented_swapped",
    "palindromic_ambiguous",
    "allele_mismatch",
]
FIELDS = [
    "accession",
    "trait",
    "source_file",
    "source_sha256",
    "source_rows",
    "hm3_rsids_mapped",
    "hm3_coordinate_match",
    "hm3_coordinate_mismatch",
    "hm3_rsids_unmapped",
    "missing_rsid",
    *ORIENTATION_FIELDS,
    "qc_interpretation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def classify_alleles(effect: str, other: str, a1: str, a2: str) -> str:
    """Return strand/orientation class; palindrome never implies a strand."""
    ea, oa = effect.strip().upper(), other.strip().upper()
    ref1, ref2 = a1.strip().upper(), a2.strip().upper()
    if any(len(a) != 1 or a not in ALLELES for a in (ea, oa, ref1, ref2)):
        return "allele_mismatch"
    if ea == oa:
        return "allele_mismatch"
    if frozenset((ea, oa)) in PALINDROMIC:
        return "palindromic_ambiguous"
    if (ea, oa) == (ref1, ref2):
        return "same_orientation"
    if (ea, oa) == (ref2, ref1):
        return "swapped_orientation"
    comp = (COMPLEMENT[ea], COMPLEMENT[oa])
    if comp == (ref1, ref2):
        return "complemented_same"
    if comp == (ref2, ref1):
        return "complemented_swapped"
    return "allele_mismatch"


def load_map(path: Path, provenance_path: Path) -> tuple[dict[str, tuple[str, int, str, str]], dict]:
    provenance = json.loads(provenance_path.read_text())
    if provenance.get("schema_version") != "atlas.hm3-grch37-variant-map.v1":
        raise ValueError("unexpected HM3 map schema")
    if provenance.get("genome_build") != "GRCh37/hg19":
        raise ValueError("HM3 map is not documented as GRCh37/hg19")
    map_hash = sha256(path)
    if provenance.get("map_sha256") != map_hash:
        raise ValueError("HM3 map checksum does not match its provenance sidecar")
    mapping: dict[str, tuple[str, int, str, str]] = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8") as stream:
        rows = csv.DictReader(stream, delimiter="\t")
        if rows.fieldnames != ["SNP", "CHR", "BP", "A1", "A2"]:
            raise ValueError(f"unexpected HM3 map header: {rows.fieldnames}")
        for row in rows:
            rsid = row["SNP"].strip().lower()
            chrom = row["CHR"].strip()
            pos = int(row["BP"])
            a1, a2 = row["A1"].strip().upper(), row["A2"].strip().upper()
            if chrom not in {str(n) for n in range(1, 23)} or pos <= 0:
                raise ValueError(f"invalid HM3 coordinate for {rsid}: {chrom}:{pos}")
            if a1 not in ALLELES or a2 not in ALLELES or a1 == a2:
                raise ValueError(f"invalid HM3 alleles for {rsid}: {a1}/{a2}")
            if rsid in mapping:
                raise ValueError(f"duplicate HM3 rsID: {rsid}")
            mapping[rsid] = (chrom, pos, a1, a2)
    if len(mapping) != int(provenance.get("mapped_rows", -1)):
        raise ValueError("HM3 row count differs from provenance")
    return mapping, provenance


def scan_source(path: Path, repo: Path, accession: str, trait: str,
                reference: dict[str, tuple[str, int, str, str]]) -> dict[str, str | int]:
    counts = {field: 0 for field in FIELDS if field not in {
        "accession", "trait", "source_file", "source_sha256", "qc_interpretation"
    }}
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"rs_id", "chromosome", "base_pair_location", "effect_allele", "other_allele"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{accession} missing columns: {sorted(required-set(reader.fieldnames or []))}")
        for row in reader:
            counts["source_rows"] += 1
            rsid = (row.get("rs_id") or "").strip().lower()
            if not rsid or rsid in {"na", "nan", "."}:
                counts["missing_rsid"] += 1
                continue
            ref = reference.get(rsid)
            if ref is None:
                counts["hm3_rsids_unmapped"] += 1
                continue
            counts["hm3_rsids_mapped"] += 1
            chrom = (row.get("chromosome") or "").strip().upper().removeprefix("CHR")
            try:
                pos = int(row["base_pair_location"])
            except (ValueError, TypeError):
                pos = -1
            ref_chrom, ref_pos, a1, a2 = ref
            if chrom != ref_chrom or pos != ref_pos:
                counts["hm3_coordinate_mismatch"] += 1
                continue
            counts["hm3_coordinate_match"] += 1
            category = classify_alleles(row["effect_allele"], row["other_allele"], a1, a2)
            counts[category] += 1
    if counts["hm3_coordinate_mismatch"] or counts["allele_mismatch"]:
        verdict = "MISMATCHES_REVIEW_REQUIRED"
    elif counts["palindromic_ambiguous"]:
        verdict = "NONPALINDROMIC_COMPATIBLE; PALINDROMIC_UNRESOLVED"
    else:
        verdict = "ALL_MAPPED_HM3_ROWS_COORDINATE_AND_ALLELE_COMPATIBLE"
    return {
        "accession": accession,
        "trait": trait,
        "source_file": path.relative_to(repo).as_posix(),
        "source_sha256": sha256(path),
        **counts,
        "qc_interpretation": verdict,
    }


def write_immutable(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text() != content:
            raise FileExistsError(f"refusing to replace non-identical report: {path}")
        return
    path.write_text(content)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--reference-map", type=Path, required=True)
    parser.add_argument("--reference-provenance", type=Path)
    parser.add_argument("--output", type=Path, default=Path("frailty_paper/manifests/gwas_latent_hm3_allele_compatibility.tsv"))
    parser.add_argument("--provenance", type=Path, default=Path("frailty_paper/manifests/gwas_latent_hm3_allele_compatibility.provenance.json"))
    args = parser.parse_args()
    repo = args.repo.resolve()
    reference_map = args.reference_map.resolve()
    provenance_path = args.reference_provenance.resolve() if args.reference_provenance else Path(str(reference_map) + ".provenance.json")
    mapping, ref_provenance = load_map(reference_map, provenance_path)
    acquisition_manifest = repo / "frailty_paper/manifests/all_acquired_resources.tsv"
    acquisition_rows = {
        row["resource_id"]: row
        for row in csv.DictReader(acquisition_manifest.open(newline=""), delimiter="\t")
    }
    acquisition_manifest_sha256 = sha256(acquisition_manifest)
    factor_index = {
        row["source_id"]: row
        for row in csv.DictReader((repo / "frailty_paper/config/secondary_frailty_factors_v1.tsv").open(newline=""), delimiter="\t")
    }
    rows = []
    source_hashes = {}
    for accession in ACCESSIONS:
        meta = factor_index.get(accession)
        if meta is None or meta.get("build") != "GRCh37":
            raise ValueError(f"{accession} lacks a frozen GRCh37 factor mapping")
        trait = meta["trait_id"]
        path = repo / "frailty_paper/data/gwas/latent_frailty_catalog" / accession / f"{accession}.tsv"
        if not path.is_file():
            raise FileNotFoundError(path)
        registered = acquisition_rows.get(f"gwas_catalog_{accession}")
        if registered is None or registered.get("genome_build") != "GRCh37":
            raise ValueError(f"{accession} is absent from the acquired-resource manifest or not GRCh37")
        if registered.get("file") != path.relative_to(repo).as_posix():
            raise ValueError(f"{accession} source path differs from the acquired-resource manifest")
        row = scan_source(path, repo, accession, trait, mapping)
        if int(registered.get("bytes", -1)) != path.stat().st_size:
            raise ValueError(f"{accession} source byte count differs from the acquired-resource manifest")
        if row["source_sha256"] != registered.get("sha256"):
            raise ValueError(f"{accession} source SHA-256 differs from the acquired-resource manifest")
        rows.append(row)
        source_hashes[accession] = row["source_sha256"]
        print(f"{accession} rows={row['source_rows']} HM3_ID={row['hm3_rsids_mapped']} coordinate_mismatch={row['hm3_coordinate_mismatch']} allele_mismatch={row['allele_mismatch']} palindromic={row['palindromic_ambiguous']} verdict={row['qc_interpretation']}", flush=True)

    output = args.output if args.output.is_absolute() else repo / args.output
    provenance_out = args.provenance if args.provenance.is_absolute() else repo / args.provenance
    import io
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    table_text = buffer.getvalue()
    write_immutable(output, table_text)
    provenance = {
        "schema_version": "frailty_latent_hm3_allele_compatibility.1",
        "analysis_scope": "Read-only allele-orientation check for the seven mapped frailty factors, restricted to the pinned EUR HapMap3 GRCh37 map.",
        "script": str(Path(__file__).resolve().relative_to(repo)),
        "script_sha256": sha256(Path(__file__).resolve()),
        "reference_map": str(reference_map),
        "reference_map_sha256": sha256(reference_map),
        "reference_map_provenance": str(provenance_path),
        "reference_map_provenance_sha256": sha256(provenance_path),
        "reference_map_schema": ref_provenance["schema_version"],
        "reference_map_build": ref_provenance["genome_build"],
        "reference_map_rows": len(mapping),
        "factor_source_sha256": source_hashes,
        "acquired_resource_manifest": str(acquisition_manifest.relative_to(repo)),
        "acquired_resource_manifest_sha256": acquisition_manifest_sha256,
        "command": [
            sys.executable,
            str(Path(__file__).resolve()),
            "--repo",
            str(repo),
            "--reference-map",
            str(reference_map),
            "--output",
            str(output),
            "--provenance",
            str(provenance_out),
        ],
        "output": str(output.relative_to(repo)),
        "output_sha256": hashlib.sha256(table_text.encode()).hexdigest(),
        "python": sys.version.split()[0],
        "limitations": [
            "Only source rows whose rsID occurs in the HapMap3 map are assessed.",
            "A1/A2 are the bundled HapMap3/EUR LD reference allele pair, not a claim about genomic REF/ALT.",
            "Palindromic A/T and C/G variants remain orientation-ambiguous; effect-allele frequency is not used to infer strand.",
            "This audit neither changes effects nor produces harmonized output.",
        ],
    }
    provenance_text = json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    write_immutable(provenance_out, provenance_text)
    print(f"LATENT_HM3_ALLELE_AUDIT_OK traits={len(rows)} map_rows={len(mapping)} output={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
