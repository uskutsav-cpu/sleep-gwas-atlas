#!/usr/bin/env python3
"""Build locus-scoped Brain6 LAVA inputs from sealed normalized sources.

Each normalized source is read through its receipt-bound SQLite index, and
each official chromosome reference index is read once. Each pair-by-locus job
then opens only the variants it needs instead of decompressing a whole
chromosome repeatedly across the 12,475 frozen tests.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import os
import sqlite3
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions" / "brain6"))
from brain6.artifacts import verify_artifact
from brain6.io import read_json, sha256

TRAITS = ("adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz")
FIELDS = ("SNP", "A1", "A2", "Z", "N")
INFO_FIELDS = ("SNP", "CHR", "POS", "A1", "A2", "NOBS", "MISS", "FREQ", "NCORRS")


def sealed_reference(reference: Path) -> dict:
    lock_path = reference / "reference.provenance.json"
    lock = read_json(lock_path)
    if lock.get("schema_version") != "sleep-atlas-lava-reference.1":
        raise ValueError("Unsupported or missing sealed LAVA reference provenance")
    if lock.get("verification") != "SHA-256 verified after official HTTPS acquisition":
        raise ValueError("Official reference provenance is not in a verified state")
    records = {r["path"]: r for r in lock.get("extracted_files", [])}
    for chromosome in range(1, 23):
        relative = f"ref/lava/ukb_v1.1/lava-ukb-v1.1_chr{chromosome}.info"
        record = records.get(relative)
        path = reference / f"lava-ukb-v1.1_chr{chromosome}.info"
        if not record or not path.is_file():
            raise ValueError(f"Reference provenance omits chromosome {chromosome} info file")
        if path.stat().st_size != int(record["bytes"]) or sha256(path) != record["sha256"]:
            raise ValueError(f"Reference info file differs from sealed provenance: chr{chromosome}")
    return lock


def load_case_control(path: Path) -> dict[str, tuple[int, int]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {}
    for row in rows:
        trait = row["trait"]
        if trait in TRAITS:
            cases, controls = int(row["cases"]), int(row["controls"])
            if cases <= 0 or controls <= 0:
                raise ValueError(f"Invalid locked case/control counts for {trait}")
            result[trait] = cases, controls
    if set(result) != set(TRAITS):
        raise ValueError("The case/control ledger does not cover the exact seven LAVA traits")
    return result


def load_reference_ids(path: Path, chromosome: int) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != INFO_FIELDS:
            raise ValueError(f"Unexpected LAVA reference index schema: chr{chromosome}")
        ids = set()
        for row in reader:
            if int(row["CHR"]) != chromosome:
                raise ValueError(f"Cross-chromosome row in LAVA chr{chromosome} index")
            snp = row["SNP"].lower()
            if not snp or snp in ids:
                raise ValueError(f"Empty or duplicate reference SNP in chromosome {chromosome}: {snp}")
            ids.add(snp)
    if not ids:
        raise ValueError(f"Empty LAVA reference index for chromosome {chromosome}")
    return ids


def open_gzip_text(path: Path):
    raw = path.open("xb")
    compressed = gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=raw, mtime=0)
    return raw, io.TextIOWrapper(compressed, encoding="utf-8", newline="")


def load_loci(path: Path, family: dict) -> list[dict[str, str]]:
    if sha256(path) != family["locus_definition"]["sha256"]:
        raise ValueError("LAVA locus-definition checksum differs from the frozen family")
    with path.open(encoding="utf-8") as handle:
        lines = [line.split() for line in handle if line.strip()]
    if not lines or lines[0] != ["LOC", "CHR", "START", "STOP"]:
        raise ValueError("LAVA locus definition header drifted")
    if any(len(row) != len(lines[0]) for row in lines[1:]):
        raise ValueError("Malformed LAVA locus definition row")
    loci = [dict(zip(lines[0], row)) for row in lines[1:]]
    if len(loci) != int(family["locus_definition"]["n_loci"]):
        raise ValueError("LAVA locus count differs from the frozen family")
    if any(not {"LOC", "CHR", "START", "STOP"} <= set(row) for row in loci):
        raise ValueError("LAVA locus definition lacks required coordinates")
    if len({row["LOC"] for row in loci}) != len(loci):
        raise ValueError("LAVA locus identifiers are duplicated")
    return loci


def build(normalized_root: Path, reference: Path, output: Path, family_lock: Path) -> Path:
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to replace existing LAVA inputs: {output}")
    lock = read_json(family_lock)
    if lock.get("analysis_id") != "brain6-lava-local-rg-v2" or lock.get("scientific_status") != "NOT_RUN":
        raise ValueError("Active Brain6 LAVA v2 family lock is missing or has already been run")
    reference_lock = sealed_reference(reference)
    cases = load_case_control(ROOT / "brain6/manifests/gwas_master.tsv")
    loci_path = ROOT / lock["locus_definition"]["path"]
    loci = load_loci(loci_path, lock)
    pair_lock_path = ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json"
    pair_lock = read_json(pair_lock_path)
    pair_body = {key: value for key, value in pair_lock.items() if key != "lock_sha256"}
    pair_digest = hashlib.sha256(json.dumps(pair_body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if pair_digest != lock["pair_lock_sha256"]:
        raise ValueError("Reviewed pair selection changed after the LAVA family was frozen")
    pair_rows = {row["pair_id"]: row for row in pair_lock["pairs"] if row["pair_id"] in lock["pairs"]}
    if set(pair_rows) != set(lock["pairs"]):
        raise ValueError("Reviewed pair lock does not cover all LAVA pairs")

    sources = {}
    for trait in TRAITS:
        artifact = normalized_root / f"normalize_{trait}"
        receipt = verify_artifact(artifact)
        if receipt.get("scientific_status") != "PASS" or receipt.get("synthetic"):
            raise ValueError(f"Rejected or synthetic normalized source: {trait}")
        sources[trait] = {
            "artifact": artifact,
            "receipt": receipt,
            "db": sqlite3.connect(f"file:{(artifact / 'variants.sqlite').resolve()}?mode=ro", uri=True),
        }

    temp_parent = output.parent
    temp_parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{output.name}.", suffix=".partial", dir=temp_parent))
    try:
        records = []
        coverage = []
        for chromosome in range(1, 23):
            ref_path = reference / f"lava-ukb-v1.1_chr{chromosome}.info"
            reference_ids = load_reference_ids(ref_path, chromosome)
            for locus in (value for value in loci if int(value["CHR"]) == chromosome):
                locus_id = str(locus["LOC"])
                locus_dir = stage / f"locus_{locus_id}"
                locus_dir.mkdir()
                ids_by_trait: dict[str, set[str]] = {}
                start, stop = int(locus["START"]), int(locus["STOP"])
                for trait in TRAITS:
                    stat_path = locus_dir / f"{trait}.sumstats.tsv.gz"
                    raw, text = open_gzip_text(stat_path)
                    writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                    writer.writerow(FIELDS)
                    ids = set()
                    cursor = sources[trait]["db"].execute(
                        "SELECT snp,a1,a2,beta,se,n FROM variants WHERE chr=? AND bp>=? AND bp<=? ORDER BY bp,snp",
                        (chromosome, start, stop),
                    )
                    count = 0
                    try:
                        for snp, a1, a2, beta, se, n in cursor:
                            normalized_snp = str(snp).lower()
                            if normalized_snp not in reference_ids:
                                continue
                            beta, se, n = float(beta), float(se), float(n)
                            if not (se > 0 and n > 0 and all(map(math.isfinite, (beta, se, n)))):
                                raise ValueError(f"Invalid normalized association statistic for {trait}:{snp}")
                            writer.writerow((normalized_snp, a1, a2, format(beta / se, ".17g"), format(n, ".17g")))
                            ids.add(normalized_snp)
                            count += 1
                    finally:
                        text.close()
                        raw.close()
                    ids_by_trait[trait] = ids
                    final_stat_path = output / stat_path.relative_to(stage)
                    records.append({
                        "kind": "sumstats", "locus_id": locus_id, "chromosome": chromosome,
                        "trait": trait, "path": str(final_stat_path), "bytes": stat_path.stat().st_size,
                        "sha256": sha256(stat_path), "rows": count,
                        "source_receipt_fingerprint": sources[trait]["receipt"]["fingerprint"],
                    })
                input_info = locus_dir / "input_info.tsv"
                with input_info.open("x", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(("phenotype", "cases", "controls", "filename"))
                    for trait in TRAITS:
                        writer.writerow((trait, *cases[trait], f"{trait}.sumstats.tsv.gz"))
                records.append({"kind": "input_info", "locus_id": locus_id,
                                "path": str(output / input_info.relative_to(stage)), "bytes": input_info.stat().st_size,
                                "sha256": sha256(input_info)})
                for pair_id in lock["pairs"]:
                    pair = pair_rows[pair_id]
                    trait1, trait2 = pair["sleep_trait"], pair["disease_trait"]
                    count1, count2 = len(ids_by_trait[trait1]), len(ids_by_trait[trait2])
                    shared = len(ids_by_trait[trait1] & ids_by_trait[trait2])
                    coverage.append({"pair_id": pair_id, "locus_id": locus_id,
                        "trait1_variants": count1, "trait2_variants": count2,
                        "shared_reference_variants": shared,
                        "status": "READY" if shared >= 2 else "NO_OVERLAP_LT_MIN_K"})
                for ids in ids_by_trait.values():
                    ids.clear()
            reference_ids.clear()

        coverage_path = stage / "locus_pair_coverage.tsv"
        with coverage_path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=("pair_id", "locus_id", "trait1_variants",
                "trait2_variants", "shared_reference_variants", "status"),
                delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(coverage)

        manifest = {
            "schema_version": "brain6-lava-locus-inputs.1",
            "analysis_id": lock["analysis_id"],
            "materializer_script_sha256": sha256(Path(__file__)),
            "family_lock_sha256": sha256(family_lock),
            "case_control_ledger_sha256": sha256(ROOT / "brain6/manifests/gwas_master.tsv"),
            "official_reference_provenance_sha256": sha256(reference / "reference.provenance.json"),
            "official_reference_archive_bytes": reference_lock["archive_bytes"],
            "normalized_sources": [
                {"trait": trait, "receipt_fingerprint": sources[trait]["receipt"]["fingerprint"],
                 "receipt_sha256": sha256(sources[trait]["artifact"] / "receipt.json")}
                for trait in TRAITS
            ],
            "n_chromosomes": 22,
            "n_traits": len(TRAITS),
            "n_loci": len(loci),
            "pair_coverage": {"path": str(output / coverage_path.relative_to(stage)), "bytes": coverage_path.stat().st_size,
                              "sha256": sha256(coverage_path), "rows": len(coverage)},
            "records": records,
        }
        with (stage / "provenance.json").open("x", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.rename(stage, output)
    finally:
        for source in sources.values():
            source["db"].close()
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--normalized-root", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=ROOT / "ref/lava/ukb_v1.1")
    parser.add_argument("--family-lock", type=Path, default=ROOT / "brain6/config/lava_family_v2.json")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    path = build(args.normalized_root.resolve(), args.reference.resolve(), args.out.resolve(), args.family_lock.resolve())
    print(f"BRAIN6_LAVA_INPUTS_READY path={path}")


if __name__ == "__main__":
    main()
