#!/usr/bin/env python3
"""Build an outcome-blind trait-level power audit from frozen LAVA v3 inputs."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "brain6/results/power_optimized_sensitivity_v1"
OUTPUT = OUTDIR / "canonical_family_trait_power_audit_v1.tsv"
PROVENANCE = OUTDIR / "canonical_family_trait_power_audit_v1.provenance.json"
MANIFEST = ROOT / "brain6/manifests/locked_dense_input_audit.tsv"
GWAS = ROOT / "brain6/manifests/gwas_master.tsv"
GWAS_QC = ROOT / "brain6/results/supplement/table_S1_gwas_metadata.tsv"
H2 = ROOT / "brain6/results/supplement/table_S3_h2_qc.tsv"
LAVA = ROOT / "work/lava-canonical-v3-production/d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/results/canonical_family_results.tsv"
STATUS = ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv"
CAUSES = ROOT / "brain6/results/lava/canonical_v3_not_run_cells.tsv"
REF_ROOT = ROOT / "ref/lava/ukb_v1.1"
REF_PROV = REF_ROOT / "reference.provenance.json"
EXPECTED_RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
REFERENCE_FIELDS = ("SNP", "CHR", "POS", "A1", "A2", "NOBS", "MISS", "FREQ", "NCORRS")
CAUSE_MAP = {
    "LOW_LOCAL_H2_UNDERPOWERED": "low_local_h2_not_run",
    "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS": "shared_reference_minK_not_run",
    "FEWER_THAN_MIN_K": "other_minK_not_run",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def reference_ids(path: Path, chromosome: int) -> set[str]:
    values: set[str] = set()
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REFERENCE_FIELDS:
            raise ValueError(f"Unexpected LAVA reference schema: {path}")
        for row in reader:
            if int(row["CHR"]) != chromosome:
                raise ValueError(f"Cross-chromosome row in {path}")
            value = row["SNP"].lower()
            if not value or value in values:
                raise ValueError(f"Empty or duplicate reference SNP ID in {path}")
            values.add(value)
    return values


def source_reference_overlap(source: Path, expected_bytes: int, expected_sha: str,
                             expected_reference_hashes: dict[str, str]) -> tuple[int, int, dict[str, str]]:
    if not source.is_file() or source.stat().st_size != expected_bytes:
        raise ValueError(f"Missing or size-mismatched locked input: {source}")
    observed_sha = sha256(source)
    if observed_sha != expected_sha:
        raise ValueError(f"Locked input checksum mismatch: {source}")
    total = overlap = 0
    source_ids: dict[int, set[str]] = defaultdict(set)
    ref_hashes: dict[str, str] = {}
    with gzip.open(source, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"SNP", "CHR"}.issubset(reader.fieldnames or ()):
            raise ValueError(f"Locked input lacks SNP/CHR columns: {source}")
        for row in reader:
            chrom = int(row["CHR"])
            snp = row["SNP"].lower()
            if not 1 <= chrom <= 22:
                raise ValueError(f"Invalid chromosome in {source}")
            if not snp or snp in source_ids[chrom]:
                raise ValueError(f"Empty or duplicate source SNP ID in {source}, chromosome {chrom}")
            source_ids[chrom].add(snp)
            total += 1
    # Locked dense inputs are not guaranteed to be chromosome sorted. Retain
    # one source trait's IDs in memory, then compare against one reference
    # chromosome at a time; this avoids sorting multi-million-row source files.
    for chrom, values in list(source_ids.items()):
        ref_file = REF_ROOT / f"lava-ukb-v1.1_chr{chrom}.info"
        if not ref_file.is_file():
            raise FileNotFoundError(ref_file)
        ref_hashes[str(chrom)] = sha256(ref_file)
        if ref_hashes[str(chrom)] != expected_reference_hashes.get(str(ref_file.relative_to(ROOT))):
            raise ValueError(f"Sealed reference checksum mismatch: {ref_file}")
        with ref_file.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if tuple(reader.fieldnames or ()) != REFERENCE_FIELDS:
                raise ValueError(f"Unexpected LAVA reference schema: {ref_file}")
            seen_reference: set[str] = set()
            for row in reader:
                if int(row["CHR"]) != chrom:
                    raise ValueError(f"Cross-chromosome reference row: {ref_file}")
                ref_id = row["SNP"].lower()
                if not ref_id or ref_id in seen_reference:
                    raise ValueError(f"Empty or duplicate reference SNP ID: {ref_file}")
                seen_reference.add(ref_id)
                overlap += ref_id in values
        del source_ids[chrom]
    return total, overlap, ref_hashes


def q(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def build() -> dict[str, Any]:
    if OUTPUT.exists() or PROVENANCE.exists():
        raise FileExistsError("Refusing to overwrite the immutable trait-power audit")
    ref_prov = json.loads(REF_PROV.read_text(encoding="utf-8"))
    if ref_prov.get("verification") != "SHA-256 verified after official HTTPS acquisition":
        raise ValueError("LAVA reference is not checksum sealed")
    manifest = {row["trait_id"]: row for row in read_tsv(MANIFEST)}
    expected_reference_hashes = {
        row["path"]: row["sha256"] for row in ref_prov.get("extracted_files", [])
        if row["path"].endswith(".info")
    }
    gwas = {row["trait"]: row for row in read_tsv(GWAS)}
    gwas_qc = {row["trait"]: row for row in read_tsv(GWAS_QC)}
    h2_rows = read_tsv(H2)
    h2_by_trait: dict[str, set[str]] = defaultdict(set)
    sleep_traits = {"insomnia", "longsleep"}
    for row in h2_rows:
        for trait, value in ((row["sleep_trait"], row["sleep_h2"]),
                             (row["brain_disorder"], row["h2_disorder"])):
            if value and value.lower() not in {"na", "nan", "not available"}:
                h2_by_trait[trait].add(value)
    h2_values = {}
    for trait, values in h2_by_trait.items():
        if len(values) != 1:
            raise ValueError(f"Expected one frozen SNP-h2 estimate for {trait}; found {values}")
        h2_values[trait] = next(iter(values))

    lava_rows = read_tsv(LAVA)
    if len(lava_rows) != 7 * 2495 or {row["phen"] for row in lava_rows} != set(manifest):
        raise ValueError("Canonical aggregate does not cover the seven frozen inputs")
    by_trait_lava: dict[str, list[dict[str, str]]] = defaultdict(list)
    local_h2: dict[str, list[float]] = defaultdict(list)
    for row in lava_rows:
        by_trait_lava[row["phen"]].append(row)
        if row["status"] == "TESTED" and row["h2.obs"] not in {"NA", ""}:
            local_h2[row["phen"]].append(float(row["h2.obs"]))
    status_rows = read_tsv(STATUS)
    status_by_trait = {row["trait_id"]: row for row in status_rows}
    cause_rows = read_tsv(CAUSES)
    cause_by_trait: dict[str, Counter[str]] = defaultdict(Counter)
    for row in cause_rows:
        cause_by_trait[row["trait_id"]][row["reason"]] += 1

    output_rows = []
    input_hashes: dict[str, str] = {}
    ref_hashes_all: dict[str, str] = {}
    for trait in sorted(manifest):
        locked = manifest[trait]
        if locked["expected_sha256"] != locked["observed_sha256"]:
            raise ValueError(f"Locked input audit has mismatched expected/observed hash: {trait}")
        source = Path(locked["dense_path"])
        total, overlap, ref_hashes = source_reference_overlap(
            source, int(locked["dense_bytes"]), locked["expected_sha256"], expected_reference_hashes)
        if total != int(gwas_qc[trait]["variant_count"]):
            raise ValueError(f"Decompressed row count differs from the frozen audit for {trait}")
        input_hashes[str(source)] = locked["expected_sha256"]
        ref_hashes_all.update({f"chr{chrom}": digest for chrom, digest in ref_hashes.items()})
        meta = gwas[trait]
        cases = int(meta["cases"]) if meta.get("cases", "").isdigit() else None
        controls = int(meta["controls"]) if meta.get("controls", "").isdigit() else None
        effective_n = 4 / (1 / cases + 1 / controls) if cases and controls else None
        aggregate = by_trait_lava[trait]
        counts = Counter(row["status"] for row in aggregate)
        ccounts = cause_by_trait[trait]
        loc_h2 = local_h2[trait]
        status = status_by_trait[trait]
        row = {
            "trait_id": trait,
            "role": "SLEEP_PHENOTYPE" if trait in sleep_traits else "BRAIN_DISORDER",
            "study": meta["study"], "PMID": meta["PMID"], "DOI": meta["DOI"],
            "source_sample_size": meta["sample_size"],
            "cases": meta["cases"], "controls": meta["controls"],
            "case_control_effective_N_approx": f"{effective_n:.1f}" if effective_n else "NA_NOT_CASE_CONTROL",
            "phenotype_definition": meta["cohorts"], "ancestry": meta["ancestry"],
            "genome_build": meta["genome_build"],
            "dense_unique_variants_after_frozen_QC": total,
            "exact_SNP_ID_reference_overlap": overlap,
            "reference_overlap_pct": f"{overlap / total * 100:.6f}",
            "SNP_h2_reported_in_locked_global_map": h2_values.get(trait, "NOT_RECORDED"),
            "planned_loci": int(status["planned_cells"]),
            "tested": int(status["tested_cells"]),
            "tested_pct": f"{int(status['tested_cells']) / int(status['planned_cells']) * 100:.4f}",
            "low_local_h2_NOT_RUN": ccounts["LOW_LOCAL_H2_UNDERPOWERED"],
            "shared_reference_minK_NOT_RUN": ccounts["FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS"],
            "other_minK_NOT_RUN": ccounts["FEWER_THAN_MIN_K"],
            "failed": counts["FAILED"], "NOT_RUN": counts["NOT_RUN"],
            "NOT_RUN_pct": f"{counts['NOT_RUN'] / len(aggregate) * 100:.4f}",
            "local_h2_median_tested": f"{statistics.median(loc_h2):.10g}" if loc_h2 else "NA",
            "local_h2_q25_tested": f"{q(loc_h2, .25):.10g}" if loc_h2 else "NA",
            "local_h2_q75_tested": f"{q(loc_h2, .75):.10g}" if loc_h2 else "NA",
            "local_h2_summary_scope": "TESTED cells only; NOT_RUN cells have no local-h2 estimate",
            "immutable_run_id": EXPECTED_RUN_ID,
        }
        if row["low_local_h2_NOT_RUN"] + row["shared_reference_minK_NOT_RUN"] + row["other_minK_NOT_RUN"] != row["NOT_RUN"]:
            raise ValueError(f"NOT_RUN causes do not sum to NOT_RUN for {trait}")
        if row["tested"] + row["NOT_RUN"] + row["failed"] != row["planned_loci"]:
            raise ValueError(f"Cell statuses do not sum to planned loci for {trait}")
        output_rows.append(row)

    fields = list(output_rows[0])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    inputs = [MANIFEST, GWAS, GWAS_QC, H2, LAVA, STATUS, CAUSES, REF_PROV]
    for path in inputs:
        input_hashes[str(path.relative_to(ROOT))] = sha256(path)
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_canonical_family_trait_power_audit_v1",
        "status": "PASS_HASH_VERIFIED_OUTCOME_BLIND_TRAIT_POWER_AUDIT",
        "scope": "All seven immutable canonical LAVA inputs; two sleep phenotypes and five brain-disorder GWAS, not seven sleep traits.",
        "decision_boundary": "Technical diagnostics only. No bivariate rg, local-rg p-values, or downstream sleep-disorder association results were used to rank traits.",
        "reference_overlap_method": "Exact lower-case SNP ID intersection, chromosome matched, using the sealed LAVA UKB v1.1 info panel; this does not establish allele/effect harmonization.",
        "local_h2_method": "Summarize canonical aggregate h2.obs across TESTED cells only; NOT_RUN cells have no reported estimate and are not imputed.",
        "canonical_run_id": EXPECTED_RUN_ID,
        "inputs_sha256": input_hashes,
        "reference_chromosome_sha256": ref_hashes_all,
        "output": {"path": str(OUTPUT.relative_to(ROOT)), "sha256": sha256(OUTPUT), "rows": len(output_rows)},
        "builder_path": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
        "rank_by_processability": sorted(
            ((row["trait_id"], float(row["tested_pct"])) for row in output_rows), key=lambda item: item[1]),
        "limitations": [
            "The canonical family contains two sleep GWAS (insomnia and long sleep) plus five brain-disorder GWAS; no other sleep traits are part of this family.",
            "Genome-wide SNP ID overlap is not allele/effect harmonization and does not replace per-locus shared-K receipt checks.",
            "SNP heritability values are inherited from locked global-map LDSC results and retain their original scale/QC limitations.",
            "Case-control effective N is the conventional 4/(1/cases+1/controls) approximation, not a replacement for source-specific power calculations.",
            "The ranking is technical processability under this frozen reference and gate, not a general statement about biological power.",
        ],
    }
    with PROVENANCE.open("x", encoding="utf-8") as stream:
        json.dump(provenance, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return {"status": provenance["status"], "rows": len(output_rows), "rank_by_processability": provenance["rank_by_processability"], "output_sha256": sha256(OUTPUT)}


def validate() -> dict[str, Any]:
    p = json.loads(PROVENANCE.read_text(encoding="utf-8"))
    if p.get("status") != "PASS_HASH_VERIFIED_OUTCOME_BLIND_TRAIT_POWER_AUDIT":
        raise ValueError("Audit status is invalid")
    if sha256(OUTPUT) != p.get("output", {}).get("sha256"):
        raise ValueError("Audit output hash mismatch")
    for name, digest in p.get("inputs_sha256", {}).items():
        path = Path(name) if name.startswith("/") else ROOT / name
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Input hash mismatch: {name}")
    if p.get("builder_sha256") != sha256(Path(__file__).resolve()):
        raise ValueError("Audit builder hash mismatch")
    for name, digest in p.get("reference_chromosome_sha256", {}).items():
        path = REF_ROOT / f"lava-ukb-v1.1_{name}.info"
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"Sealed reference hash mismatch: {path}")
    rows = read_tsv(OUTPUT)
    if len(rows) != 7 or sum(row["role"] == "SLEEP_PHENOTYPE" for row in rows) != 2:
        raise ValueError("Expected seven canonical inputs, including exactly two sleep phenotypes")
    for row in rows:
        if int(row["tested"]) + int(row["NOT_RUN"]) + int(row["failed"]) != 2495:
            raise ValueError(f"Incomplete status accounting: {row['trait_id']}")
    if sorted(((row["trait_id"], float(row["tested_pct"])) for row in rows),
              key=lambda item: item[1]) != [tuple(item) for item in p["rank_by_processability"]]:
        raise ValueError("Processability ranking differs from the published audit table")
    return {"status": "PASS", "rows": len(rows), "output_sha256": sha256(OUTPUT), "rank_by_processability": p["rank_by_processability"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--validate", action="store_true")
    args = parser.parse_args()
    print(json.dumps(build() if args.write else validate(), indent=2, sort_keys=True))
