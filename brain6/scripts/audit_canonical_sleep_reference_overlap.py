#!/usr/bin/env python3
"""Audit canonical sleep GWAS variant-ID overlap with sealed LAVA reference.

This is a read-only Brain6 diagnostic. It verifies each harmonized dense source
against ``locked_dense_input_audit.tsv`` and counts exact, chromosome-matched
SNP-ID intersections with the sealed UKB LAVA v1.1 reference files.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/power_optimized_sensitivity_v1"
TRAITS = ("insomnia", "longsleep")
REFERENCE_FIELDS = ("SNP", "CHR", "POS", "A1", "A2", "NOBS", "MISS", "FREQ", "NCORRS")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def reference_ids(path: Path, chromosome: int) -> set[str]:
    ids: set[str] = set()
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REFERENCE_FIELDS:
            raise ValueError(f"unexpected reference schema: {path}")
        for row in reader:
            if int(row["CHR"]) != chromosome:
                raise ValueError(f"cross-chromosome reference row: {path}")
            snp = row["SNP"].lower()
            if not snp or snp in ids:
                raise ValueError(f"empty or duplicate reference SNP ID: {path}")
            ids.add(snp)
    return ids


def audit_trait(trait: str, manifest_row: dict[str, str], reference_root: Path) -> dict[str, Any]:
    source = Path(manifest_row["dense_path"])
    if not source.is_file() or source.stat().st_size != int(manifest_row["dense_bytes"]):
        raise ValueError(f"locked dense source missing or size mismatch for {trait}: {source}")
    source_hash = sha256(source)
    if source_hash != manifest_row["expected_sha256"] or source_hash != manifest_row["observed_sha256"]:
        raise ValueError(f"locked dense source checksum mismatch for {trait}: {source}")

    total_rows = 0
    overlap_rows = 0
    last_chromosome = 0
    seen_on_chromosome: set[str] = set()
    reference_snp_ids: set[str] = set()
    reference_hashes: dict[str, str] = {}
    with gzip.open(source, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"SNP", "CHR"}.issubset(reader.fieldnames or ()):
            raise ValueError(f"dense source lacks SNP/CHR fields: {source}")
        for row in reader:
            chromosome = int(row["CHR"])
            snp = row["SNP"].lower()
            if not 1 <= chromosome <= 22 or chromosome < last_chromosome:
                raise ValueError(f"source rows are not ordered by chromosome: {source}")
            if chromosome != last_chromosome:
                reference_snp_ids.clear()
                seen_on_chromosome.clear()
                reference_file = reference_root / f"lava-ukb-v1.1_chr{chromosome}.info"
                reference_hashes[str(chromosome)] = sha256(reference_file)
                reference_snp_ids = reference_ids(reference_file, chromosome)
                last_chromosome = chromosome
            if not snp or snp in seen_on_chromosome:
                raise ValueError(f"empty or duplicate source SNP ID in {trait}, chromosome {chromosome}")
            seen_on_chromosome.add(snp)
            total_rows += 1
            overlap_rows += snp in reference_snp_ids

    return {
        "trait": trait,
        "source_path": str(source),
        "source_bytes": source.stat().st_size,
        "source_sha256": source_hash,
        "source_unique_variant_rows": total_rows,
        "reference_overlap_unique_variant_ids": overlap_rows,
        "reference_nonoverlap_variant_ids": total_rows - overlap_rows,
        "reference_overlap_fraction": overlap_rows / total_rows,
        "reference_overlap_percent": 100 * overlap_rows / total_rows,
        "reference_chromosome_sha256": reference_hashes,
    }


def build(output_dir: Path = BASE) -> dict[str, Any]:
    manifest_path = ROOT / "brain6/manifests/locked_dense_input_audit.tsv"
    reference_root = ROOT / "ref/lava/ukb_v1.1"
    reference_provenance_path = reference_root / "reference.provenance.json"
    base_table = BASE / "current_sleep_trait_audit.tsv"
    overlap_path = output_dir / "canonical_sleep_reference_overlap_v1.tsv"
    complete_table_path = output_dir / "current_sleep_trait_audit_complete.tsv"
    provenance_path = output_dir / "canonical_sleep_reference_overlap_v1.provenance.json"
    outputs = (overlap_path, complete_table_path, provenance_path)
    exists = [path for path in outputs if path.exists()]
    if exists:
        raise FileExistsError(f"Refusing to overwrite overlap audit outputs: {exists}")

    reference_provenance = json.loads(reference_provenance_path.read_text(encoding="utf-8"))
    if reference_provenance.get("verification") != "SHA-256 verified after official HTTPS acquisition":
        raise ValueError("LAVA reference provenance is not sealed as verified")
    with manifest_path.open(encoding="utf-8", newline="") as stream:
        manifest = {row["trait_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
    with base_table.open(encoding="utf-8", newline="") as stream:
        base_rows = list(csv.DictReader(stream, delimiter="\t"))
    if {row["trait"] for row in base_rows} != set(TRAITS):
        raise ValueError("current Brain6 sleep audit must contain insomnia and longsleep only")

    audits = [audit_trait(trait, manifest[trait], reference_root) for trait in TRAITS]
    by_trait = {row["trait"]: row for row in audits}
    for row in base_rows:
        expected = int(row["usable_variant_count"])
        if by_trait[row["trait"]]["source_unique_variant_rows"] != expected:
            raise ValueError(f"source variant count differs from locked sleep audit for {row['trait']}")
        row["frozen_reference_overlap_count"] = str(by_trait[row["trait"]]["reference_overlap_unique_variant_ids"])
        row["frozen_reference_nonoverlap_count"] = str(by_trait[row["trait"]]["reference_nonoverlap_variant_ids"])
        row["frozen_reference_overlap_pct"] = f"{by_trait[row['trait']]['reference_overlap_percent']:.6f}"
        row["reference_overlap_method"] = "Exact lower-case SNP ID intersection, chromosome matched; not allele harmonization"

    overlap_fields = (
        "trait", "source_path", "source_bytes", "source_sha256", "source_unique_variant_rows",
        "reference_overlap_unique_variant_ids", "reference_nonoverlap_variant_ids",
        "reference_overlap_fraction", "reference_overlap_percent",
    )
    overlap_path.parent.mkdir(parents=True, exist_ok=True)
    with overlap_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=overlap_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows([{key: row[key] for key in overlap_fields} for row in audits])
    complete_fields = list(base_rows[0])
    with complete_table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=complete_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(base_rows)

    inputs = {
        str(manifest_path.relative_to(ROOT)): sha256(manifest_path),
        str(base_table.relative_to(ROOT)): sha256(base_table),
        str(reference_provenance_path.relative_to(ROOT)): sha256(reference_provenance_path),
    }
    for audit in audits:
        source = Path(audit["source_path"])
        inputs[audit["source_path"]] = audit["source_sha256"]
        for chromosome, digest in audit["reference_chromosome_sha256"].items():
            path = reference_root / f"lava-ukb-v1.1_chr{chromosome}.info"
            inputs[str(path.relative_to(ROOT))] = digest
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_canonical_sleep_reference_overlap_v1",
        "scope": "Brain6 canonical sleep traits only; exact SNP-ID overlap with the frozen LAVA UKB v1.1 panel",
        "status": "PASS_HASH_VERIFIED_EXACT_ID_OVERLAP",
        "method": "Stream chromosome-ordered locked harmonized sources; intersect lower-case SNP IDs with the matching sealed LAVA chromosome info file; validate no duplicate IDs within chromosome and verify all locked source/reference hashes.",
        "limitations": [
            "ID overlap does not establish allele, effect, or frequency harmonization.",
            "Aggregate genome-wide overlap does not replace per-locus LAVA shared-reference/min-K eligibility, which remains sourced from canonical v3 receipts and Tables S25/S26.",
            "Brain6 v3 contains two sleep GWAS (insomnia and longsleep) among seven total inputs; the other five inputs are disorders.",
        ],
        "reference_provenance_sha256": sha256(reference_provenance_path),
        "inputs_sha256": inputs,
        "trait_results": audits,
        "outputs_sha256": {
            "canonical_sleep_reference_overlap_v1.tsv": sha256(overlap_path),
            "current_sleep_trait_audit_complete.tsv": sha256(complete_table_path),
        },
        "builder_path": str(Path(__file__).resolve().relative_to(ROOT)),
        "builder_sha256": sha256(Path(__file__).resolve()),
    }
    with provenance_path.open("x", encoding="utf-8") as stream:
        json.dump(provenance, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return {
        "status": provenance["status"],
        "traits": [{key: row[key] for key in (
            "trait", "source_unique_variant_rows", "reference_overlap_unique_variant_ids",
            "reference_overlap_percent")} for row in audits],
        "output_sha256": {path.name: sha256(path) for path in outputs},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=BASE)
    args = parser.parse_args()
    result = build(args.output_dir if args.output_dir.is_absolute() else ROOT / args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
