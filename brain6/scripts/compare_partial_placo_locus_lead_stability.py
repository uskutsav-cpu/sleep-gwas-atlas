#!/usr/bin/env python3
"""Compare partial PLACO lead sets across the existing LD diagnostics.

This is a diagnostic comparison only. It does not promote candidate loci,
claim genotype-validated LD, or assign cross-trait locus tiers.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6/results/loci"
METHODS = {
    "baseline_partial": OUT / "placo_candidate_loci_partial.tsv",
    "range_edges_excluded": OUT / "placo_candidate_loci_valid_ld_partial_v2.tsv",
    "factor_normalized_sensitivity": OUT / "placo_factor_normalized_sensitivity_v2/placo_candidate_loci.tsv",
}
VARIANTS = {
    "baseline_partial": OUT / "placo_candidate_variants_partial.tsv",
    "range_edges_excluded": OUT / "placo_candidate_variants_valid_ld_partial_v2.tsv",
    "factor_normalized_sensitivity": OUT / "placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv",
}
RECEIPTS = {
    "baseline_partial": OUT / "placo_candidate_loci_partial.provenance.json",
    "range_edges_excluded": OUT / "placo_candidate_valid_ld_partial_v2.provenance.json",
    "factor_normalized_sensitivity": OUT / "placo_factor_normalized_sensitivity_v2/provenance.json",
}
RECEIPT_STATUSES = {
    "baseline_partial": "PASS_PARTIAL_FAMILY",
    "range_edges_excluded": "PASS_PARTIAL_WITH_LD_RANGE_EXCEPTIONS",
    "factor_normalized_sensitivity": "DIAGNOSTIC_SENSITIVITY_NOT_PROMOTED",
}
OUT_PATH = OUT / "placo_candidate_locus_lead_stability_v1.tsv"
PROVENANCE_PATH = OUT / "placo_candidate_locus_lead_stability_v1.provenance.json"
PAIR_SET = {"insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz"}
FIELDS = [
    "pair_id", "lead_variant", "CHR", "BP", "P_PLACO",
    "baseline_partial", "range_edges_excluded", "factor_normalized_sensitivity",
    "stability_class", "interpretation",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate_output_receipt(receipt_path: Path, output_paths: tuple[Path, Path], expected_status: str) -> None:
    """Require each candidate table to match its producing stage receipt."""
    document = json.loads(receipt_path.read_text(encoding="utf-8"))
    if document.get("status") != expected_status:
        raise ValueError(f"unexpected source-stage status in {receipt_path}: {document.get('status')}")
    records = document.get("outputs")
    if records is None:
        records = document.get("outputs", [])
    recorded = {}
    for record in records:
        path = record.get("path")
        digest = record.get("sha256")
        if path and digest:
            recorded[Path(path).as_posix()] = digest
    for path in output_paths:
        relpath = path.relative_to(ROOT).as_posix()
        if recorded.get(relpath) != sha256(path):
            raise ValueError(f"candidate source receipt does not bind the expected output: {relpath}")


def extract_leads(locus_rows: list[dict[str, str]], variant_rows: list[dict[str, str]], method: str
                 ) -> dict[tuple[str, str], dict[str, str]]:
    positions: dict[tuple[str, str], dict[str, str]] = {}
    for row in variant_rows:
        key = row["pair_id"], row["SNP"]
        if row.get("candidate_status") == "LEAD":
            if key in positions:
                raise ValueError(f"duplicate lead metadata for {method}: {key}")
            positions[key] = row

    leads: dict[tuple[str, str], dict[str, str]] = {}
    for locus in locus_rows:
        pair = locus["pair_id"]
        for lead in filter(None, locus["lead_variants"].split(";")):
            key = pair, lead
            variant = positions.get(key)
            if variant is None:
                raise ValueError(f"missing lead metadata for {method}: {key}")
            if key in leads:
                raise ValueError(f"lead is assigned to more than one interval for {method}: {key}")
            if (variant["CHR"] != locus["CHR"] or
                    not int(locus["START"]) - 500_000 <= int(variant["BP"]) <= int(locus["STOP"]) + 500_000):
                raise ValueError(f"lead coordinate does not agree with its locus row for {method}: {key}")
            leads[key] = variant
    if len(locus_rows) != 19 or {pair for pair, _ in leads} != PAIR_SET:
        raise ValueError(f"unexpected pair coverage for {method}: {sorted({p for p, _ in leads})}")
    return leads


def build() -> dict[str, Any]:
    leads_by_method: dict[str, dict[tuple[str, str], dict[str, str]]] = {}
    source_hashes: dict[str, str] = {}
    for method in METHODS:
        receipt = RECEIPTS[method]
        if not receipt.is_file():
            raise FileNotFoundError(receipt)
        validate_output_receipt(receipt, (METHODS[method], VARIANTS[method]), RECEIPT_STATUSES[method])
        source_hashes[receipt.relative_to(ROOT).as_posix()] = sha256(receipt)
        for path in (METHODS[method], VARIANTS[method]):
            if not path.is_file():
                raise FileNotFoundError(path)
            source_hashes[path.relative_to(ROOT).as_posix()] = sha256(path)
        lead_map = extract_leads(read_tsv(METHODS[method]), read_tsv(VARIANTS[method]), method)
        leads_by_method[method] = lead_map

    baseline = set(leads_by_method["baseline_partial"])
    range_filtered = set(leads_by_method["range_edges_excluded"])
    normalized = set(leads_by_method["factor_normalized_sensitivity"])
    if baseline != normalized:
        raise ValueError("baseline and factor-normalized sensitivity lead sets differ")
    if not baseline <= range_filtered:
        raise ValueError("range-edge-excluded lead set is not a superset of the baseline")
    if len(baseline) != 21 or len(range_filtered) != 25 or len(range_filtered - baseline) != 4:
        raise ValueError("observed lead-set sizes differ from the verified diagnostic comparison")

    rows: list[dict[str, str]] = []
    all_keys = sorted(baseline | range_filtered, key=lambda x: (x[0], x[1]))
    for key in all_keys:
        pair, lead = key
        flags = {method: key in leads_by_method[method] for method in METHODS}
        classification = (
            "STABLE_ALL_THREE" if all(flags.values()) else
            "BASELINE_AND_FACTOR_ONLY" if flags["baseline_partial"] and flags["factor_normalized_sensitivity"] else
            "RANGE_FILTER_ONLY" if flags["range_edges_excluded"] else
            "UNEXPECTED_METHOD_PATTERN"
        )
        metadata = next(leads_by_method[method][key] for method in METHODS if key in leads_by_method[method])
        comparable_metadata = [leads_by_method[method][key] for method in METHODS if key in leads_by_method[method]]
        if len({(r["CHR"], r["BP"], r["P_PLACO"]) for r in comparable_metadata}) != 1:
            raise ValueError(f"candidate coordinates or PLACO statistic differ across methods: {key}")
        rows.append({
            "pair_id": pair,
            "lead_variant": lead,
            "CHR": metadata["CHR"],
            "BP": metadata["BP"],
            "P_PLACO": metadata["P_PLACO"],
            "baseline_partial": str(flags["baseline_partial"]).upper(),
            "range_edges_excluded": str(flags["range_edges_excluded"]).upper(),
            "factor_normalized_sensitivity": str(flags["factor_normalized_sensitivity"]).upper(),
            "stability_class": classification,
            "interpretation": "Diagnostic PLACO candidate lead only; incomplete five-track family and no genotype-level LD validation",
        })
    if any(row["stability_class"] == "UNEXPECTED_METHOD_PATTERN" for row in rows):
        raise ValueError("unexpected lead-set method pattern")

    with OUT_PATH.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_placo_candidate_locus_lead_stability_v1",
        "status": "PASS_DIAGNOSTIC_ONLY_NO_PROMOTION",
        "scope": {
            "pairs": sorted(PAIR_SET), "baseline_leads": len(baseline),
            "range_edges_excluded_leads": len(range_filtered),
            "factor_normalized_sensitivity_leads": len(normalized),
            "stable_all_three": sum(r["stability_class"] == "STABLE_ALL_THREE" for r in rows),
            "baseline_and_factor_only": sum(r["stability_class"] == "BASELINE_AND_FACTOR_ONLY" for r in rows),
            "range_filter_only": sum(r["stability_class"] == "RANGE_FILTER_ONLY" for r in rows),
            "candidate_intervals_each_method": 19,
            "full_five_track_family": False,
            "genotype_level_ld_validation": False,
            "fine_mapping_or_colocalization": False,
            "tier_promotion": False,
        },
        "method": {
            "candidate_threshold": 1e-8,
            "r2_threshold": 0.1,
            "window_kb": 500,
            "comparison": "Lead-ID set comparison across baseline partial clumping, raw out-of-range-edge exclusion, and factor-normalized internal sensitivity",
            "interpretation_boundary": "The four range-filter-only leads are method-sensitive candidates, not confirmed independent loci.",
        },
        "sources": source_hashes,
        "builder": {"path": Path(__file__).resolve().relative_to(ROOT).as_posix(),
                    "sha256": sha256(Path(__file__).resolve())},
        "output": {"path": OUT_PATH.relative_to(ROOT).as_posix(), "sha256": sha256(OUT_PATH), "rows": len(rows)},
    }
    PROVENANCE_PATH.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return provenance


if __name__ == "__main__":
    result = build()
    print(json.dumps(result, indent=2, sort_keys=True))
