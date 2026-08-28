#!/usr/bin/env python3
"""Apply pre-correlation eligibility and redundancy filters to the Pan-UKB universe."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Iterable


CORE_OR_SLEEP_PATTERNS = [
    r"\bsleep\b",
    r"\binsomnia\b",
    r"\bsnoring\b",
    r"\bchronotype\b",
    r"\basthma\b",
    r"\brheumatoid arthritis\b",
    r"\btype 2 diabetes\b",
    r"\bnon insulin dependent diabetes\b",
    r"\bdepressive episode\b",
    r"\bischemic heart disease\b",
    r"\bischaemic heart disease\b",
    r"\bcoronary atherosclerosis\b",
    r"\bldl direct\b",
]

# Pan-UKB contains both named-drug rows and duplicate pharmacologic-class rows.
# Only the interpretable named-drug representations enter the candidate pool.
NAMED_PRESCRIPTIONS = {
    "gabapentin",
    "paracetamol",
    "mebeverine",
    "amoxicillin_clavulanic_acid",
    "metformin",
    "amlodipine",
    "orlistat",
    "clarithromycin",
    "tramadol",
    "latanoprost",
    "citalopram",
    "gliclazide",
    "bendroflumethiazide",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def exclusion_reason(row: dict[str, str], max_source_intercept: float) -> str:
    try:
        intercept = float(row["LDSC_intercept"])
    except (KeyError, ValueError):
        return "missing_source_ldsc_intercept"
    if intercept > max_source_intercept:
        return "source_ldsc_intercept_above_1.2"

    name = " ".join([row["phenotype_name"], row["phenotype_definition"]]).lower()
    name = name.replace("-", " ")
    category = row["phenotype_category"].lower()
    if any(re.search(pattern, name) for pattern in CORE_OR_SLEEP_PATTERNS):
        return "core_or_sleep_overlap"
    administrative_names = [
        "methods of admission",
        "methods of discharge",
        "sources of admission",
        "destinations on discharge",
        "patient classification",
        "intended management",
    ]
    if (
        "summary administration" in category
        or "consultant" in name
        or any(value in row["phenotype_name"].lower() for value in administrative_names)
    ):
        return "administrative_process_phenotype"
    if "family history" in category or row["extension_trait_id"].startswith("panukbb_icd10__z82"):
        return "family_history_proxy"
    if "none of the above" in row["phenotype_definition"].lower():
        return "negative_composite_definition"
    if row["trait_type"] == "prescriptions":
        token = row["extension_trait_id"].split("__")[1]
        if token not in NAMED_PRESCRIPTIONS:
            return "duplicate_pharmacologic_class_row"
    return "PASS"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--annotated-universe",
        type=Path,
        default=Path("discovery_extension/results/panukbb_novelty_prescreen.tsv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("discovery_extension/results/candidate_pool.tsv"),
    )
    parser.add_argument(
        "--audit-out",
        type=Path,
        default=Path("discovery_extension/results/candidate_pool_exclusion_audit.tsv"),
    )
    parser.add_argument(
        "--run-metadata",
        type=Path,
        default=Path("discovery_extension/provenance/candidate_pool_build.json"),
    )
    parser.add_argument("--max-source-intercept", type=float, default=1.2)
    args = parser.parse_args()

    rows = read_tsv(args.annotated_universe)
    kept: list[dict[str, str]] = []
    audit: list[dict[str, str]] = []
    counts: Counter[str] = Counter()
    for row in rows:
        reason = exclusion_reason(row, args.max_source_intercept)
        counts[reason] += 1
        audit.append(
            {
                "extension_trait_id": row["extension_trait_id"],
                "phenotype_name": row["phenotype_name"],
                "candidate_pool_status": "INCLUDE" if reason == "PASS" else "EXCLUDE",
                "exclusion_reason": reason,
                "source_h2_z": row["h2_z"],
                "source_LDSC_intercept": row["LDSC_intercept"],
                "prior_sleep_screen_coverage": row["prior_sleep_screen_coverage"],
            }
        )
        if reason == "PASS":
            kept.append(row)

    if not 150 <= len(kept) <= 300:
        raise SystemExit(f"ERROR: candidate pool must contain 150-300 traits; found {len(kept)}")
    if len({row["extension_trait_id"] for row in kept}) != len(kept):
        raise SystemExit("ERROR: duplicate extension_trait_id in candidate pool")
    write_tsv(args.out, list(kept[0]), kept)
    write_tsv(args.audit_out, list(audit[0]), audit)

    metadata = {
        "schema_version": "1.0.0",
        "selection_timing": "before_any_extension_genetic_correlation",
        "source_rows": len(rows),
        "candidate_pool_rows": len(kept),
        "filters": {
            "source_PanUKBB_h2_z_minimum": 4.0,
            "source_LDSC_intercept_maximum": args.max_source_intercept,
            "sample_size_and_case_control_minima": "inherited_from_01_build_panukbb_universe.py",
            "exclusions": [
                "direct sleep phenotypes and close core phenotype representations",
                "hospital administrative/process phenotypes",
                "family-history proxies",
                "negative composite definitions",
                "duplicate pharmacologic-class rows when a named-drug row is available",
            ],
        },
        "status_counts": dict(sorted(counts.items())),
        "input_sha256": sha256(args.annotated_universe),
        "candidate_pool_sha256": sha256(args.out),
        "exclusion_audit_sha256": sha256(args.audit_out),
    }
    args.run_metadata.parent.mkdir(parents=True, exist_ok=True)
    args.run_metadata.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(f"CANDIDATE_POOL_OK rows={len(kept)} exclusions={len(rows) - len(kept)}")


if __name__ == "__main__":
    main()
