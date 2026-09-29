#!/usr/bin/env python3
"""Audit whether sleep-trait replacements can satisfy the frozen family missingness gate."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
CANONICAL_RUN = ROOT / "work/lava-canonical-v3-production" / RUN_ID
CANONICAL_DECISION = CANONICAL_RUN / "canonical_family_decision.json"
CANONICAL_BY_TRAIT = ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv"
CANDIDATE_SUMMARY = ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/summary.json"
PROMOTION_REASSESSMENT = ROOT / "brain6/results/power_optimized_sensitivity_v1/promotion_rule_reassessment_2026-09-26.json"
DEFAULT_OUTPUT = ROOT / "brain6/results/power_optimized_sensitivity_v1/sensitivity_family_feasibility_v1.json"
EXPECTED_CELLS = 17_465
FAMILY_LIMIT = 873


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def audit(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    decision = json.loads(CANONICAL_DECISION.read_text(encoding="utf-8"))
    candidate = json.loads(CANDIDATE_SUMMARY.read_text(encoding="utf-8"))
    promotion = json.loads(PROMOTION_REASSESSMENT.read_text(encoding="utf-8"))
    trait_rows = read_tsv(CANONICAL_BY_TRAIT)
    not_run_by_trait = {row["trait_id"]: int(row["not_run_cells"]) for row in trait_rows}
    if (decision.get("run_id") != RUN_ID or decision.get("overall_status") != "FAILED_QC_NOT_PROMOTED" or
            decision.get("status_counts", {}).get("NOT_RUN") != 3720 or
            decision.get("planned_cells") != EXPECTED_CELLS or
            candidate.get("rows") != 2495 or candidate.get("not_run") != 876 or
            promotion.get("association_results_consulted") is not False or
            set(not_run_by_trait) != {"adhd", "bipolar", "insomnia", "longsleep", "mdd", "parkinson", "scz"} or
            sum(not_run_by_trait.values()) != 3720):
        raise ValueError("Canonical or candidate source counts do not match the frozen feasibility inputs")

    unchanged_after_long_sleep = sum(
        value for trait, value in not_run_by_trait.items() if trait != "longsleep"
    )
    disorder_not_run = sum(
        value for trait, value in not_run_by_trait.items()
        if trait in {"adhd", "bipolar", "mdd", "parkinson", "scz"}
    )
    candidate_not_run = int(candidate["not_run"])
    family_with_candidate = unchanged_after_long_sleep + candidate_not_run
    scenarios = {
        "observed_candidate_replacing_longsleep": family_with_candidate,
        "perfect_longsleep_replacement_zero_not_run": unchanged_after_long_sleep,
        "perfect_replacements_for_both_sleep_traits_zero_not_run": disorder_not_run,
    }
    if (family_with_candidate != 3305 or unchanged_after_long_sleep != 2429 or
            disorder_not_run != 1758 or candidate_not_run <= FAMILY_LIMIT or
            unchanged_after_long_sleep - FAMILY_LIMIT != 1556 or disorder_not_run - FAMILY_LIMIT != 885):
        raise ValueError("Family feasibility arithmetic differs from expected frozen counts")
    sources = [CANONICAL_DECISION, CANONICAL_BY_TRAIT, CANDIDATE_SUMMARY, PROMOTION_REASSESSMENT]
    result = {
        "schema_version": 1,
        "analysis_id": "brain6_power_optimized_sensitivity_family_feasibility_v1",
        "status": "PASS_SLEEP_ONLY_REPLACEMENTS_CANNOT_PASS_FROZEN_FAMILY_QC",
        "association_results_consulted": False,
        "assumptions": [
            "The seven-trait family size and 5% missingness ceiling remain frozen.",
            "Replacing one sleep trait does not alter eligibility counts for retained traits.",
            "The zero-NOT_RUN cases are theoretical lower bounds, not observed analyses.",
        ],
        "family": {
            "cells": EXPECTED_CELLS,
            "maximum_not_run": FAMILY_LIMIT,
            "canonical_not_run": 3720,
            "candidate_longsleep_replacement_not_run": candidate_not_run,
            "candidate_family_not_run": family_with_candidate,
        },
        "canonical_not_run_by_trait": not_run_by_trait,
        "best_case_lower_bounds": {
            "one_longsleep_replacement_zero_not_run": {
                "not_run": unchanged_after_long_sleep,
                "excess_over_frozen_limit": unchanged_after_long_sleep - FAMILY_LIMIT,
                "passes": unchanged_after_long_sleep <= FAMILY_LIMIT,
            },
            "both_sleep_traits_zero_not_run": {
                "retained_five_disorder_traits_not_run": disorder_not_run,
                "excess_over_frozen_limit": disorder_not_run - FAMILY_LIMIT,
                "passes": disorder_not_run <= FAMILY_LIMIT,
            },
        },
        "scenarios_not_run_cells": scenarios,
        "implication": (
            "No single sleep-trait replacement can pass the frozen seven-trait family QC. "
            "Even perfect local-h2 eligibility for both sleep traits would leave the five "
            "unchanged disorder inputs above the family ceiling. Do not launch bivariate LAVA "
            "or relax the threshold on the basis of another sleep-only candidate."
        ),
        "source_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in sources},
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        if output.read_text(encoding="utf-8") != serialized:
            raise FileExistsError(f"Existing feasibility audit differs; refusing to overwrite {output}")
    else:
        output.write_text(serialized, encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    print(json.dumps(audit(output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
