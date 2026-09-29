#!/usr/bin/env python3
"""Reconcile candidate technical advancement with the frozen family QC gate."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "brain6/results/power_optimized_sensitivity_v1"
INPUT_TABLE = BASE / "candidate_replacement_screen_complete_v2.tsv"
PROMOTION = BASE / "promotion_rule_reassessment_2026-09-26.json"
OVERLAP = BASE / "overlap_integrity_audit_v1.json"
FEASIBILITY = BASE / "sensitivity_family_feasibility_v1.json"
OUTPUT = BASE / "candidate_replacement_screen_complete_v3.tsv"
PROVENANCE = BASE / "candidate_replacement_screen_complete_v3.provenance.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def build(output: Path = OUTPUT, provenance_path: Path = PROVENANCE) -> dict[str, Any]:
    fields, rows = read_tsv(INPUT_TABLE)
    promotion = json.loads(PROMOTION.read_text(encoding="utf-8"))
    overlap = json.loads(OVERLAP.read_text(encoding="utf-8"))
    feasibility = json.loads(FEASIBILITY.read_text(encoding="utf-8"))
    if (len(rows) != 6 or promotion.get("decision_status") != "FULL_SENSITIVITY_TECHNICAL_ADVANCEMENT_JUSTIFIED" or
            promotion.get("association_results_consulted") is not False or
            overlap.get("status") != "PASS_OVERLAP_AUDIT_FAMILY_FAILS_FROZEN_QC" or
            overlap.get("family_decision", {}).get("bivariate_lava_started") is not False or
            feasibility.get("status") != "PASS_SLEEP_ONLY_REPLACEMENTS_CANNOT_PASS_FROZEN_FAMILY_QC"):
        raise ValueError("Candidate decision inputs do not match the current blinded technical/QC state")
    candidate = next(row for row in rows if row["candidate"].startswith("Dashti 2019"))
    criteria = promotion.get("qualifying_criteria", {})
    family = overlap.get("family_decision", {})
    overlap_pairs = family.get("candidate_pair_gate_locus_intersections", {})
    if not all(criteria.get(key) is True for key in (
            "finite_local_h2_processability_gain", "low_local_h2_not_run_reduction")):
        raise ValueError("Candidate no longer meets both predeclared technical advancement examples")
    if (family.get("family_not_run") != 3305 or family.get("family_cells") != 17465 or
            family.get("maximum_not_run_cells") != 873 or family.get("passes_frozen_5_percent_qc") is not False or
            family.get("bivariate_lava_started") is not False or set(overlap_pairs) != {"bipolar", "parkinson", "scz"}):
        raise ValueError("Candidate family QC or pair-specific overlap preparation has changed")
    candidate["decision"] = (
        "TECHNICAL ADVANCEMENT CRITERIA MET: +16.63 percentage points finite local-h2 processability "
        "and 31.08% fewer low-local-h2 NOT_RUN loci. Candidate-specific LDSC overlap preparation passed. "
        "Bivariate LAVA remains blocked because the substituted seven-trait family has 3,305/17,465 "
        "NOT_RUN cells versus the frozen 873-cell ceiling; no replacement local-rg result is promoted."
    )
    serialized_rows = []
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() or provenance_path.exists():
        raise FileExistsError("Refusing to overwrite a versioned candidate decision artifact")
    with output.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    sources = (INPUT_TABLE, PROMOTION, OVERLAP, FEASIBILITY, Path(__file__).resolve())
    try:
        output_name = str(output.resolve().relative_to(ROOT))
    except ValueError:
        output_name = str(output.resolve())
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_candidate_replacement_decision_reassessment_v1",
        "status": "PASS_TECHNICAL_ADVANCEMENT_FAMILY_QC_BLOCKED",
        "decision_rule": "Phase 5 advancement is met when either predeclared objective technical criterion qualifies; this is separate from the frozen seven-trait family QC gate.",
        "association_results_consulted": False,
        "technical_advancement_candidates": 1,
        "candidate_specific_overlap_matrices_passed": len(overlap_pairs),
        "substituted_family_not_run": family["family_not_run"],
        "substituted_family_cells": family["family_cells"],
        "substituted_family_maximum_not_run": family["maximum_not_run_cells"],
        "bivariate_lava_started": False,
        "full_family_qc_passing_candidates": 0,
        "downstream_loci_eligible": 0,
        "sleep_only_best_case_both_traits_excess": feasibility["best_case_lower_bounds"]["both_sleep_traits_zero_not_run"]["excess_over_frozen_limit"],
        "input_sha256": {str(path.relative_to(ROOT)): sha256(path) for path in sources},
        "output_sha256": {output_name: sha256(output)},
        "candidate_count": len(rows),
    }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE)
    args = parser.parse_args()
    print(json.dumps(build(args.output, args.provenance), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
