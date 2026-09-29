#!/usr/bin/env python3
"""Freeze pre-annotation gate decisions for partial PLACO candidate regions.

This ledger is not the final shared-locus table. It records why the current
four-pair, factor-normalized candidate regions cannot receive evidence tiers,
so later annotations cannot influence the gate decision.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
LOCI_PATH = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_loci.tsv"
SENSITIVITY_PROVENANCE = ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/provenance.json"
TIER_POLICY = ROOT / "brain6/config/shared_locus_evidence_tiers_v1.json"
TIER_POLICY_SHA = ROOT / "brain6/config/shared_locus_evidence_tiers_v1.json.sha256"
LAVA_RUN_ID = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
LAVA_DECISION = ROOT / "work/lava-canonical-v3-production" / LAVA_RUN_ID / "canonical_family_decision.json"
OUT_DIR = ROOT / "brain6/results/annotation"
OUT_PATH = OUT_DIR / "partial_candidate_region_gate_v1.tsv"
PROVENANCE_PATH = OUT_DIR / "partial_candidate_region_gate_v1.provenance.json"
EXPECTED_PAIRS = {
    "insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz",
}
FIELDS = [
    "candidate_region_id", "pair_id", "sleep_trait", "brain_disorder", "CHR",
    "START", "STOP", "lead_variants", "n_lead_signals", "lead_P_PLACO",
    "PLACO_evidence_scope", "canonical_LAVA_family_status",
    "five_track_PLACO_complete", "external_genotype_LD_validated",
    "evidence_tier", "tier_decision_reason", "final_shared_locus",
    "annotation_scope",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(out_path: Path = OUT_PATH, provenance_path: Path = PROVENANCE_PATH) -> dict[str, Any]:
    policy = json.loads(TIER_POLICY.read_text(encoding="utf-8"))
    policy_hash = TIER_POLICY_SHA.read_text(encoding="utf-8").split()[0]
    if (policy.get("status") != "FROZEN_PRE_ANNOTATION" or
            policy.get("functional_annotation_reviewed_before_freeze") is not False or
            sha256(TIER_POLICY) != policy_hash):
        raise ValueError("Frozen pre-annotation tier policy is missing, changed, or invalid")

    sensitivity = json.loads(SENSITIVITY_PROVENANCE.read_text(encoding="utf-8"))
    if (sensitivity.get("analysis_id") != "brain6_placo_factor_normalized_sensitivity_v2" or
            sensitivity.get("status") != "DIAGNOSTIC_SENSITIVITY_NOT_PROMOTED" or
            sensitivity.get("scope", {}).get("full_five_pair_family") is not False or
            sensitivity.get("method", {}).get("external_genotype_validation") != "NOT_AVAILABLE" or
            sha256(LOCI_PATH) != next(item["sha256"] for item in sensitivity["outputs"]
                                     if item["name"] == "loci")):
        raise ValueError("Candidate-region source is not the expected unpromoted sensitivity")

    decision = json.loads(LAVA_DECISION.read_text(encoding="utf-8"))
    if (decision.get("overall_status") != "FAILED_QC_NOT_PROMOTED" or
            decision.get("status_counts", {}).get("NOT_RUN") != 3720 or
            decision.get("status_counts", {}).get("TESTED") != 13745):
        raise ValueError("Canonical LAVA gate no longer matches the frozen failed-QC evidence")

    loci = read_tsv(LOCI_PATH)
    if len(loci) != 19 or {row["pair_id"] for row in loci} != EXPECTED_PAIRS:
        raise ValueError("Expected 19 exploratory candidate regions across four admitted pairs")
    keys: set[str] = set()
    output: list[dict[str, str]] = []
    for row in loci:
        key = row["locus_id"]
        if key in keys:
            raise ValueError(f"Duplicate candidate-region key: {key}")
        keys.add(key)
        if row["locus_status"] != "PLACO_ONLY_CANDIDATE_LOCUS":
            raise ValueError(f"Unexpected candidate evidence status: {row['locus_status']}")
        output.append({
            "candidate_region_id": key,
            "pair_id": row["pair_id"],
            "sleep_trait": row["sleep_trait"],
            "brain_disorder": row["brain_disorder"],
            "CHR": row["CHR"], "START": row["START"], "STOP": row["STOP"],
            "lead_variants": row["lead_variants"],
            "n_lead_signals": row["n_lead_signals"],
            "lead_P_PLACO": row["lead_P_PLACO"],
            "PLACO_evidence_scope": "FOUR_PAIR_DIAGNOSTIC_FACTOR_NORMALIZED_SENSITIVITY",
            "canonical_LAVA_family_status": "FAILED_QC_NOT_PROMOTED",
            "five_track_PLACO_complete": "False",
            "external_genotype_LD_validated": "False",
            "evidence_tier": "NOT_TIERED",
            "tier_decision_reason": (
                "Incomplete five-track PLACO family; canonical LAVA family failed frozen QC; "
                "factor-normalized LD is not genotype-validated"
            ),
            "final_shared_locus": "False",
            "annotation_scope": "DESCRIPTIVE_CANDIDATE_REGION_ONLY_AFTER_TIER_GATE",
        })

    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(out_path, output)
    provenance = {
        "schema_version": 1,
        "analysis_id": "brain6_partial_candidate_region_preannotation_gate_v1",
        "status": "PASS_NOT_TIERED_PREREQUISITES_BLOCKED",
        "scope": {
            "candidate_regions": len(output), "pairs": sorted(EXPECTED_PAIRS),
            "five_track_PLACO_complete": False,
            "canonical_LAVA_family_status": "FAILED_QC_NOT_PROMOTED",
            "final_shared_loci_assigned": 0,
            "functional_annotations_used_for_tier_decision": False,
            "descriptive_annotation_permitted_after_gate": True,
        },
        "frozen_rule": {"path": str(TIER_POLICY.relative_to(ROOT)), "sha256": sha256(TIER_POLICY)},
        "sources": {
            "candidate_regions": {"path": str(LOCI_PATH.relative_to(ROOT)), "sha256": sha256(LOCI_PATH)},
            "candidate_sensitivity_provenance": {
                "path": str(SENSITIVITY_PROVENANCE.relative_to(ROOT)),
                "sha256": sha256(SENSITIVITY_PROVENANCE),
            },
            "canonical_lava_decision": {
                "path": str(LAVA_DECISION.relative_to(ROOT)), "sha256": sha256(LAVA_DECISION),
            },
        },
        "builder": {"path": str(Path(__file__).resolve().relative_to(ROOT)),
                    "sha256": sha256(Path(__file__).resolve())},
        "output": {"path": str(out_path.relative_to(ROOT)), "sha256": sha256(out_path),
                   "rows": len(output)},
    }
    provenance_path.parent.mkdir(parents=True, exist_ok=True)
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n",
                               encoding="utf-8")
    return provenance


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT_PATH)
    parser.add_argument("--provenance", type=Path, default=PROVENANCE_PATH)
    args = parser.parse_args()
    print(json.dumps(build(args.output, args.provenance), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
