#!/usr/bin/env python3
"""Independently validate the fixed MDD2025 pilot selection and local locks."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from build_mdd2025_trait_pilot_v1 import (LOCI_FILE, PILOT_CONFIG, PROVENANCE,
                                           RUNNER, SOURCE_CONFIG, loci_by_chromosome)
from screen_mdd2025_lava_source_coverage_v1 import digest

ROOT = Path(__file__).resolve().parents[2]
RULE = ROOT / "brain6/config/lava_multitrait_feasibility_v1/mdd2025_pilot_decision_rule_v1.json"
SCREEN = ROOT / "brain6/results/lava_multitrait_feasibility_v1/mdd2025_source_coverage_v1.json"
RESCUE_MANIFEST = ROOT / "brain6/results/lava_rescue_v1/lava_rescue_v1_manifest.json"
OUTPUT = ROOT / "brain6/results/lava_multitrait_feasibility_v1/mdd2025_pilot_lock_validation_v1.json"


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(OUTPUT)
    pilot = json.loads(PILOT_CONFIG.read_text())
    source = json.loads(SOURCE_CONFIG.read_text())
    screen = json.loads(SCREEN.read_text())
    rule = json.loads(RULE.read_text())
    rescue = json.loads(RESCUE_MANIFEST.read_text())
    loci = loci_by_chromosome()
    selected = []
    for chrom, rows in loci.items():
        ranked = sorted(rows, key=lambda row: hashlib.sha256(f"brain6-mdd2025-pilot-v1:{row[2]}".encode()).digest())
        selected.extend(row[2] for row in ranked[:5])
    selected = sorted(set(selected + ["950"]), key=int)
    checks = {
        "selection_exact": pilot["locus_ids"] == selected and pilot["expected_loci"] == len(selected) == 111,
        "source_config_hash": pilot["source_config_sha256"] == digest(SOURCE_CONFIG),
        "source_identity": pilot["source_sha256"] == source["source_sha256"] == screen["source_sha256"],
        "source_coverage_receipt": screen["status"] == "PASS_SOURCE_COVERAGE_DIAGNOSTIC_ONLY"
                                   and screen["config_sha256"] == digest(SOURCE_CONFIG),
        "cohort_sidecar_identity": pilot["cohort_sidecar_md5"] == source["cohort_sidecar_md5"],
        "locus_hash": pilot["loci_file_sha256"] == digest(LOCI_FILE) == screen["locus_file_sha256"],
        "reference_manifest_identity": pilot["reference_provenance_sha256"] == rescue["reference_provenance_sha256"],
        "materializer_hash": pilot["materializer_sha256"] == digest(ROOT / "brain6/scripts/build_mdd2025_trait_pilot_v1.py"),
        "r_runner_hash": pilot["r_runner_sha256"] == digest(RUNNER),
        "fixed_scientific_parameters": pilot["strict_gate_p"] == 0.05 / 17465
                                       and pilot["worker_count"] == 4
                                       and pilot["locus_processing"] == {
                                           "min_K": 2, "prune_thresh": 99, "max_prop_K": 0.75,
                                           "drop_failed": True, "max_block_size": 3000, "cap_estimates": True},
        "decision_rule_identity": rule["analysis_id"] == pilot["analysis_id"]
                                  and rule["advance_if_any"] == {
                                      "absolute_tested_fraction_gain_at_least": 0.10,
                                      "relative_low_local_h2_not_run_reduction_at_least": 0.25},
        "no_family_promotion": "no canonical/family promotion" in pilot["scope"]
                               and "never family promotion" in rule["decision_scope"],
    }
    status = "PASS_LOCAL_LOCKS_SOURCE_ARCHIVE_NOT_RECHECKED" if all(checks.values()) else "FAIL_LOCAL_LOCKS"
    reference_live = PROVENANCE.is_file()
    if reference_live and digest(PROVENANCE) != pilot["reference_provenance_sha256"]:
        status = "FAIL_LIVE_REFERENCE_HASH"
    result = {"analysis_id": pilot["analysis_id"], "status": status, "checks": checks,
              "pilot_config_sha256": digest(PILOT_CONFIG), "rule_sha256": digest(RULE),
              "source_screen_receipt_sha256": digest(SCREEN),
              "source_archive_currently_available": Path(source["source_path"]).is_file(),
              "reference_provenance_currently_available": reference_live,
              "scope": "Local manifest and code validation only; worker and source archive bytes require mounted SSD"}
    with OUTPUT.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"status": status, "checks_passed": sum(checks.values()),
                      "checks_total": len(checks), "source_archive_available": result["source_archive_currently_available"]},
                     sort_keys=True))
    if not all(checks.values()) or status.startswith("FAIL"):
        raise ValueError("MDD2025 pilot lock validation failed")


if __name__ == "__main__":
    main()
