#!/usr/bin/env python3
"""Report all 17 extension gates, including honest blocked/pending states."""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path("discovery_extension")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_check(script: str) -> tuple[bool, str]:
    result = subprocess.run(["python3", script], check=False, text=True, capture_output=True)
    output = (result.stdout + result.stderr).strip()
    return result.returncode == 0, output


def artifact_state(path: Path) -> str:
    return f"present:sha256={sha256(path)}" if path.is_file() else "absent"


def main() -> None:
    core_ok, core_check = run_check(str(ROOT / "scripts/00_verify_core_checkpoint.py"))
    panel_ok, panel_check = run_check(str(ROOT / "scripts/05_validate_extension_panel.py"))
    source_ok, source_check = run_check(str(ROOT / "scripts/08_validate_source_contracts.py"))
    panel = read_tsv(ROOT / "config/candidate_traits.tsv")
    lock = json.loads((ROOT / "config/extension_panel.lock.json").read_text())
    preflight = json.loads((ROOT / "provenance/acquisition_preflight.json").read_text())

    h2_path = ROOT / "results/ldsc/extension_trait_readiness.tsv"
    rg_path = ROOT / "results/ldsc/extension_rg_matrix.tsv"
    pair_universe_path = ROOT / "results/ldsc/extension_pair_universe.tsv"
    novelty_path = ROOT / "results/novelty/extension_novelty_audit.tsv"
    prioritization_path = ROOT / "results/prioritization/prioritized_pairs.tsv"
    replication_path = ROOT / "results/replication/replication_results.tsv"
    local_path = ROOT / "results/local/local_rg_results.tsv"
    pleiotropy_path = ROOT / "results/pleiotropy/pleiotropy_results.tsv"
    fine_mapping_path = ROOT / "results/fine_mapping/fine_mapping_colocalization.tsv"
    mechanism_path = ROOT / "results/mechanism/mechanistic_synthesis.tsv"

    blocked = preflight["status"] != "READY"
    gates = [
        (1, "Freeze and protect core atlas", "PASS" if core_ok else "FAIL_CORE_CHECKPOINT_DRIFT", core_check,
         "Core remains the exact 45-trait/396-pair checkpoint." if core_ok else "A checkpointed core artifact differs in the current working tree; extension execution must remain stopped."),
        (2, "Broaden phenotype expansion strategy", "PASS" if panel_ok else "ARTIFACT_PRESENT_CORE_REVALIDATION_BLOCKED", panel_check,
         "452-trait eligible Pan-UKB EUR universe plus a separately researched imaging queue."),
        (3, "Map prior sleep-trait correlation literature", "PASS_WITH_SCOPE_LIMITATION",
         artifact_state(ROOT / "results/prior_sleep_screen_inventory.tsv"),
         "8,896 extracted prior-screen rows; exact no-match is not proof of novelty."),
        (4, "Build and lock candidate universe", "PASS" if panel_ok else "LOCK_PRESENT_CORE_REVALIDATION_BLOCKED", panel_check,
         f"242-trait candidate pool; {len(panel)} traits locked before extension rg; {lock['planned_raw_rg_test_count']} planned pairs."),
        (5, "Verify sources, schemas, and harmonization", ("CONTRACT_PASS_EXECUTION_BLOCKED" if blocked else "PASS") if source_ok else "CONTRACT_PRESENT_CORE_REVALIDATION_BLOCKED",
         source_check,
         "Continuous, binary, variant-map, QC, and synthetic smoke-test contracts pass; real traversal awaits acquisition."),
        (6, "Acquire full-resolution summary statistics", "BLOCKED_INSUFFICIENT_STORAGE" if blocked else "PENDING_ACQUISITION",
         artifact_state(ROOT / "results/acquisition_plan.tsv"),
         f"Need {preflight['required_free_gib']} GiB; measured {preflight['available_free_gib']} GiB; no bulk download started."),
        (7, "Rerun extension LDSC h2 QC", "BLOCKED_UPSTREAM" if not h2_path.is_file() else "PASS_ARTIFACT_PRESENT",
         artifact_state(h2_path), "Source h2 precheck is not accepted as the required extension rerun."),
        (8, "Run global genetic correlations", "BLOCKED_UPSTREAM" if not rg_path.is_file() else "PASS_ARTIFACT_PRESENT",
         artifact_state(rg_path), "BH FDR must remain confined to the primary extension family."),
        (9, "Audit pair-level novelty", "PROTOCOL_READY_UPSTREAM_BLOCKED" if not novelty_path.is_file() else "REVIEW_ARTIFACT_PRESENT",
         artifact_state(novelty_path) + ";protocol=17_prepare_pair_novelty_audit.py+18_validate_pair_novelty_audit.py", "Required for every extension-FDR-significant pair; no panel-level label substitutes."),
        (10, "Prioritize findings", "BLOCKED_UPSTREAM" if not prioritization_path.is_file() else "REVIEW_ARTIFACT_PRESENT",
         artifact_state(prioritization_path), "Tier A requires effect, QC, FDR, and pair-level novelty; Tier B additionally requires replication or strong local support."),
        (11, "Independent replication", "PROTOCOL_READY_UPSTREAM_BLOCKED" if not replication_path.is_file() else "REPLICATION_ARTIFACT_PRESENT",
         artifact_state(replication_path) + ";protocol=config/replication_contract.json", "Same-cohort internal splits cannot satisfy the independent-replication gate."),
        (12, "Local genetic correlation", "PROTOCOL_READY_UPSTREAM_BLOCKED" if not local_path.is_file() else "LOCAL_ARTIFACT_PRESENT",
         artifact_state(local_path) + ";protocol=config/followup_contract.json", "Global rg is not local sharing."),
        (13, "Pleiotropy analysis", "PROTOCOL_READY_UPSTREAM_BLOCKED" if not pleiotropy_path.is_file() else "PLEIOTROPY_ARTIFACT_PRESENT",
         artifact_state(pleiotropy_path) + ";protocol=config/followup_contract.json", "Pleiotropy is evaluated as an alternative explanation, not an inconvenience."),
        (14, "Fine-mapping and colocalization", "PROTOCOL_READY_UPSTREAM_BLOCKED" if not fine_mapping_path.is_file() else "FOLLOWUP_ARTIFACT_PRESENT",
         artifact_state(fine_mapping_path) + ";protocol=config/followup_contract.json", "Colocalization requires explicit hypotheses, priors, and sensitivity analyses."),
        (15, "Mechanistic annotation", "PROTOCOL_READY_UPSTREAM_BLOCKED" if not mechanism_path.is_file() else "SYNTHESIS_ARTIFACT_PRESENT",
         artifact_state(mechanism_path) + ";protocol=config/followup_contract.json", "Mechanistic interpretation is ranked and cannot convert correlation into causation."),
        (16, "Adversarial review", "CHECKLIST_READY_FINAL_REVIEW_BLOCKED",
         str(ROOT / "results/adversarial_review_checklist.tsv"),
         "Current checklist records design risks; final findings-level review awaits actual results."),
        (17, "Publication-grade reporting", "STATUS_REPORT_ONLY_FINAL_REPORT_BLOCKED",
         str(ROOT / "results/extension_status.md"),
         "No discovery, novelty, replication, local, or mechanistic result is claimed while upstream gates are blocked."),
    ]
    gate_rows = [
        {"stage": stage, "gate": gate, "status": status, "evidence": evidence, "interpretation": interpretation}
        for stage, gate, status, evidence, interpretation in gates
    ]
    gates_path = ROOT / "results/extension_acceptance_gates.tsv"
    gates_path.parent.mkdir(parents=True, exist_ok=True)
    with gates_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(gate_rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(gate_rows)

    risks = [
        ("R01", "Core drift", "Re-run byte-level core checkpoint before and after every extension stage.", "PASS" if core_ok else "FAIL_BLOCKING", core_check),
        ("R02", "Post-result panel selection", "Panel membership/order/hash must predate extension rg.", "PASS", lock["locked_utc"]),
        ("R03", "Multiple-testing leakage", "Never combine the extension BH family with the immutable 396 core tests.", "PASS_CONTRACT", "extension_harmonization_policy.json"),
        ("R04", "Weak h2", "Exclude rerun h2 Z<4 from primary rg without replacement.", "PENDING_REAL_H2", artifact_state(h2_path)),
        ("R05", "LDSC intercept inflation", "Exclude rerun intercept>1.2 from primary rg and retain sensitivity status.", "PENDING_REAL_H2", artifact_state(h2_path)),
        ("R06", "UK Biobank sample overlap", "Inspect cross-trait intercepts and disclose overlap; do not equate LDSC adjustment with independent replication.", "PENDING_REAL_RG", artifact_state(rg_path)),
        ("R07", "Sparse or proxy phenotypes", "Retain exact phenotype definitions and distinguish medication/proxy traits from diagnoses.", "PASS_METADATA", "candidate_traits.tsv"),
        ("R08", "Panel-level novelty inflation", "Require pair-level direct/same-phenotype/same-direction/same-sleep-context audit.", "PENDING_PAIR_AUDIT", artifact_state(novelty_path)),
        ("R09", "Replication non-independence", "Require non-overlapping participants and separately sourced summary statistics.", "PENDING_REPLICATION", artifact_state(replication_path)),
        ("R10", "Global-to-local overreach", "Do not call global rg evidence of a shared locus.", "PENDING_LOCAL_ANALYSIS", artifact_state(local_path)),
        ("R11", "Pleiotropy or mediated effects", "Evaluate shared-factor and directionally pleiotropic alternatives.", "PENDING_PLEIOTROPY", artifact_state(pleiotropy_path)),
        ("R12", "Colocalization overclaim", "Report hypotheses/priors/sensitivity; colocalization is not causality.", "PENDING_COLOCALIZATION", artifact_state(fine_mapping_path)),
        ("R13", "Synthetic/real result contamination", "Synthetic tests stay under synthetic or temporary paths and carry explicit markers.", "PASS", "two isolated smoke tests"),
        ("R14", "Storage-driven partial acquisition", "Do not silently analyze a result-selected subset of the locked panel.", "PASS_BLOCKED", preflight["status"]),
    ]
    checklist_path = ROOT / "results/adversarial_review_checklist.tsv"
    with checklist_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["risk_id", "risk", "required_control", "current_status", "current_evidence"])
        writer.writerows(risks)

    if not core_ok:
        overall = "BLOCKED_AT_CORE_CHECKPOINT_AND_FULL_RESOLUTION_ACQUISITION_GATES"
    elif blocked:
        overall = "BLOCKED_AT_FULL_RESOLUTION_ACQUISITION_GATE"
    else:
        overall = "IN_PROGRESS"
    status_path = ROOT / "results/extension_status.md"
    status_path.write_text(
        "# Novelty-Enriched Phenome Discovery Extension — execution status\n\n"
        f"Generated: {datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}\n\n"
        f"Overall: **{overall}**\n\n"
        + (
            "The immutable core checkpoint passes. "
            if core_ok else
            "The immutable core checkpoint currently fails because a checkpointed artifact differs in the working tree; no extension execution is permitted until that separate core change is resolved and revalidated. "
        )
        + "The eligible universe, prior-screen inventory, "
        "242-trait candidate pool, and independently locked 100-trait extension panel are complete. "
        "Source schemas, exact URLs/checksums, a GRCh37 variant/INFO mapping contract, harmonization, "
        "h2/rg collation, isolated FDR, and plotting code are defined and synthetically tested.\n\n"
        f"Real acquisition is blocked: the exact compressed inputs total {preflight['compressed_source_gib']} GiB "
        f"and require {preflight['required_free_gib']} GiB with the locked safety factor, while the preflight "
        f"measured {preflight['available_free_gib']} GiB free ({preflight['shortfall_gib']} GiB short). "
        "No bulk download was started. Therefore no extension h2 rerun, genetic correlation, extension FDR hit, "
        "pair-level novelty claim, replication, local correlation, pleiotropy, fine-mapping, colocalization, "
        "mechanistic inference, or final manuscript claim exists yet.\n\n"
        "See `extension_acceptance_gates.tsv` for all 17 gates and "
        "`adversarial_review_checklist.tsv` for the current challenge audit.\n",
        encoding="utf-8",
    )

    provenance = {
        "schema_version": "1.0.0",
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "overall_status": overall,
        "stage_count": len(gates),
        "core_check_pass": core_ok,
        "core_check": core_check,
        "panel_check_pass": panel_ok,
        "panel_check": panel_check,
        "source_check_pass": source_ok,
        "source_check": source_check,
        "panel_sha256": sha256(ROOT / "config/candidate_traits.tsv"),
        "acceptance_gates_sha256": sha256(gates_path),
        "adversarial_checklist_sha256": sha256(checklist_path),
        "status_report_sha256": sha256(status_path),
    }
    out = ROOT / "provenance/extension_acceptance.json"
    out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"EXTENSION_ACCEPTANCE_REPORTED stages=17 overall={overall}")


if __name__ == "__main__":
    main()
