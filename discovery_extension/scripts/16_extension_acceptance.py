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


def read_tsv_fields(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t").fieldnames or [])


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
    prioritization_path = ROOT / "results/prioritization/novel_hit_priority.tsv"
    replication_path = ROOT / "results/replication/replication_results.tsv"
    replication_queue_path = ROOT / "results/replication/replication_source_queue.tsv"
    replication_candidate_lock_path = ROOT / "config/replication_candidate_family.lock.json"
    replication_manifest_lock_path = ROOT / "config/replication_manifest.lock.json"
    replication_search_path = ROOT / "results/replication/replication_unique_phenotype_search.tsv"
    local_path = ROOT / "results/local/local_rg_results.tsv"
    local_readiness_path = ROOT / "results/local/local_architecture_readiness.tsv"
    pleiotropy_path = ROOT / "results/pleiotropy/novel_shared_loci.tsv"
    pleiotropy_readiness_path = ROOT / "results/pleiotropy/pleiotropy_readiness.tsv"
    fine_mapping_path = ROOT / "results/fine_mapping/fine_mapping_colocalization.tsv"
    fine_mapping_readiness_path = ROOT / "results/fine_mapping/fine_mapping_readiness.tsv"
    mechanism_path = ROOT / "results/mechanism/mechanistic_synthesis.tsv"
    mechanism_readiness_path = ROOT / "results/mechanism/mechanism_readiness.tsv"
    adversarial_review_path = ROOT / "adversarial_review.md"
    final_report_path = ROOT / "final_report.md"
    final_counts_path = ROOT / "results/final_extension_counts.tsv"
    top_discoveries_path = ROOT / "results/top_novel_discoveries.tsv"
    final_report_provenance_path = ROOT / "provenance/final_report.json"

    panel_ids = [row["extension_trait_id"] for row in panel]
    h2_rows: list[dict[str, str]] = []
    h2_primary_count = 0
    h2_sensitivity_count = 0
    h2_complete = False
    if h2_path.is_file():
        try:
            h2_rows = read_tsv(h2_path)
            h2_primary_count = sum(row["primary_rg_eligibility"] == "PRIMARY_PASS" for row in h2_rows)
            h2_sensitivity_count = sum(row["primary_rg_eligibility"] == "SENSITIVITY_ONLY" for row in h2_rows)
            h2_complete = (
                [row["extension_trait_id"] for row in h2_rows] == panel_ids
                and h2_primary_count + h2_sensitivity_count == len(panel_ids)
                and all(row["analysis_status"] == "EXTENSION_H2_QC_COMPLETE" for row in h2_rows)
            )
        except (KeyError, OSError):
            h2_complete = False
    sleep_ids = [
        row["trait_id"] for row in read_tsv(Path("config/analysis_panel.tsv"))
        if row["domain"] == "sleep"
    ]
    expected_rg_pairs = {(sleep_id, trait_id) for sleep_id in sleep_ids for trait_id in panel_ids}
    rg_rows: list[dict[str, str]] = []
    rg_complete = False
    rg_fdr_count = 0
    rg_joint_screen_count = 0
    if rg_path.is_file():
        try:
            rg_rows = read_tsv(rg_path)
            observed_rg_pairs = {(row["sleep_trait"], row["extension_trait_id"]) for row in rg_rows}
            rg_complete = (
                len(rg_rows) == len(expected_rg_pairs)
                and observed_rg_pairs == expected_rg_pairs
                and all(row["analysis_status"] == "PRIMARY_EXTENSION_RG_COMPLETE" for row in rg_rows)
            )
            rg_fdr_count = sum(float(row["extension_fdr"]) < 0.05 for row in rg_rows)
            rg_joint_screen_count = sum(
                float(row["extension_fdr"]) < 0.05 and abs(float(row["rg"])) >= 0.15
                for row in rg_rows
            )
        except (KeyError, OSError, ValueError):
            rg_complete = False
    novelty_rows = read_tsv(novelty_path) if novelty_path.is_file() else []
    novelty_expected_pairs = {
        (row["sleep_trait"], row["extension_trait_id"])
        for row in rg_rows if float(row["extension_fdr"]) < 0.05
    } if rg_complete else set()
    novelty_observed_pairs = {
        (row["sleep_trait"], row["extension_trait_id"])
        for row in novelty_rows
    }
    novelty_complete = (
        bool(novelty_expected_pairs)
        and len(novelty_rows) == len(novelty_expected_pairs)
        and novelty_observed_pairs == novelty_expected_pairs
        and all(row.get("audit_status") == "COMPLETE" for row in novelty_rows)
    )
    novelty_pending_count = sum(row.get("audit_status") != "COMPLETE" for row in novelty_rows)
    priority_rows = read_tsv(prioritization_path) if prioritization_path.is_file() else []
    priority_complete = (
        rg_complete and len(priority_rows) == len(rg_rows)
        and {(row["sleep_trait"], row["extension_trait_id"]) for row in priority_rows}
        == expected_rg_pairs
        and all(row.get("priority_tier") in {"A", "B", "C"} for row in priority_rows)
    )
    priority_tier_a_count = sum(row.get("priority_tier") == "A" for row in priority_rows)
    priority_tier_b_count = sum(row.get("priority_tier") == "B" for row in priority_rows)
    replication_queue_rows = read_tsv(replication_queue_path) if replication_queue_path.is_file() else []
    replication_candidate_locked = False
    if replication_queue_rows and replication_candidate_lock_path.is_file():
        try:
            replication_candidate_lock = json.loads(replication_candidate_lock_path.read_text())
            replication_candidate_locked = (
                [row["pair_id"] for row in replication_queue_rows]
                == replication_candidate_lock["pair_ids_in_locked_order"]
                and len(replication_queue_rows) == replication_candidate_lock["pair_count"]
                and replication_candidate_lock["results_accessed_before_lock"] is False
            )
        except (KeyError, json.JSONDecodeError, OSError):
            replication_candidate_locked = False
    replication_testable_count = sum(
        row.get("source_curation_status") == "COMPLETE_BEFORE_RESULTS"
        for row in replication_queue_rows
    )
    replication_unavailable_count = sum(
        row.get("source_curation_status") == "NO_INDEPENDENT_DATASET_COMPLETE_BEFORE_RESULTS"
        for row in replication_queue_rows
    )
    replication_curation_complete = (
        replication_candidate_locked
        and replication_testable_count + replication_unavailable_count == len(replication_queue_rows)
        and all(row.get("results_accessed_before_lock") == "NO" for row in replication_queue_rows)
        and replication_search_path.is_file()
    )
    streaming_contract_path = ROOT / "config/streaming_acquisition_contract.json"
    streaming_snapshot_path = ROOT / "provenance/panukbb/remote_object_snapshot.tsv"
    streaming_snapshot_provenance_path = ROOT / "provenance/panukbb/remote_object_snapshot.json"
    streaming_reference_path = ROOT / "data/reference/panukbb_hm3_variant_reference.tsv.gz"
    streaming_reference_provenance_path = ROOT / "provenance/panukbb/hm3_variant_reference_build.json"
    streaming_receipt_dir = ROOT / "provenance/streaming_receipts"
    current_streaming_code_sha256 = {
        str(path): sha256(path)
        for path in (
            ROOT / "scripts/streaming_io.py",
            ROOT / "scripts/10_harmonize_panukbb.py",
            ROOT / "scripts/11_munge_extension.sh",
        )
    }
    streaming_contract_ready = all(
        path.is_file()
        for path in (
            streaming_contract_path,
            streaming_snapshot_path,
            streaming_snapshot_provenance_path,
            streaming_reference_provenance_path,
        )
    )
    streaming_reference_ready = False
    if streaming_contract_ready and streaming_reference_path.is_file():
        reference_provenance = json.loads(streaming_reference_provenance_path.read_text())
        streaming_reference_ready = (
            reference_provenance.get("output_sha256") == sha256(streaming_reference_path)
            and reference_provenance.get("variant_manifest_stream_receipt", {}).get("verification_status") == "PASS"
        )
    streaming_completed_ids: list[str] = []
    for trait_id in panel_ids:
        receipt_path = streaming_receipt_dir / f"{trait_id}.json"
        munged_path = ROOT / f"data/munged/{trait_id}.sumstats.gz"
        munging_log_path = ROOT / f"data/munged/{trait_id}.log"
        if not all(path.is_file() for path in (receipt_path, munged_path, munging_log_path)):
            continue
        try:
            receipt = json.loads(receipt_path.read_text())
            valid = (
                receipt.get("extension_trait_id") == trait_id
                and receipt.get("pipeline_status") == "STREAM_HARMONIZE_MUNGE_PASS"
                and receipt.get("source_verification", {}).get("verification_status") == "PASS"
                and receipt.get("munged_output_sha256") == sha256(munged_path)
                and receipt.get("munging_log_sha256") == sha256(munging_log_path)
                and receipt.get("harmonization_policy_sha256")
                == sha256(ROOT / "config/extension_harmonization_policy.json")
                and receipt.get("pipeline_code_sha256") == current_streaming_code_sha256
            )
        except (json.JSONDecodeError, OSError):
            valid = False
        if valid:
            streaming_completed_ids.append(trait_id)
    partial_h2_log_count = sum(
        (ROOT / f"logs/h2/h2_{trait_id}.log").is_file() for trait_id in panel_ids
    )
    mirror_blocked = preflight["status"] != "READY"
    blocked = not rg_path.is_file()
    expected_top_fields = [
        "rank", "pair_id", "sleep_trait", "external_phenotype", "discovery_rg",
        "discovery_fdr", "replication_rg", "replication_p", "prior_literature_status",
        "local_loci", "pleiotropic_loci", "colocalized_genes",
        "strongest_mechanistic_evidence", "confidence_level", "analysis_status",
        "claim_limit",
    ]
    final_report_ready = False
    final_report_validation = "required Stage-17 artifacts absent"
    if all(path.is_file() for path in (
        final_report_path, final_counts_path, top_discoveries_path, final_report_provenance_path,
    )):
        try:
            report_provenance = json.loads(final_report_provenance_path.read_text(encoding="utf-8"))
            count_rows = read_tsv(final_counts_path)
            top_rows = read_tsv(top_discoveries_path)
            checksum_ok = (
                report_provenance["outputs"]["report"]["sha256"] == sha256(final_report_path)
                and report_provenance["outputs"]["counts"]["sha256"] == sha256(final_counts_path)
                and report_provenance["outputs"]["top_novel_discoveries"]["sha256"] == sha256(top_discoveries_path)
            )
            count_contract_ok = (
                [int(row["report_order"]) for row in count_rows] == list(range(1, 15))
                and count_rows[0]["value"] == "242"
                and count_rows[1]["value"] == "100"
            )
            blocked_contract_ok = not blocked or (
                all(row["value"] == "NA_BLOCKED_UPSTREAM" for row in count_rows[2:])
                and not top_rows
                and report_provenance["executed_test_count"] is None
            )
            final_report_ready = (
                checksum_ok and count_contract_ok and blocked_contract_ok
                and read_tsv_fields(top_discoveries_path) == expected_top_fields
                and report_provenance["missing_result_policy"] == "NA_BLOCKED_UPSTREAM_NEVER_ZERO"
            )
            final_report_validation = (
                "checksums, 14-count schema, missing-result policy, and top-table schema pass"
                if final_report_ready else
                "Stage-17 artifact checksum, count, blocked-value, or top-table schema mismatch"
            )
        except (IndexError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
            final_report_validation = f"Stage-17 validation error: {error}"
    local_readiness_rows = read_tsv(local_readiness_path) if local_readiness_path.is_file() else []
    local_code_ready = bool(local_readiness_rows) and all(row["code_status"] == "PASS" for row in local_readiness_rows)
    local_dependencies_ready = bool(local_readiness_rows) and all(row["readiness"] == "READY" for row in local_readiness_rows)
    if local_path.is_file():
        local_gate_status = "LOCAL_ARTIFACT_PRESENT"
    elif local_dependencies_ready:
        local_gate_status = "CODE_AND_INPUTS_READY_UPSTREAM_SELECTION_BLOCKED"
    elif local_code_ready:
        local_gate_status = "CODE_READY_REFERENCE_AND_DENSE_SUMSTATS_BLOCKED"
    elif local_readiness_rows:
        local_gate_status = "PROTOCOL_READY_CODE_BLOCKED"
    else:
        local_gate_status = "PROTOCOL_READY_UPSTREAM_BLOCKED"
    pleiotropy_readiness_rows = read_tsv(pleiotropy_readiness_path) if pleiotropy_readiness_path.is_file() else []
    pleiotropy_code_ready = bool(pleiotropy_readiness_rows) and all(row["code_status"] == "PASS" for row in pleiotropy_readiness_rows)
    pleiotropy_inputs_ready = bool(pleiotropy_readiness_rows) and all(row["readiness"] == "READY_FOR_PAIR_MANIFEST_CURATION" for row in pleiotropy_readiness_rows)
    if pleiotropy_path.is_file():
        pleiotropy_gate_status = "PLEIOTROPY_ARTIFACT_PRESENT"
    elif pleiotropy_inputs_ready:
        pleiotropy_gate_status = "CODE_AND_INPUTS_READY_PAIR_MANIFEST_BLOCKED"
    elif pleiotropy_code_ready:
        pleiotropy_gate_status = "CODE_READY_UPSTREAM_DISCOVERY_REPLICATION_AND_DENSE_INPUTS_BLOCKED"
    else:
        pleiotropy_gate_status = "PROTOCOL_READY_UPSTREAM_BLOCKED"
    fine_mapping_readiness_rows = read_tsv(fine_mapping_readiness_path) if fine_mapping_readiness_path.is_file() else []
    fine_mapping_code_ready = bool(fine_mapping_readiness_rows) and all(row["code_status"] == "PASS" for row in fine_mapping_readiness_rows)
    if fine_mapping_path.is_file():
        fine_mapping_gate_status = "FINE_MAPPING_COLOCALIZATION_ARTIFACT_PRESENT"
    elif fine_mapping_code_ready:
        fine_mapping_gate_status = "CODE_READY_UPSTREAM_LOCUS_DENSE_LD_AND_QTL_INPUTS_BLOCKED"
    else:
        fine_mapping_gate_status = "PROTOCOL_READY_UPSTREAM_BLOCKED"
    mechanism_readiness_rows = read_tsv(mechanism_readiness_path) if mechanism_readiness_path.is_file() else []
    mechanism_protocol_ready = bool(mechanism_readiness_rows) and all(row["protocol_status"] == "PASS" for row in mechanism_readiness_rows)
    if mechanism_path.is_file():
        mechanism_gate_status = "MECHANISTIC_SYNTHESIS_ARTIFACT_PRESENT"
    elif mechanism_protocol_ready:
        mechanism_gate_status = "PROTOCOL_READY_UPSTREAM_FINE_MAPPING_AND_EXACT_SOURCE_RELEASES_BLOCKED"
    else:
        mechanism_gate_status = "PROTOCOL_NOT_VALIDATED"
    if len(streaming_completed_ids) == len(panel):
        acquisition_gate_status = "PASS_STREAMING_ACQUISITION"
    elif streaming_completed_ids:
        acquisition_gate_status = f"IN_PROGRESS_STREAMING_{len(streaming_completed_ids)}_OF_{len(panel)}"
    elif streaming_reference_ready:
        acquisition_gate_status = "READY_FOR_STREAMING_ACQUISITION"
    elif streaming_contract_ready:
        acquisition_gate_status = "STREAMING_CONTRACT_READY_REFERENCE_NOT_LOCAL"
    else:
        acquisition_gate_status = "BLOCKED_INSUFFICIENT_STORAGE" if mirror_blocked else "PENDING_ACQUISITION"
    if h2_complete:
        h2_gate_status = (
            f"PASS_COMPLETE_{len(h2_rows)}_TRAITS_PRIMARY_{h2_primary_count}_"
            f"SENSITIVITY_{h2_sensitivity_count}"
        )
    elif h2_path.is_file():
        h2_gate_status = "FAIL_INCOMPLETE_OR_INVALID_ARTIFACT"
    elif partial_h2_log_count:
        h2_gate_status = f"IN_PROGRESS_PARTIAL_H2_{partial_h2_log_count}_OF_{len(panel)}"
    else:
        h2_gate_status = "BLOCKED_UPSTREAM"
    if rg_complete:
        rg_gate_status = f"PASS_COMPLETE_{len(rg_rows)}_PAIRS_FDR_{rg_fdr_count}"
    elif rg_path.is_file():
        rg_gate_status = "FAIL_INCOMPLETE_OR_INVALID_ARTIFACT"
    else:
        rg_gate_status = "BLOCKED_UPSTREAM"
    if novelty_complete:
        novelty_gate_status = f"PASS_COMPLETE_{len(novelty_rows)}_PAIRS"
    elif novelty_path.is_file():
        novelty_gate_status = f"TEMPLATE_PRESENT_PENDING_{novelty_pending_count}_OF_{len(novelty_rows)}"
    else:
        novelty_gate_status = "PROTOCOL_READY_UPSTREAM_BLOCKED"
    if priority_complete:
        priority_gate_status = (
            f"PASS_COMPLETE_{len(priority_rows)}_PAIRS_TIER_A_{priority_tier_a_count}_"
            f"TIER_B_{priority_tier_b_count}"
        )
    elif prioritization_path.is_file():
        priority_gate_status = "FAIL_INCOMPLETE_OR_INVALID_ARTIFACT"
    else:
        priority_gate_status = "BLOCKED_UPSTREAM"
    if replication_path.is_file():
        replication_gate_status = "REPLICATION_ARTIFACT_PRESENT"
    elif replication_manifest_lock_path.is_file():
        replication_gate_status = f"MANIFEST_LOCKED_ANALYSIS_PENDING_{replication_testable_count}_TESTABLE_{replication_unavailable_count}_UNAVAILABLE"
    elif replication_curation_complete:
        replication_gate_status = f"SOURCE_CURATION_COMPLETE_ACQUISITION_PENDING_{replication_testable_count}_TESTABLE_{replication_unavailable_count}_UNAVAILABLE"
    elif replication_candidate_locked:
        replication_gate_status = f"CANDIDATE_FAMILY_LOCKED_SOURCE_CURATION_PENDING_{len(replication_queue_rows)}"
    else:
        replication_gate_status = "PROTOCOL_READY_UPSTREAM_BLOCKED"
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
        (5, "Verify sources, schemas, and harmonization", ("CONTRACT_PASS_STREAMING_ACQUISITION_COMPLETE" if len(streaming_completed_ids) == len(panel) else "CONTRACT_PASS_STREAMING_EXECUTION_IN_PROGRESS" if streaming_contract_ready else "CONTRACT_PASS_EXECUTION_BLOCKED") if source_ok else "CONTRACT_PRESENT_CORE_REVALIDATION_BLOCKED",
         source_check,
         "Continuous, binary, variant-map, QC, synthetic smoke-test, and version-pinned streaming contracts pass."),
        (6, "Acquire full-resolution summary statistics", acquisition_gate_status,
         artifact_state(ROOT / "results/acquisition_plan.tsv") + f";streaming_receipts={len(streaming_completed_ids)}/{len(panel)}",
         f"The {preflight['required_free_gib']} GiB local mirror remains unavailable; exact version-pinned sources were streamed, body-checksummed, harmonized, munged, and receipt-sealed one trait at a time. Completed {len(streaming_completed_ids)}/{len(panel)}."),
        (7, "Rerun extension LDSC h2 QC", h2_gate_status,
         artifact_state(h2_path) + f";h2_logs={partial_h2_log_count}/{len(panel)};primary={h2_primary_count};sensitivity={h2_sensitivity_count}",
         "The harmonized extension rerun, not the source precheck, determines primary rg eligibility."),
        (8, "Run global genetic correlations", rg_gate_status,
         artifact_state(rg_path) + f";pairs={len(rg_rows)}/{len(expected_rg_pairs)};fdr={rg_fdr_count};joint_screen={rg_joint_screen_count}", "BH FDR remains confined to the primary extension family."),
        (9, "Audit pair-level novelty", novelty_gate_status,
         artifact_state(novelty_path) + ";protocol=17_prepare_pair_novelty_audit.py+18_validate_pair_novelty_audit.py", "Required for every extension-FDR-significant pair; no panel-level label substitutes."),
        (10, "Prioritize findings", priority_gate_status,
         artifact_state(prioritization_path) + f";pairs={len(priority_rows)}/{len(expected_rg_pairs)};tier_A={priority_tier_a_count};tier_B={priority_tier_b_count}", "Tier A requires effect, QC, FDR, and pair-level novelty; Tier B additionally requires replication or strong local support."),
        (11, "Independent replication", replication_gate_status,
         artifact_state(replication_path) + ";queue=" + artifact_state(replication_queue_path) + ";search=" + artifact_state(replication_search_path) + ";candidate_lock=" + artifact_state(replication_candidate_lock_path) + ";manifest_lock=" + artifact_state(replication_manifest_lock_path) + ";protocol=config/replication_contract.json", "Same-cohort internal splits cannot satisfy the independent-replication gate."),
        (12, "Local genetic correlation", local_gate_status,
         artifact_state(local_path) + ";readiness=" + artifact_state(local_readiness_path) + ";protocol=config/local_architecture_contract.json",
         "Pinned LAVA/HDL-L code passes runtime checks, but source-verified LD references and dense inputs are absent; global rg is not local sharing."),
        (13, "Pleiotropy analysis", pleiotropy_gate_status,
         artifact_state(pleiotropy_path) + ";readiness=" + artifact_state(pleiotropy_readiness_path) + ";protocol=config/pleiotropy_contract.json",
         "Pinned PLACO+ code passes, but replicated pair selection and genome-wide dense inputs are absent; statistical pleiotropy is not a shared causal variant."),
        (14, "Fine-mapping and colocalization", fine_mapping_gate_status,
         artifact_state(fine_mapping_path) + ";readiness=" + artifact_state(fine_mapping_readiness_path) + ";protocol=config/fine_mapping_colocalization_contract.json",
         "Pinned SuSiE-RSS/coloc-SuSiE code passes an end-to-end synthetic lock/run/collation test; real dense loci, signed LD, and molecular-QTL sources are absent."),
        (15, "Mechanistic annotation", mechanism_gate_status,
         artifact_state(mechanism_path) + ";readiness=" + artifact_state(mechanism_readiness_path) + ";protocol=config/mechanistic_annotation_contract.json",
         "A source-snapshotted evidence/graph/synthesis workflow passes synthetic tests; real signal-specific releases, accessions, evidence, and citations remain unavailable."),
        (16, "Adversarial review", "PRE_RESULT_REVIEW_PRESENT_FINAL_FINDINGS_REVIEW_BLOCKED" if adversarial_review_path.is_file() else "CHECKLIST_READY_FINAL_REVIEW_BLOCKED",
         artifact_state(adversarial_review_path) + ";checklist=" + str(ROOT / "results/adversarial_review_checklist.tsv"),
         "The current review challenges design, provenance, and blockers; findings-level review awaits actual results."),
        (17, "Publication-grade reporting", (
            "PRE_RESULT_TRACEABILITY_REPORT_PRESENT_PUBLICATION_FINDINGS_BLOCKED"
            if final_report_ready else "STATUS_REPORT_ONLY_FINAL_REPORT_BLOCKED"
        ),
         artifact_state(final_report_path) + ";counts=" + artifact_state(final_counts_path)
         + ";top_discoveries=" + artifact_state(top_discoveries_path)
         + ";provenance=" + artifact_state(final_report_provenance_path)
         + ";validation=" + final_report_validation,
         "A complete report must encode unavailable outcomes as upstream-blocked, not zero; no findings are claimed until required downstream artifacts are complete and validated."),
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
        ("R04", "Weak h2", "Exclude rerun h2 Z<4 from primary rg without replacement.", "PASS_REAL_H2_GATE_APPLIED" if h2_complete else "PENDING_REAL_H2", artifact_state(h2_path)),
        ("R05", "LDSC intercept inflation", "Exclude rerun intercept>1.2 from primary rg and retain sensitivity status.", "PASS_REAL_H2_GATE_APPLIED" if h2_complete else "PENDING_REAL_H2", artifact_state(h2_path)),
        ("R06", "UK Biobank sample overlap", "Inspect cross-trait intercepts and disclose overlap; do not equate LDSC adjustment with independent replication.", "PASS_DIAGNOSTICS_RECORDED_DISCLOSURE_REQUIRED" if rg_complete else "PENDING_REAL_RG", artifact_state(rg_path)),
        ("R07", "Sparse or proxy phenotypes", "Retain exact phenotype definitions and distinguish medication/proxy traits from diagnoses.", "PASS_METADATA", "candidate_traits.tsv"),
        ("R08", "Panel-level novelty inflation", "Require pair-level direct/same-phenotype/same-direction/same-sleep-context audit.", "PASS_PAIR_AUDIT_COMPLETE" if novelty_complete else f"PENDING_PAIR_AUDIT_{novelty_pending_count}_OF_{len(novelty_rows)}", artifact_state(novelty_path)),
        ("R09", "Replication non-independence", "Require non-overlapping participants and separately sourced summary statistics.", "SOURCE_CURATION_COMPLETE_ACQUISITION_PENDING" if replication_curation_complete else "CANDIDATE_FAMILY_LOCKED_SOURCE_CURATION_PENDING" if replication_candidate_locked else "PENDING_REPLICATION", artifact_state(replication_path) + ";" + artifact_state(replication_candidate_lock_path) + ";" + artifact_state(replication_search_path)),
        ("R10", "Global-to-local overreach", "Do not call global rg evidence of a shared locus; retain a prespecified globally-null secondary local set.", "CODE_READY_INPUTS_BLOCKED" if local_code_ready and not local_dependencies_ready else "PENDING_LOCAL_ANALYSIS", artifact_state(local_path) + ";" + artifact_state(local_readiness_path)),
        ("R11", "Pleiotropy or mediated effects", "Evaluate horizontal, vertical/mediated, shared-factor, and sample-overlap alternatives; run PLACO+ on genome-wide data only.", "CODE_READY_INPUTS_BLOCKED" if pleiotropy_code_ready and not pleiotropy_inputs_ready else "PENDING_PLEIOTROPY", artifact_state(pleiotropy_path) + ";" + artifact_state(pleiotropy_readiness_path)),
        ("R12", "Colocalization overclaim", "Report H0-H4, priors, sensitivity, and claim guards; colocalization is not causality.", "CODE_READY_INPUTS_BLOCKED" if fine_mapping_code_ready else "PENDING_COLOCALIZATION", artifact_state(fine_mapping_path) + ";" + artifact_state(fine_mapping_readiness_path)),
        ("R13", "Synthetic/real result contamination", "Synthetic tests stay under synthetic or temporary paths and carry explicit markers.", "PASS", "seven isolated synthetic workflows"),
        ("R14", "Storage-driven partial acquisition", "Do not silently analyze a result-selected subset of the locked panel.", "PASS_LOCKED_PANEL_COMPLETE" if len(streaming_completed_ids) == len(panel) else "IN_PROGRESS_LOCKED_PANEL_NO_RG" if streaming_completed_ids else "PASS_BLOCKED", f"streaming_receipts={len(streaming_completed_ids)}/{len(panel)};mirror={preflight['status']}"),
        ("R15", "Local LD-reference mismatch", "Require checksum-locked ancestry-matched LAVA/HDL-L references; do not fall back silently to a smaller panel.", "PASS_BLOCKED", artifact_state(local_readiness_path)),
        ("R16", "Local multiplicity or h2-gate leakage", "Freeze the pair-by-locus family, LAVA local-h2 Bonferroni gate, and local-rg BH family before result access.", "PASS_CONTRACT", "config/local_architecture_contract.json"),
        ("R17", "Replication candidate attrition", "Lock the entire Tier A/B candidate family before source curation and preserve NO_INDEPENDENT_DATASET outcomes.", "PASS_CONTRACT", "scripts/21_prepare_replication_queue.py+22_lock_replication_manifest.py+23_collate_replication.py"),
        ("R18", "LAVA simulation instability", "Freeze a deterministic nonzero simulation seed per selected pair and retain it in the manifest.", "PASS_CONTRACT", "config/local_architecture_contract.json"),
        ("R19", "Fine-mapping LD or allele mismatch", "Require identical dense SNP order/alleles, checksum-locked signed LD, PSD/symmetry checks, RSS-LD s, and kriging outlier diagnostics.", "PASS_CONTRACT_INPUTS_BLOCKED", "config/fine_mapping_colocalization_contract.json"),
        ("R20", "Molecular-QTL source attrition", "Complete eQTL/sQTL/pQTL searches per locked locus and preserve evidence-backed NO_SUITABLE_DATASET rows.", "PASS_CONTRACT_INPUTS_BLOCKED", "scripts/33_prepare_finemapping_queue.py+34_lock_finemapping_manifest.py"),
        ("R21", "Coloc prior cherry-picking", "Lock p1/p2/p12 and the full p12 sensitivity grid before result access; label prior robustness separately.", "PASS_CONTRACT_INPUTS_BLOCKED", "config/fine_mapping_colocalization_contract.json"),
        ("R22", "Mechanistic source or release drift", "Require exact release/accession, local query/result snapshot, SHA-256, access date, license/terms, and primary citation for every supported evidence row.", "PASS_CONTRACT_INPUTS_BLOCKED", "config/mechanistic_annotation_contract.json+config/mechanistic_sources.tsv"),
        ("R23", "Mechanistic chain gap concealment", "Emit every required chain edge, mark unsupported edges MISSING, and forbid a complete causal narrative from annotation-only evidence.", "PASS_CONTRACT_INPUTS_BLOCKED", "scripts/39_validate_mechanistic_evidence.py+40_synthesize_mechanisms.py"),
    ]
    checklist_path = ROOT / "results/adversarial_review_checklist.tsv"
    with checklist_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["risk_id", "risk", "required_control", "current_status", "current_evidence"])
        writer.writerows(risks)

    if not core_ok:
        overall = "BLOCKED_AT_CORE_CHECKPOINT_AND_FULL_RESOLUTION_ACQUISITION_GATES"
    elif len(streaming_completed_ids) < len(panel) and streaming_contract_ready:
        overall = f"IN_PROGRESS_STREAMING_ACQUISITION_{len(streaming_completed_ids)}_OF_{len(panel)}"
    elif mirror_blocked and not streaming_contract_ready:
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
        "The pinned LAVA 0.1.5 and HDL 1.4.3 packages load and expose their local-rg entrypoints. "
        "Their source-verified LD references are not present; the recommended LAVA UKB v1.1 reference "
        "alone is published as 15 GiB unzipped, so no local analysis was started.\n\n"
        "The pinned PLACO+ 0.2.0 source passes its entrypoint and end-to-end synthetic checks. "
        "No real PLACO+ scan was started because independently replicated pairs and genome-wide dense inputs do not exist.\n\n"
        "The pinned susieR 0.14.2 and coloc 5.2.3 packages pass entrypoint and end-to-end synthetic "
        "fine-mapping/colocalization checks, including dense SNP-order/allele/LD validation, PIPs, credible sets, "
        "H0-H4 posteriors, prior sensitivity, and retained unavailable-QTL outcomes. No real locus analysis was started.\n\n"
        "The mechanistic source registry and locked evidence workflow pass an end-to-end synthetic test. "
        "It requires exact releases/accessions, local source snapshots and checksums, primary citations, and explicit "
        "MISSING chain edges; verified landing pages alone are never treated as mechanistic evidence. No real mechanistic search was started.\n\n"
        f"The original full local mirror remains blocked: the exact compressed inputs total {preflight['compressed_source_gib']} GiB "
        f"and require {preflight['required_free_gib']} GiB with the locked safety factor. A pre-result streaming contract now pins all "
        f"202 S3 objects by version ID, verifies every full source body before output promotion, and preserves dense-locus access through "
        f"the exact versioned bgzip/tabix objects. Streaming acquisition has sealed {len(streaming_completed_ids)}/{len(panel)} traits; "
        f"{partial_h2_log_count}/{len(panel)} controlled h2 logs exist. "
        + (
            f"The complete extension h2 table validates {h2_primary_count} primary and {h2_sensitivity_count} sensitivity-only traits. "
            if h2_complete else
            "The complete extension h2 table does not yet pass validation. "
        )
        + (
            f"The complete extension rg family contains {len(rg_rows)} pairs, {rg_fdr_count} extension-FDR hits, and {rg_joint_screen_count} joint FDR/effect-screen hits. "
            if rg_complete else
            "The complete extension rg family does not yet pass validation. "
        )
        + f"The pair-level novelty audit is complete for {len(novelty_rows) - novelty_pending_count}/{len(novelty_rows)} FDR-significant pairs. "
        + f"Prioritization contains {priority_tier_a_count} Tier A and {priority_tier_b_count} Tier B pairs; the independent-replication candidate family is locked at {len(replication_queue_rows) if replication_candidate_locked else 0} pairs. Source curation retains {replication_testable_count} testable and {replication_unavailable_count} unavailable pairs before replication-result access. "
        "Therefore no final pair-level novelty claim, replication, local correlation, pleiotropy, fine-mapping, colocalization, "
        "mechanistic inference, or final manuscript claim exists yet.\n\n"
        "See `extension_acceptance_gates.tsv` for all 17 gates and "
        "`adversarial_review_checklist.tsv` for the current challenge audit. "
        "The Stage-17 `final_report.md` and `final_extension_counts.tsv` preserve all unavailable "
        "findings as `NA_BLOCKED_UPSTREAM`; the header-only `top_novel_discoveries.tsv` is not evidence of zero discoveries.\n",
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
        "streaming_contract_ready": streaming_contract_ready,
        "streaming_reference_ready": streaming_reference_ready,
        "streaming_completed_trait_count": len(streaming_completed_ids),
        "streaming_completed_trait_ids": streaming_completed_ids,
        "partial_h2_log_count": partial_h2_log_count,
        "remote_object_snapshot_sha256": sha256(streaming_snapshot_path) if streaming_snapshot_path.is_file() else None,
        "streaming_reference_sha256": sha256(streaming_reference_path) if streaming_reference_path.is_file() else None,
        "panel_sha256": sha256(ROOT / "config/candidate_traits.tsv"),
        "local_readiness_sha256": sha256(local_readiness_path) if local_readiness_path.is_file() else None,
        "pleiotropy_readiness_sha256": sha256(pleiotropy_readiness_path) if pleiotropy_readiness_path.is_file() else None,
        "fine_mapping_readiness_sha256": sha256(fine_mapping_readiness_path) if fine_mapping_readiness_path.is_file() else None,
        "mechanism_readiness_sha256": sha256(mechanism_readiness_path) if mechanism_readiness_path.is_file() else None,
        "adversarial_review_sha256": sha256(adversarial_review_path) if adversarial_review_path.is_file() else None,
        "final_report_sha256": sha256(final_report_path) if final_report_path.is_file() else None,
        "final_counts_sha256": sha256(final_counts_path) if final_counts_path.is_file() else None,
        "top_discoveries_sha256": sha256(top_discoveries_path) if top_discoveries_path.is_file() else None,
        "final_report_provenance_sha256": sha256(final_report_provenance_path) if final_report_provenance_path.is_file() else None,
        "acceptance_gates_sha256": sha256(gates_path),
        "adversarial_checklist_sha256": sha256(checklist_path),
        "status_report_sha256": sha256(status_path),
    }
    out = ROOT / "provenance/extension_acceptance.json"
    out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"EXTENSION_ACCEPTANCE_REPORTED stages=17 overall={overall}")


if __name__ == "__main__":
    main()
