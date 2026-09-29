#!/usr/bin/env python3
"""Validate the provisional Brain6 report against its locked evidence."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

from validate_pgc_mdd2025_rg_sensitivity import validate as validate_pgc_mdd2025_rg_sensitivity

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "brain6/FINAL_BRAIN6_REPORT.md"
PROVENANCE = ROOT / "brain6/FINAL_BRAIN6_REPORT.provenance.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def report_significant_rows(report: str) -> list[tuple[str, ...]]:
    section = report.split("Complete inherited significant set", 1)
    if len(section) != 2:
        raise ValueError("Report is missing the complete inherited significant table")
    section = section[1].split("## Replication", 1)[0]
    rows = []
    for line in section.splitlines():
        if not line.startswith("|") or line.startswith("|---"):
            continue
        cells = tuple(cell.strip() for cell in line.strip().strip("|").split("|"))
        if cells and cells[0] != "Sleep trait ID":
            rows.append(cells)
    if any(len(row) != 6 for row in rows):
        raise ValueError("Significant-result table has a malformed row")
    return rows


def validate_report_table(report: str, global_rows: list[dict[str, str]]) -> int:
    significant = sorted(
        (row for row in global_rows if row["significance_under_original_396_family"].lower() == "true"),
        key=lambda row: float(row["original_BH_FDR_q"]),
    )
    fmt = lambda value: f"{float(value):.3g}"
    expected = [
        (row["sleep_trait"], row["brain_disorder_label"], fmt(row["rg"]),
         fmt(row["se"]), fmt(row["p"]), fmt(row["original_BH_FDR_q"]))
        for row in significant
    ]
    observed = report_significant_rows(report)
    if len(global_rows) != 72 or len(significant) != 35:
        raise ValueError("Locked global map no longer has the report's 72-row/35-hit dimensions")
    if observed != expected:
        raise ValueError("Report significant-result table differs from the locked global map")
    if any(row["brain_disorder"] == "alz" and
           row["significance_under_original_396_family"].lower() == "true"
           for row in global_rows):
        raise ValueError("Report's no-Alzheimer-hit statement conflicts with the locked map")
    return len(observed)


def validate_report(root: Path = ROOT) -> dict[str, Any]:
    report_path = root / "brain6/FINAL_BRAIN6_REPORT.md"
    provenance_path = root / "brain6/FINAL_BRAIN6_REPORT.provenance.json"
    report = report_path.read_text(encoding="utf-8")
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("status") != "PROVISIONAL_REPORT_NOT_FROZEN":
        raise ValueError("The report must remain provisional while final-family work is pending")
    if provenance.get("report_path") != str(report_path.relative_to(root)):
        raise ValueError("Report provenance points to a different report path")
    if sha256(report_path) != provenance.get("report_sha256"):
        raise ValueError("Report checksum differs from its provenance")
    inventory = root / provenance.get("latest_inventory_path", "")
    if sha256(inventory) != provenance.get("inventory_sha256"):
        raise ValueError("Project-state inventory checksum differs from report provenance")
    for source, expected_hash in provenance.get("source_hashes", {}).items():
        path = root / source
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"Report source is missing or changed: {source}")

    latest_inventory_path = root / provenance.get("latest_inventory_path", "")
    latest_inventory = json.loads(latest_inventory_path.read_text(encoding="utf-8"))
    if (provenance.get("latest_inventory_path") != str(latest_inventory_path.relative_to(root)) or
            provenance.get("latest_inventory_sha256") != sha256(latest_inventory_path) or
            latest_inventory.get("inventory_id") != "brain6-only-project-state-20260925" or
            len(latest_inventory.get("stages", [])) != 12):
        raise ValueError("Refreshed project-state inventory does not match its report provenance")
    evidence_names = set(latest_inventory.get("evidence_files", {}))
    if any(name not in evidence_names for stage in latest_inventory.get("stages", [])
           for name in stage.get("evidence", [])):
        raise ValueError("Project-state inventory contains a dangling stage evidence reference")
    for evidence in latest_inventory.get("evidence_files", {}).values():
        path = Path(evidence["path"])
        if not path.is_file() or path.stat().st_size != evidence["bytes"] or sha256(path) != evidence["sha256"]:
            raise ValueError(f"Project-state inventory evidence changed: {evidence['path']}")
    stage_by_id = {stage["stage"]: stage for stage in latest_inventory["stages"]}
    screen_stage = stage_by_id.get("sleep_power_optimized_screen", {})
    screen_metrics = screen_stage.get("metrics", {})
    screen_state_valid = (
        (screen_stage.get("status") == "RUNNING" and
         screen_metrics.get("candidate_screen_workers_configured") == 4 and
         screen_metrics.get("running_workers") == 4 and
         screen_metrics.get("worker_liveness") == "VERIFIED_RUNNING") or
        (screen_stage.get("status") == "PARTIAL" and
         screen_metrics.get("candidate_screen_workers_configured") == 4 and
         screen_metrics.get("worker_liveness") == "NOT_VERIFIED")
    )
    if (stage_by_id.get("canonical_lava_v3", {}).get("status") != "COMPLETE_FAILED_QC" or
            stage_by_id.get("protected_track_b", {}).get("status") != "BLOCKED_EXTERNAL" or
            not screen_state_valid or
            screen_metrics.get("technical_advancement_trigger_met") is not True or
            screen_metrics.get("substitution_family_not_run") != 3305 or
            screen_metrics.get("substitution_family_cells") != 17465 or
            screen_metrics.get("substitution_family_maximum_not_run_cells") != 873 or
            screen_metrics.get("substitution_family_passes_frozen_5_percent_qc") is not False or
            screen_metrics.get("candidate_specific_overlap_matrices") != 3 or
            screen_metrics.get("candidate_bivariate_lava_started") is not False or
            screen_metrics.get("sleep_replacement_family_feasibility_status") != "PASS_SLEEP_ONLY_REPLACEMENTS_CANNOT_PASS_FROZEN_FAMILY_QC" or
            screen_metrics.get("candidate_decision_status") != "PASS_TECHNICAL_ADVANCEMENT_FAMILY_QC_BLOCKED" or
            screen_metrics.get("candidate_technical_advancement_count") != 1 or
            screen_metrics.get("candidate_full_family_qc_passing_count") != 0 or
            screen_metrics.get("best_case_one_longsleep_replacement_not_run") != 2429 or
            screen_metrics.get("best_case_one_longsleep_replacement_excess") != 1556 or
            screen_metrics.get("best_case_both_sleep_replacements_not_run") != 1758 or
            screen_metrics.get("best_case_both_sleep_replacements_excess") != 885):
        raise ValueError("Brain6-only inventory has stale or misclassified Brain6 stages")
    annotation_stage = stage_by_id.get("partial_candidate_gene_context", {})
    annotation_metrics = annotation_stage.get("metrics", {})
    if (annotation_stage.get("status") != "COMPLETE_DESCRIPTIVE_DIAGNOSTIC" or
            annotation_metrics.get("candidate_regions_gated") != 19 or
            annotation_metrics.get("regions_not_tiered") != 19 or
            annotation_metrics.get("lead_signals_annotated") != 21 or
            annotation_metrics.get("cross_pair_overlapping_region_pairs") != 3 or
            annotation_metrics.get("exact_shared_lead_variant_ids") != ["rs1051168"] or
            annotation_metrics.get("lead_signals_stable_across_ld_methods") != 21 or
            annotation_metrics.get("range_filtered_method_specific_leads") != 4 or
            annotation_metrics.get("method_specific_lead_ids") != ["rs199534", "rs199535", "rs4938023", "rs6430538"] or
            annotation_metrics.get("candidate_placo_pvalues_recomputed") != 2297 or
            annotation_metrics.get("candidate_placo_pvalues_matching") != 2297 or
            annotation_metrics.get("candidate_placo_pvalue_max_relative_error", 1) >= 1e-12 or
            annotation_metrics.get("causal_gene_claims") is not False or
            annotation_metrics.get("fine_mapping_performed") is not False or
            annotation_metrics.get("colocalization_performed") is not False or
            "Tables S36–S39" not in report or
            "2,297 candidate p-values" not in report or
            any(lead not in report for lead in ("rs6430538", "rs199535", "rs199534", "rs4938023")) or
            "rs1051168" not in report or "NOT_TIERED" not in report):
        raise ValueError("Report or inventory overstates partial candidate gene context")
    placo_replay_path = root / "brain6/results/loci/placo_candidate_pvalue_validation_v1.tsv"
    placo_replay_provenance_path = root / "brain6/results/loci/placo_candidate_pvalue_validation_v1.provenance.json"
    placo_replay_summary_path = root / "brain6/results/supplement/table_S39_partial_placo_pvalue_replay.tsv"
    placo_replay = read_tsv(placo_replay_path)
    placo_replay_provenance = json.loads(placo_replay_provenance_path.read_text(encoding="utf-8"))
    placo_replay_summary = read_tsv(placo_replay_summary_path)
    if (placo_replay_provenance.get("status") != "PASS_PINNED_SOURCE_CANDIDATE_PVALUES_RECOMPUTED" or
            placo_replay_provenance.get("output", {}).get("sha256") != sha256(placo_replay_path) or
            placo_replay_provenance.get("summary_output", {}).get("sha256") != sha256(placo_replay_summary_path) or
            placo_replay_provenance.get("output", {}).get("rows") != 2297 or
            len(placo_replay) != 2297 or len(placo_replay_summary) != 4 or
            any(row.get("validation_status") != "PASS" for row in placo_replay) or
            sum(int(row["candidate_variants"]) for row in placo_replay_summary) != 2297 or
            sum(int(row["matching_p_values"]) for row in placo_replay_summary) != 2297 or
            any(float(row["maximum_relative_error"]) >= 1e-8 for row in placo_replay_summary) or
            "2,297 candidate p-values" not in report or "Table S39" not in report):
        raise ValueError("Candidate PLACO p-value replay is missing, stale, or overinterpreted")
    if any("frailty" in stage["stage"].lower() or "scheduler_state" in stage.get("evidence", [])
           for stage in latest_inventory["stages"]):
        raise ValueError("Brain6 inventory contains cross-project FI/frailty scheduler scope")
    if latest_inventory_path.name not in report:
        raise ValueError("Report omits the refreshed project-state inventory")
    if "even perfect coverage for both sleep traits leaves 885 NOT_RUN cells above the frozen family limit" not in report:
        raise ValueError("Report omits the sleep-replacement family feasibility lower bound")

    global_path = root / "brain6/results/global/brain6_72_locked.tsv"
    global_rows = read_tsv(global_path)
    significant_count = validate_report_table(report, global_rows)

    novelty_path = root / "brain6/results/novelty/brain6_significant_pair_novelty_crosswalk.tsv"
    novelty_provenance_path = novelty_path.with_suffix(".provenance.json")
    novelty_rows = read_tsv(novelty_path)
    novelty_provenance = json.loads(novelty_provenance_path.read_text(encoding="utf-8"))
    parent_audit_path = root / "results/analysis/literature_novelty_audit.tsv"
    parent_rows = read_tsv(parent_audit_path)
    parent_by_pair = {(row["sleep_trait"], row["external_trait"]): row for row in parent_rows}
    expected_novelty = {
        (row["sleep_trait"], row["brain_disorder"])
        for row in global_rows
        if row["significance_under_original_396_family"].lower() == "true"
    }
    observed_novelty = {(row["sleep_trait"], row["brain_disorder"]) for row in novelty_rows}
    classes = Counter(row["novelty_classification"] for row in novelty_rows)
    novelty_inputs = novelty_provenance.get("inputs", {})
    if (observed_novelty != expected_novelty or len(novelty_rows) != 35 or
            any(key not in parent_by_pair for key in observed_novelty) or
            novelty_provenance.get("status") != "PASS_DATED_LITERATURE_CROSSWALK_NO_NOVELTY_CLAIM" or
            novelty_provenance.get("output", {}).get("sha256") != sha256(novelty_path) or
            novelty_provenance.get("summary", {}).get("parent_audit_search_as_of_date") != "2026-08-28" or
            classes != Counter({"DIRECT_RG_PREVIOUSLY_REPORTED": 25,
                                "DIRECT_RG_REPLICATION_DIFFERENT_DATASET": 9, "MR_ONLY": 1}) or
            any(not (root / path).is_file() or sha256(root / path) != digest
                for path, digest in novelty_inputs.items())):
        raise ValueError("Brain6 novelty crosswalk is incomplete or not source-bound")
    for row in novelty_rows:
        source = parent_by_pair[(row["sleep_trait"], row["brain_disorder"])]
        if any(row[field] != source[field] for field in (
                "novelty_classification", "representative_pmid", "representative_doi",
                "phenotype_definition_match", "dataset_relation", "comparison_independence",
                "search_as_of_date")):
            raise ValueError("Brain6 novelty crosswalk differs from the parent audit")
    if any(value not in report for value in (
            "2026-08-28", "25 are classified as previously reported direct rg",
            "9 as direct rg reported in another dataset", "none is classified as apparently novel",
            "Table S29")):
        raise ValueError("Report omits the dated, qualified novelty crosswalk")

    pgc_sensitivity = validate_pgc_mdd2025_rg_sensitivity()
    if (pgc_sensitivity.get("status") != "PASS_PGC_MDD2025_RG_SENSITIVITY" or
            any(value not in report for value in ("0.4771", "0.0228", "4.65×10⁻⁹⁷", "Table S34",
                                                   "nine cohort entries", "not independent replication"))):
        raise ValueError("Report omits or misstates the qualified PGC MDD2025 global sensitivity")

    replication = read_tsv(root / "brain6/results/replication/replication_master.tsv")
    adhd = next(row for row in replication if row["brain_disorder"] == "adhd")
    expected_replication = (adhd["replication_rg"], adhd["replication_SE"], adhd["replication_P"])
    if expected_replication != ("0.3817", "0.0521", "2.32742962810806e-13"):
        raise ValueError("The archived insomnia–ADHD replication row changed")
    for value in ("0.3817", "0.0521", "2.327×10⁻¹³"):
        if value not in report:
            raise ValueError("Report omits a checked insomnia–ADHD replication value")

    decision = json.loads((root / (
        "work/lava-canonical-v3-production/"
        "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/"
        "canonical_family_decision.json")).read_text(encoding="utf-8"))
    if (decision.get("overall_status") != "FAILED_QC_NOT_PROMOTED" or
            decision.get("status_counts", {}).get("NOT_RUN") != 3720):
        raise ValueError("Canonical LAVA decision no longer matches the report")
    if "3,720" not in report or "FAILED_QC_NOT_PROMOTED" not in report:
        raise ValueError("Report omits the canonical LAVA frozen-QC failure")
    cause_provenance_path = root / "brain6/results/lava/canonical_v3_cause_breakdown_v2.provenance.json"
    cause_provenance = json.loads(cause_provenance_path.read_text(encoding="utf-8"))
    trait_summary_path = root / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv"
    locus_summary_path = root / "brain6/results/lava/canonical_v3_status_by_locus_v1.tsv"
    trait_summary = read_tsv(trait_summary_path)
    locus_summary = read_tsv(locus_summary_path)
    source_hashes = cause_provenance.get("source_hashes", {})
    cause_sources = {
        "canonical_family_decision": root / (
            "work/lava-canonical-v3-production/"
            "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/"
            "canonical_family_decision.json"),
        "canonical_aggregate": root / (
            "work/lava-canonical-v3-production/"
            "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/"
            "results/canonical_family_results.tsv"),
        "canonical_latest_audit": root / (
            "work/lava-canonical-v3-production/"
            "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b/"
            "latest_audit.json"),
        "cell_level_not_run_table": root / "brain6/results/lava/canonical_v3_not_run_cells.tsv",
        "cell_level_not_run_provenance": root / "brain6/results/lava/canonical_v3_not_run_cells.provenance.json",
    }
    if set(source_hashes) != set(cause_sources):
        raise ValueError("Canonical v3 cause summary has an incomplete source-hash map")
    if (cause_provenance.get("builder_sha256") != sha256(root / cause_provenance.get("builder_path", "")) or
            any(sha256(cause_sources[name]) != expected for name, expected in source_hashes.items()) or
            cause_provenance.get("status") != "COMPLETE_DIAGNOSTIC_ONLY" or
            cause_provenance.get("summary", {}).get("not_run_cells") != 3720 or
            cause_provenance.get("summary", {}).get("not_run_loci") != 2020 or
            cause_provenance.get("summary", {}).get("loci_with_at_least_five_not_run") != 41 or
            cause_provenance.get("summary", {}).get("max_not_run_cells_at_one_locus") != 7 or
            cause_provenance.get("trait_table", {}).get("sha256") != sha256(trait_summary_path) or
            cause_provenance.get("locus_table", {}).get("sha256") != sha256(locus_summary_path) or
            len(trait_summary) != 7 or len(locus_summary) != 2495 or
            sum(int(row["not_run_cells"]) for row in trait_summary) != 3720 or
            sum(int(row["planned_cells"]) for row in locus_summary) != 17465):
        raise ValueError("Canonical v3 trait/locus cause summaries are incomplete or hash-invalid")
    if any(value not in report for value in ("Table S25", "Table S26", "2,020 loci", "long sleep (1,291/2,495; 51.74%)")):
        raise ValueError("Report omits the validated canonical v3 trait/locus cause breakdown")

    candidate_loci = read_tsv(root / "brain6/results/loci/placo_candidate_loci_partial.tsv")
    candidate_variants = read_tsv(root / "brain6/results/loci/placo_candidate_variants_partial.tsv")
    if len(candidate_loci) != 19 or len(candidate_variants) != 2297:
        raise ValueError("Partial PLACO candidate counts differ from their machine-readable tables")
    if "2,297" not in report or "19 UKB-reference intervals" not in report:
        raise ValueError("Report omits the partial-family PLACO candidate scope/counts")

    direction_path = root / "brain6/results/loci/placo_candidate_directions_partial.tsv"
    direction_provenance_path = direction_path.with_suffix(".provenance.json")
    direction_provenance = json.loads(direction_provenance_path.read_text(encoding="utf-8"))
    direction_rows = read_tsv(direction_path)
    if (direction_provenance.get("status") != "PASS_PARTIAL_FAMILY_DESCRIPTIVE" or
            direction_provenance.get("scope", {}).get("complete_five_track_family") is not False or
            direction_provenance.get("output", {}).get("sha256") != sha256(direction_path) or
            len(direction_rows) != 2297 or
            direction_provenance.get("direction_counts", {}).get("all") !=
            {"CONCORDANT": 2296, "OPPOSING": 1}):
        raise ValueError("Partial PLACO direction summary differs from its source-bound artifact")
    if "2,296" not in report or "rs4790841" not in report or "not signed-LD analysis" not in report:
        raise ValueError("Report omits the qualified partial PLACO direction description")

    valid_ld_provenance_path = root / "brain6/results/loci/placo_candidate_valid_ld_partial_v2.provenance.json"
    valid_ld_provenance = json.loads(valid_ld_provenance_path.read_text(encoding="utf-8"))
    valid_ld_summary = read_tsv(root / "brain6/results/supplement/table_S22_partial_placo_valid_ld_rebuild.tsv")
    if (valid_ld_provenance.get("status") != "PASS_PARTIAL_WITH_LD_RANGE_EXCEPTIONS" or
            valid_ld_provenance.get("scope", {}).get("out_of_range_ld_edges") != 12317 or
            valid_ld_provenance.get("scope", {}).get("pair_candidate_exception_rows") != 12321 or
            sum(int(row["original_lead_variants"]) for row in valid_ld_summary) != 21 or
            sum(int(row["rebuild_valid_ld_lead_variants"]) for row in valid_ld_summary) != 25 or
            sum(int(row["rebuild_intervals"]) for row in valid_ld_summary) != 19):
        raise ValueError("Range-filtered partial LD rebuild differs from its provenance/summary table")
    if any(value not in report for value in ("12,317", "1,526,086", "25 lead signals", "Table S22")):
        raise ValueError("Report omits the qualified partial LD range-QC finding")

    signed_effect_path = root / "brain6/results/loci/placo_candidate_signed_effect_ld_partial_v1.tsv"
    signed_effect_provenance_path = signed_effect_path.with_suffix(".provenance.json")
    signed_effect_provenance = json.loads(signed_effect_provenance_path.read_text(encoding="utf-8"))
    signed_effect_rows = read_tsv(signed_effect_path)
    signed_summary = read_tsv(root / "brain6/results/supplement/table_S23_partial_placo_signed_effect_ld.tsv")
    matched_effects = [row for row in signed_effect_rows if row["ld_qc_status"] == "PASS_VALID_SIGNED_R"]
    if (signed_effect_provenance.get("status") != "PASS_PARTIAL_FAMILY_DESCRIPTIVE_WITH_LD_RANGE_EXCEPTIONS" or
            signed_effect_provenance.get("scope", {}).get("complete_five_track_family") is not False or
            len(signed_effect_rows) != 2297 or len(matched_effects) != 2273 or
            any(row["sleep_effect_vs_lead_ld"] != "CONSISTENT_WITH_LD_SIGN" or
                row["disorder_effect_vs_lead_ld"] != "CONSISTENT_WITH_LD_SIGN" for row in matched_effects) or
            sum(int(row["both_effects_consistent_with_ld_sign"]) for row in signed_summary) != 2273):
        raise ValueError("Partial effect-versus-signed-LD summary differs from its source-bound artifact")
    if any(value not in report for value in ("2,273 exact-reference candidates", "Table S23", "not causal evidence")):
        raise ValueError("Report omits the qualified partial effect-versus-LD description")

    roundoff_path = root / provenance["roundoff_audit_path"]
    roundoff_audit = json.loads(roundoff_path.read_text(encoding="utf-8"))
    roundoff = roundoff_audit["lava_roundoff"]
    baseline = roundoff_audit["lava"]
    if (roundoff["run_id"] != "af43fde06c4c6c0d57a6ab104bea5aa33f24fb2ff2b89b53163529160452d0cf" or
            roundoff["receipt_errors"] != 0 or roundoff["complete"] is not True or
            roundoff["verified_loci"] != 2495 or roundoff["planned_loci"] != 2495 or
            baseline["verified_loci"] != 2495 or baseline["receipt_errors"] != 0 or
            baseline["manifest_completed_loci"] != 743 or baseline["manifest_matches_receipts"] is not False):
        raise ValueError("Roundoff audit is not a complete, receipt-valid family snapshot")
    if ("2,495/2,495" not in report or "12,475" not in report or "insomnia/2" not in report or
            "manifest counter still says 743" not in report):
        raise ValueError("Report roundoff completion/comparison/aggregation failure is not documented")
    manifest = json.loads((root / "brain6/manifests/lava_roundoff_v1_run.json").read_text(encoding="utf-8"))
    if (manifest.get("run_id") != roundoff["run_id"] or
            manifest.get("status") != "AGGREGATION_FAILED_COMPLETE_RECEIPTS" or
            manifest.get("family_decision_status") != "NOT_CREATED_AGGREGATION_FAILED"):
        raise ValueError("Roundoff run state changed; regenerate the report snapshot")
    comparison_path = root / "brain6/results/lava/lava_baseline_roundoff_pair_locus_comparison_v1.tsv"
    comparison_provenance_path = comparison_path.with_suffix(".provenance.json")
    comparison_provenance = json.loads(comparison_provenance_path.read_text(encoding="utf-8"))
    comparison_rows = read_tsv(comparison_path)
    summary_path = root / "brain6/results/supplement/table_S24_lava_roundoff_comparison.tsv"
    summary_rows = read_tsv(summary_path)
    if (comparison_provenance.get("status") != "DESCRIPTIVE_COMPLETE_FAMILY_QC_UNRESOLVED" or
            comparison_provenance.get("planned_pair_locus_slots") != 12475 or
            comparison_provenance.get("comparison_table_sha256") != sha256(comparison_path) or
            comparison_provenance.get("summary_table_sha256") != sha256(summary_path) or
            len(comparison_rows) != 12475 or len(summary_rows) != 5 or
            sum(int(row["loci"]) for row in summary_rows) != 12475 or
            comparison_provenance.get("both_tested_slots") != 1 or
            comparison_provenance.get("numeric_difference_summary", {}).get("max_abs_delta_p") != 0.0 or
            comparison_provenance.get("numeric_difference_summary", {}).get("max_abs_delta_local_rg") != 0.0):
        raise ValueError("Roundoff comparison is not a complete, hash-valid descriptive artifact")
    if "Table S24" not in report:
        raise ValueError("Report omits the roundoff comparison supplement table")

    context_path = root / "brain6/results/lava/lava_pair_context_univariates_v1.tsv"
    context_summary_path = root / "brain6/results/supplement/table_S27_lava_pair_context_univariates.tsv"
    context_provenance_path = context_path.with_suffix(".provenance.json")
    context_provenance = json.loads(context_provenance_path.read_text(encoding="utf-8"))
    context_summary = read_tsv(context_summary_path)
    if (context_provenance.get("status") != "DESCRIPTIVE_COMPLETE_NO_FAMILY_PROMOTION" or
            context_provenance.get("run_ids") != {
                "baseline": "13b80e6b64a5179fec70ba510a89c51241c9d58ec05833cc52161b435aee44fa",
                "roundoff": "af43fde06c4c6c0d57a6ab104bea5aa33f24fb2ff2b89b53163529160452d0cf",
            } or
            context_provenance.get("receipt_verified_loci") != {"baseline": 2495, "roundoff": 2495} or
            context_provenance.get("pair_context_rows") != 36819 or
            context_provenance.get("output_sha256") != sha256(context_path) or
            context_provenance.get("summary_sha256") != sha256(context_summary_path) or
            len(read_tsv(context_path)) != 36819 or len(context_summary) != 14):
        raise ValueError("Pair-context univariate audit is incomplete or hash-invalid")
    for label, expected_repeated, expected_conflict, expected_missing in (
        ("baseline", 2804, 2802, 2833), ("roundoff", 2806, 2804, 2832),
    ):
        selected = [row for row in context_summary if row["run"] == label]
        if (sum(int(row["repeated_keys"]) for row in selected) != expected_repeated or
                sum(int(row["repeated_keys_disagree_gt_1e-12"]) for row in selected) != expected_conflict or
                sum(int(row["missing_unique_trait_locus"]) for row in selected) != expected_missing):
            raise ValueError(f"Pair-context audit counts changed for {label}")
    if any(value not in report for value in ("Table S27", "2,802/2,804", "2,804/2,806", "14,632/17,465", "14,633/17,465")):
        raise ValueError("Report omits the source-verified pair-context discrepancy results")

    probe_path = root / "brain6/results/lava/lava_unique_gate_probe_v1.tsv"
    probe_provenance_path = probe_path.with_suffix(".provenance.json")
    probe_provenance = json.loads(probe_provenance_path.read_text(encoding="utf-8"))
    probe_rows = read_tsv(probe_path)
    if (probe_provenance.get("status") != "DIAGNOSTIC_ONLY_NOT_PROMOTED" or
            probe_provenance.get("canonical_v3_status") != "FAILED_QC_NOT_PROMOTED" or
            len(probe_rows) != 2 or
            len(probe_provenance.get("candidate_pair_loci", [])) != 2):
        raise ValueError("Pair-independent gate probe is missing or not explicitly diagnostic-only")
    for record in probe_provenance.get("files", []):
        path = Path(record["path"])
        if not path.is_file() or path.stat().st_size != record["bytes"] or sha256(path) != record["sha256"]:
            raise ValueError(f"Pair-independent gate probe provenance mismatch: {record['path']}")
    by_pair = {row["pair_id"]: row for row in probe_rows}
    if (set(by_pair) != {"longsleep__scz", "longsleep__parkinson"} or
            any(row["status"] != "DIAGNOSTIC_ONLY" or row["locus_id"] != "2207" for row in probe_rows) or
            by_pair["longsleep__scz"]["p"] != "0.120592" or
            by_pair["longsleep__scz"]["q_bh_full_12475"] != "1" or
            by_pair["longsleep__parkinson"]["p"] != "2.21035e-05" or
            abs(float(by_pair["longsleep__parkinson"]["q_bh_full_12475"]) - 0.2757411625) > 1e-12 or
            any(float(row["q_bh_full_12475"]) <= 0.05 for row in probe_rows)):
        raise ValueError("Pair-independent gate probe estimates or full-family correction changed")
    if any(value not in report for value in (
            "diagnostic-only probe", "12,475-slot BH q=1.000", "q=0.276",
            "lava_unique_gate_probe_v1.tsv")):
        raise ValueError("Final report omits qualified pair-independent gate probe results")

    return {
        "status": "PASS_PROVISIONAL_REPORT",
        "global_rows": len(global_rows),
        "inherited_significant_rows": significant_count,
        "roundoff_verified_loci": roundoff["verified_loci"],
        "roundoff_pair_locus_slots_compared": len(comparison_rows),
        "roundoff_comparison_summary_rows": len(summary_rows),
        "unique_gate_probe_rows": len(probe_rows),
        "current_inventory_stages": len(latest_inventory["stages"]),
        "partial_candidate_regions_not_tiered": annotation_metrics["regions_not_tiered"],
        "partial_candidate_lead_nearest_gene_rows": annotation_metrics["lead_signals_annotated"],
        "partial_candidate_cross_pair_region_overlaps": annotation_metrics["cross_pair_overlapping_region_pairs"],
        "partial_candidate_exact_shared_lead_variant_ids": annotation_metrics["exact_shared_lead_variant_ids"],
        "partial_candidate_leads_stable_across_ld_methods": annotation_metrics["lead_signals_stable_across_ld_methods"],
        "partial_candidate_range_filter_specific_leads": annotation_metrics["method_specific_lead_ids"],
        "canonical_v3_cause_summary_trait_rows": len(trait_summary),
        "canonical_v3_cause_summary_locus_rows": len(locus_summary),
        "partial_placo_valid_ld_status": valid_ld_provenance["status"],
        "partial_placo_valid_ld_exceptions": valid_ld_provenance["scope"]["pair_candidate_exception_rows"],
        "roundoff_receipt_errors": roundoff["receipt_errors"],
        "baseline_verified_loci": baseline["verified_loci"],
        "baseline_manifest_counter": baseline["manifest_completed_loci"],
        "baseline_manifest_matches_receipts": baseline["manifest_matches_receipts"],
        "report_sha256": sha256(report_path),
        "provenance_sha256": sha256(provenance_path),
    }


def main() -> None:
    print(json.dumps(validate_report(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
