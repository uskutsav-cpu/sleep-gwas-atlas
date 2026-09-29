"""Write an immutable, hash-bound snapshot of current Brain6 stage status."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
ROUND_AUDIT = ROOT / "brain6/qc/lava_production_checkpoint_audit_20260924T2252Z.json"
V3_RUN = ROOT / (
    "work/lava-canonical-v3-production/"
    "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite project-state inventory: {output}")

    screen_root = ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1"
    worker_configs = sorted(screen_root.glob("worker_*.config.json"))
    worker_config_data = [read_json(path) for path in worker_configs]
    worker_result_paths = [ROOT / item["output_tsv"] for item in worker_config_data]
    worker_summary_paths = [ROOT / item["summary_json"] for item in worker_config_data]
    durable_worker_results = [path for path in worker_result_paths if path.is_file()]
    durable_worker_summaries = [path for path in worker_summary_paths if path.is_file()]
    aggregate_config = read_json(screen_root / "config.json")
    aggregate_tsv = ROOT / aggregate_config["output_tsv"]
    aggregate_summary_path = ROOT / aggregate_config["summary_json"]
    aggregate_provenance_path = screen_root / "aggregate.provenance.json"
    aggregate_summary: dict[str, Any] | None = None
    if any(path.exists() for path in
           (aggregate_tsv, aggregate_summary_path, aggregate_provenance_path)):
        if not all(path.is_file() for path in
                   (aggregate_tsv, aggregate_summary_path, aggregate_provenance_path)):
            raise ValueError("Candidate LAVA aggregate is only partially durable")
        aggregate_provenance = read_json(aggregate_provenance_path)
        aggregate_script = ROOT / "brain6/scripts/aggregate_power_optimized_sleep_screen_v1.py"
        if (aggregate_provenance.get("status") != "COMPLETE_SCREEN_NO_FAILED_ROWS" or
                aggregate_provenance.get("base_config_sha256") != sha256(screen_root / "config.json") or
                aggregate_provenance.get("aggregator_sha256") != sha256(aggregate_script)):
            raise ValueError("Candidate LAVA aggregate provenance/configuration mismatch")
        # The aggregate is outcome-blind and binds upstream source/output hashes
        # directly. Do not create a provenance cycle by making this aggregate
        # sidecar hash the later top-level audit provenance that itself records
        # the aggregate outputs.
        for rel, expected_hash in aggregate_provenance.get("upstream_audit_inputs", {}).items():
            path = ROOT / rel
            if not path.is_file() or sha256(path) != expected_hash:
                raise ValueError(f"Candidate aggregate upstream audit input changed: {rel}")
        upstream_outputs = aggregate_provenance.get("upstream_audit_outputs", {})
        if set(upstream_outputs) != {"README.md", "candidate_replacement_screen.tsv", "current_sleep_trait_audit.tsv"}:
            raise ValueError("Candidate aggregate upstream output map is malformed")
        for name, expected_hash in upstream_outputs.items():
            path = ROOT / "brain6/results/power_optimized_sensitivity_v1" / name
            if not path.is_file() or sha256(path) != expected_hash:
                raise ValueError(f"Candidate aggregate upstream audit output changed: {rel}")
        for rel, expected_hash in aggregate_provenance.get("inputs", {}).items():
            path = ROOT / rel
            if not path.is_file() or sha256(path) != expected_hash:
                raise ValueError(f"Candidate LAVA aggregate input changed: {rel}")
        for rel, expected_hash in aggregate_provenance.get("outputs", {}).items():
            path = ROOT / rel
            if not path.is_file() or sha256(path) != expected_hash:
                raise ValueError(f"Candidate LAVA aggregate output changed: {rel}")
        aggregate_summary = read_json(aggregate_summary_path)
        if (aggregate_summary.get("status") != "COMPLETE_SCREEN_NO_FAILED_ROWS" or
                aggregate_summary.get("rows") != 2495 or aggregate_summary.get("failed") != 0 or
                aggregate_summary.get("canonical_decision_unchanged") != "FAILED_QC_NOT_PROMOTED"):
            raise ValueError("Candidate LAVA aggregate summary is incomplete or misclassified")
        aggregate_rows = read_tsv(aggregate_tsv)
        if len(aggregate_rows) != 2495:
            raise ValueError("Candidate LAVA aggregate row count is incomplete")
        status_counts = {status: sum(row["status"] == status for row in aggregate_rows)
                         for status in ("TESTED", "NOT_RUN", "FAILED")}
        reason_counts: dict[str, int] = {}
        for row in aggregate_rows:
            if row["status"] == "NOT_RUN":
                reason_counts[row["reason"]] = reason_counts.get(row["reason"], 0) + 1
        if (status_counts["TESTED"] != aggregate_summary.get("tested") or
                status_counts["NOT_RUN"] != aggregate_summary.get("not_run") or
                status_counts["FAILED"] != aggregate_summary.get("failed") or
                sum(reason_counts.values()) != status_counts["NOT_RUN"]):
            raise ValueError("Candidate LAVA aggregate rows disagree with its summary")
        aggregate_summary["reason_counts"] = reason_counts
        completion_provenance_path = ROOT / "brain6/results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.provenance.json"
        completion_provenance = read_json(completion_provenance_path)
        if completion_provenance.get("status") != "COMPLETE_SCREEN_NO_PROMOTION":
            raise ValueError("Candidate completion report has an invalid promotion status")
        completion_inputs = {
            "current_trait_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_sleep_trait_audit.tsv",
            "canonical_cause_summary": ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv",
            "upstream_power_audit_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/provenance.json",
            "candidate_aggregate": aggregate_tsv,
            "candidate_aggregate_summary": aggregate_summary_path,
            "candidate_aggregate_provenance": aggregate_provenance_path,
            "aggregate_builder": ROOT / "brain6/scripts/aggregate_power_optimized_sleep_screen_v1.py",
        }
        if set(completion_provenance.get("inputs", {})) != set(completion_inputs):
            raise ValueError("Candidate completion report input manifest is incomplete")
        for name, expected_hash in completion_provenance.get("inputs", {}).items():
            path = completion_inputs[name]
            if not path.is_file() or sha256(path) != expected_hash:
                raise ValueError(f"Candidate completion report input changed: {name}")
        for rel, expected_hash in completion_provenance.get("outputs", {}).items():
            path = ROOT / rel
            if not path.is_file() or sha256(path) != expected_hash:
                raise ValueError(f"Candidate completion report output changed: {rel}")
        aggregate_summary["completion_comparison"] = completion_provenance["summary"]
    runtime_audits = sorted(
        (ROOT / "brain6/qc").glob("power_screen_runtime_*.json"),
        key=lambda path: path.stat().st_mtime,
    )
    runtime_audit_path = runtime_audits[-1] if runtime_audits else None
    runtime_audit = read_json(runtime_audit_path) if runtime_audit_path else {}
    runtime_age_seconds = (
        max(0, int((datetime.now(timezone.utc) - datetime.fromisoformat(
            runtime_audit["generated_utc"]).astimezone(timezone.utc)).total_seconds()))
        if runtime_audit else None
    )
    runtime_is_fresh = runtime_age_seconds is not None and runtime_age_seconds <= 900
    runtime_status = (
        "RUNNING" if runtime_is_fresh and runtime_audit.get("status") == "RUNNING"
        and runtime_audit.get("active_workers") == 4 else "PARTIAL"
    )

    paths = {
        "global_map": ROOT / "brain6/results/global/brain6_72_locked.tsv",
        "global_source": ROOT / "results/atlas/trait_pairs.tsv",
        "novelty_crosswalk": ROOT / "brain6/results/novelty/brain6_significant_pair_novelty_crosswalk.tsv",
        "novelty_crosswalk_provenance": ROOT / "brain6/results/novelty/brain6_significant_pair_novelty_crosswalk.provenance.json",
        "novelty_crosswalk_builder": ROOT / "brain6/scripts/build_brain6_novelty_crosswalk.py",
        "novelty_crosswalk_tests": ROOT / "brain6/scripts/tests/test_build_brain6_novelty_crosswalk.py",
        "published_rg_context": ROOT / "brain6/results/novelty/jia_2025_sleep_psychiatric_rg_context.tsv",
        "published_rg_context_provenance": ROOT / "brain6/results/novelty/jia_2025_sleep_psychiatric_rg_context.provenance.json",
        "published_rg_context_builder": ROOT / "brain6/scripts/build_brain6_published_rg_context.py",
        "published_rg_context_tests": ROOT / "brain6/scripts/tests/test_build_brain6_published_rg_context.py",
        "targeted_literature_update": ROOT / "brain6/results/novelty/literature_update_20260926.md",
        "targeted_literature_update_table": ROOT / "brain6/results/novelty/literature_update_20260926.tsv",
        "targeted_literature_update_provenance": ROOT / "brain6/results/novelty/literature_update_20260926.provenance.json",
        "grover_sleep_pd_ldsc_context": ROOT / "brain6/results/novelty/grover_2022_sleep_pd_ldsc_context.tsv",
        "grover_sleep_pd_ldsc_context_provenance": ROOT / "brain6/results/novelty/grover_2022_sleep_pd_ldsc_context.provenance.json",
        "grover_sleep_pd_ldsc_context_builder": ROOT / "brain6/scripts/build_grover_2022_sleep_pd_ldsc_context.py",
        "grover_sleep_pd_article_xml": ROOT / "brain6/qc/replication_sources/grover_2022_article_fullTextXML.xml",
        "sleepchart_longsleep_rg_context": ROOT / "brain6/results/novelty/sleepchart_2026_long_sleep_rg_context.tsv",
        "sleepchart_longsleep_rg_context_provenance": ROOT / "brain6/results/novelty/sleepchart_2026_long_sleep_rg_context.provenance.json",
        "sleepchart_longsleep_rg_context_builder": ROOT / "brain6/scripts/build_sleepchart_2026_long_sleep_rg_context.py",
        "sleepchart_longsleep_rg_context_tests": ROOT / "brain6/scripts/tests/test_build_sleepchart_2026_long_sleep_rg_context.py",
        "partial_annotation_tier_policy": ROOT / "brain6/config/shared_locus_evidence_tiers_v1.json",
        "partial_annotation_gate": ROOT / "brain6/results/annotation/partial_candidate_region_gate_v1.tsv",
        "partial_annotation_gate_provenance": ROOT / "brain6/results/annotation/partial_candidate_region_gate_v1.provenance.json",
        "partial_annotation_gate_builder": ROOT / "brain6/scripts/build_partial_annotation_gate.py",
        "partial_nearest_gene_context": ROOT / "brain6/results/annotation/partial_candidate_lead_nearest_gene_grch37_v1.tsv",
        "partial_nearest_gene_provenance": ROOT / "brain6/results/annotation/partial_candidate_lead_nearest_gene_grch37_v1.provenance.json",
        "partial_nearest_gene_raw_responses": ROOT / "brain6/results/annotation/partial_candidate_lead_ensembl_grch37_responses_v1.json",
        "partial_nearest_gene_builder": ROOT / "brain6/scripts/annotate_partial_candidate_leads_ensembl_grch37.py",
        "partial_region_overlap": ROOT / "brain6/results/annotation/partial_candidate_region_overlap_v1.tsv",
        "partial_region_overlap_provenance": ROOT / "brain6/results/annotation/partial_candidate_region_overlap_v1.provenance.json",
        "partial_region_overlap_builder": ROOT / "brain6/scripts/build_partial_candidate_region_overlap.py",
        "partial_locus_lead_stability": ROOT / "brain6/results/loci/placo_candidate_locus_lead_stability_v1.tsv",
        "partial_locus_lead_stability_provenance": ROOT / "brain6/results/loci/placo_candidate_locus_lead_stability_v1.provenance.json",
        "partial_locus_lead_stability_builder": ROOT / "brain6/scripts/compare_partial_placo_locus_lead_stability.py",
        "partial_candidate_pvalue_validation": ROOT / "brain6/results/loci/placo_candidate_pvalue_validation_v1.tsv",
        "partial_candidate_pvalue_validation_provenance": ROOT / "brain6/results/loci/placo_candidate_pvalue_validation_v1.provenance.json",
        "partial_candidate_pvalue_validation_builder": ROOT / "brain6/scripts/validate_partial_placo_pvalues.py",
        "partial_candidate_pvalue_validation_r_helper": ROOT / "brain6/scripts/recompute_partial_placo_candidate_pvalues.R",
        "partial_candidate_pvalue_validation_summary": ROOT / "brain6/results/supplement/table_S39_partial_placo_pvalue_replay.tsv",
        "partial_placo_factor_sensitivity_provenance": ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/provenance.json",
        "parent_novelty_audit": ROOT / "results/analysis/literature_novelty_audit.tsv",
        "parent_novelty_queries": ROOT / "results/analysis/literature_search_queries.tsv",
        "supplement_index": ROOT / "brain6/paper/supplement.md",
        "hostile_review": ROOT / "brain6/paper/REVIEWER_2_AUDIT.md",
        "replication": ROOT / "brain6/results/replication/replication_master.tsv",
        "pgc_mdd2025_rg_sensitivity": ROOT / "brain6/results/replication/pgc_mdd2025_insomnia_rg_sensitivity_v1.tsv",
        "pgc_mdd2025_rg_sensitivity_provenance": ROOT / "brain6/results/replication/pgc_mdd2025_insomnia_rg_sensitivity_v1.provenance.json",
        "pgc_mdd2025_rg_sensitivity_runner": ROOT / "brain6/scripts/run_pgc_mdd2025_insomnia_rg_sensitivity.py",
        "pgc_mdd2025_rg_sensitivity_materializer": ROOT / "brain6/scripts/materialize_pgc_mdd2025_rg_record.py",
        "pgc_mdd2025_rg_sensitivity_validator": ROOT / "brain6/scripts/validate_pgc_mdd2025_rg_sensitivity.py",
        "replication_source_feasibility": ROOT / "brain6/qc/replication_source_feasibility_20260923.md",
        "mdd2025_no_ukbb_cohort_sidecar": ROOT / "brain6/qc/replication_sources/pgc-mdd2025_no23andMe-noUKBB_eur_v3.49.24.11.txt",
        "mdd2025_no_ukbb_cohort_overlap_audit": ROOT / "brain6/qc/replication_sources/pgc_mdd2025_noUKBB_overlap_audit.json",
        "scz_cohort_overlap_audit": ROOT / "brain6/qc/cohort_overlap_scz_followup_20260925.md",
        "scz_cohort_overlap_provenance": ROOT / "brain6/qc/cohort_overlap_scz_followup_20260925.provenance.json",
        "adhd_replication_overlap_audit": ROOT / "brain6/qc/cohort_overlap_insomnia_adhd_replication_20260925.md",
        "adhd_replication_overlap_provenance": ROOT / "brain6/qc/cohort_overlap_insomnia_adhd_replication_20260925.provenance.json",
        "finngen_r9_manifest": ROOT / "brain6/qc/replication_sources/finngen_r9_manifest.tsv",
        "finngen_r9_candidate_table": ROOT / "brain6/qc/finngen_r9_replication_candidates.tsv",
        "finngen_r9_candidate_provenance": ROOT / "brain6/qc/finngen_r9_replication_candidates.provenance.json",
        "finngen_r9_feasibility_report": ROOT / "brain6/qc/finngen_r9_replication_feasibility_20260925.md",
        "finngen_r9_validator": ROOT / "brain6/scripts/validate_finngen_r9_replication_candidates.py",
        "finngen_r9_validator_tests": ROOT / "brain6/scripts/tests/test_validate_finngen_r9_replication_candidates.py",
        "finngen_pd_r4_fallback_candidate": ROOT / "brain6/qc/finngen_pd_r4_fallback_candidate_20260925.tsv",
        "finngen_pd_r4_fallback_provenance": ROOT / "brain6/qc/finngen_pd_r4_fallback_candidate_20260925.provenance.json",
        "validation_makefile": ROOT / "brain6/Makefile",
        "track_b_archive_audit": ROOT / "brain6/qc/protected_track_b_archive_audit.md",
        "canonical_decision": V3_RUN / "canonical_family_decision.json",
        "canonical_aggregate": V3_RUN / "results/canonical_family_results.tsv",
        "canonical_causes": ROOT / "brain6/results/lava/canonical_v3_not_run_cells.tsv",
        "canonical_cause_summary": ROOT / "brain6/results/lava/canonical_v3_not_run_by_trait_v1.tsv",
        "roundoff_manifest": ROOT / "brain6/manifests/lava_roundoff_v1_run.json",
        "roundoff_audit": ROUND_AUDIT,
        "roundoff_comparison": ROOT / "brain6/results/lava/lava_baseline_roundoff_pair_locus_comparison_v1.tsv",
        "pair_context_audit": ROOT / "brain6/results/lava/lava_pair_context_univariates_v1.tsv",
        "pair_gate_probe": ROOT / "brain6/results/lava/lava_unique_gate_probe_v1.tsv",
        "pair_gate_probe_provenance": ROOT / "brain6/results/lava/lava_unique_gate_probe_v1.provenance.json",
        "placo_master": ROOT / "brain6/results/placo/placo_master.tsv",
        "placo_qc": ROOT / "brain6/results/placo/placo_v3_pair_qc_validation.json",
        "placo_candidates": ROOT / "brain6/results/loci/placo_candidate_loci_partial.tsv",
        "placo_valid_ld": ROOT / "brain6/results/loci/placo_candidate_valid_ld_partial_v2.provenance.json",
        "placo_ld_factor_diagnostic": ROOT / "brain6/qc/placo_ld_factor_normalization_v3.provenance.json",
        "placo_ld_block_psd_diagnostic": ROOT / "brain6/qc/placo_ld_primary_block_psd_v2.provenance.json",
        "placo_factor_sensitivity_provenance": ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/provenance.json",
        "placo_factor_sensitivity_variants": ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_variants.tsv",
        "placo_factor_sensitivity_loci": ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/placo_candidate_loci.tsv",
        "placo_factor_sensitivity_edges": ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/normalized_range_exceptions.tsv",
        "placo_factor_sensitivity_summary": ROOT / "brain6/results/loci/placo_factor_normalized_sensitivity_v2/pair_summary.tsv",
        "placo_factor_sensitivity_builder": ROOT / "brain6/scripts/build_placo_factor_normalized_sensitivity_v1.py",
        "placo_factor_sensitivity_validator": ROOT / "brain6/scripts/validate_placo_factor_normalized_sensitivity_v1.py",
        "placo_factor_sensitivity_tests": ROOT / "brain6/scripts/tests/test_validate_placo_factor_normalized_sensitivity_v1.py",
        "cross_layer": ROOT / "brain6/results/evidence/brain6_cross_layer_evidence.tsv",
        "readiness_report": ROOT / "brain6/FINAL_READINESS_REPORT.md",
        "scope_audit": ROOT / "brain6/SCOPE_AUDIT.md",
        "probe_script": ROOT / "brain6/scripts/probe_lava_unique_gate_v1.R",
        "collector_script": Path(__file__),
        "all_locked_sleep_power_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/all_locked_sleep_trait_power_audit_v1.tsv",
        "all_locked_sleep_power_audit_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/all_locked_sleep_trait_power_audit_v1.provenance.json",
        "all_locked_sleep_power_audit_builder": ROOT / "brain6/scripts/build_locked_sleep_panel_power_audit.py",
        "all_locked_sleep_power_audit_tests": ROOT / "brain6/scripts/tests/test_build_locked_sleep_panel_power_audit.py",
        "canonical_family_trait_power_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/canonical_family_trait_power_audit_v1.tsv",
        "canonical_family_trait_power_audit_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/canonical_family_trait_power_audit_v1.provenance.json",
        "canonical_family_trait_power_audit_builder": ROOT / "brain6/scripts/build_canonical_family_trait_power_audit.py",
        "power_screen": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_sleep_trait_audit.tsv",
        "power_screen_complete_trait_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_sleep_trait_audit_complete_v2.tsv",
        "power_screen_historical_trait_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_sleep_trait_audit_complete.tsv",
        "power_screen_canonical_metrics_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_sleep_trait_audit_complete_v2.provenance.json",
        "power_screen_canonical_metrics_builder": ROOT / "brain6/scripts/build_current_sleep_trait_audit_canonical_metrics.py",
        "power_screen_reference_overlap": ROOT / "brain6/results/power_optimized_sensitivity_v1/canonical_sleep_reference_overlap_v1.tsv",
        "power_screen_reference_overlap_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/canonical_sleep_reference_overlap_v1.provenance.json",
        "power_screen_reference_overlap_builder": ROOT / "brain6/scripts/audit_canonical_sleep_reference_overlap.py",
        "power_candidates": ROOT / "brain6/results/power_optimized_sensitivity_v1/candidate_replacement_screen_complete_v2.tsv",
        "power_candidates_historical": ROOT / "brain6/results/power_optimized_sensitivity_v1/candidate_replacement_screen.tsv",
        "power_candidates_v2_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/candidate_replacement_screen_complete_v2.provenance.json",
        "power_candidates_v2_builder": ROOT / "brain6/scripts/build_candidate_replacement_screen_complete_v2.py",
        "power_candidates_v3": ROOT / "brain6/results/power_optimized_sensitivity_v1/candidate_replacement_screen_complete_v3.tsv",
        "power_candidates_v3_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/candidate_replacement_screen_complete_v3.provenance.json",
        "power_candidates_v3_builder": ROOT / "brain6/scripts/reassess_candidate_replacement_decisions_v1.py",
        "power_candidates_v3_tests": ROOT / "brain6/scripts/tests/test_reassess_candidate_replacement_decisions_v1.py",
        "power_candidate_decision_addendum": ROOT / "brain6/results/power_optimized_sensitivity_v1/continuation_decision_addendum_2026-09-26.md",
        "power_screen_report": ROOT / "brain6/results/power_optimized_sensitivity_v1/README.md",
        "power_screen_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/provenance.json",
        "power_screen_builder": ROOT / "brain6/scripts/build_power_optimized_sleep_audit.py",
        "power_screen_validator": ROOT / "brain6/scripts/validate_power_optimized_sleep_audit.py",
        "power_screen_literature_addendum": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_literature_addendum_2026-09-26.md",
        "power_screen_literature_table": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_literature_addendum_2026-09-26.tsv",
        "power_screen_literature_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/current_literature_addendum_2026-09-26.provenance.json",
        **({"power_screen_runtime_audit": runtime_audit_path} if runtime_audit_path else {}),
        "power_screen_runtime_audit_builder": ROOT / "brain6/scripts/audit_power_screen_runtime.py",
        "power_candidate_source_config": ROOT / "brain6/config/power_optimized_sensitivity_v1/source_sleep_duration_continuous_dashti2019.json",
        "power_candidate_normalization_receipt": ROOT / "brain6/results/power_optimized_sensitivity_v1/normalized/sleep_duration_continuous_dashti_2019/receipt.json",
        "power_candidate_materialization": ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_inputs_v1/materialization.provenance.json",
        "power_candidate_reference_coverage": ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_inputs_v1/reference_coverage.tsv",
        "power_candidate_lava_worker": ROOT / "brain6/scripts/run_power_optimized_sleep_univariate_v1.R",
        "power_candidate_lava_config": ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/config.json",
        **{f"power_candidate_worker_{i}_config": screen_root / f"worker_{i}.config.json"
           for i in range(1, 5)},
        "power_candidate_lava_worker_1_config": ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/worker_1.config.json",
        "power_candidate_lava_worker_2_config": ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/worker_2.config.json",
        "power_candidate_lava_worker_3_config": ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/worker_3.config.json",
        "power_candidate_lava_worker_4_config": ROOT / "brain6/results/power_optimized_sensitivity_v1/lava_trait_screen_v1/worker_4.config.json",
        "power_promotion_rule_reassessment": ROOT / "brain6/results/power_optimized_sensitivity_v1/promotion_rule_reassessment_2026-09-26.json",
        "power_overlap_integrity_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/overlap_integrity_audit_v1.json",
        "power_overlap_readiness_report": ROOT / "brain6/results/power_optimized_sensitivity_v1/overlap_readiness_v1.md",
        "power_overlap_external_provenance": Path("/Volumes/Extreme SSD/brain6-work/lava-power-optimized-sensitivity-v1/overlap_v1_retry/overlap.provenance.json"),
        "power_overlap_validator": ROOT / "brain6/scripts/validate_power_optimized_lava_overlap_v1.py",
        "power_overlap_validator_tests": ROOT / "brain6/scripts/tests/test_validate_power_optimized_lava_overlap_v1.py",
        "power_family_feasibility_audit": ROOT / "brain6/results/power_optimized_sensitivity_v1/sensitivity_family_feasibility_v1.json",
        "power_family_feasibility_builder": ROOT / "brain6/scripts/audit_power_optimized_family_feasibility_v1.py",
        "power_family_feasibility_tests": ROOT / "brain6/scripts/tests/test_audit_power_optimized_family_feasibility_v1.py",
    }
    if aggregate_summary is not None:
        paths.update({
            "power_candidate_aggregate_tsv": aggregate_tsv,
            "power_candidate_aggregate_summary": aggregate_summary_path,
            "power_candidate_aggregate_provenance": aggregate_provenance_path,
            "power_candidate_aggregate_builder": ROOT / "brain6/scripts/aggregate_power_optimized_sleep_screen_v1.py",
            "power_screen_completion_table": ROOT / "brain6/results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.tsv",
            "power_screen_completion_report": ROOT / "brain6/results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.md",
            "power_screen_completion_provenance": ROOT / "brain6/results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.provenance.json",
            "power_screen_completion_builder": ROOT / "brain6/scripts/build_power_screen_completion_report_v1.py",
        })
        for index, (result_path, summary_path) in enumerate(
                zip(worker_result_paths, worker_summary_paths, strict=True), start=1):
            if result_path.is_file() and summary_path.is_file():
                paths[f"power_candidate_worker_{index}_result"] = result_path
                paths[f"power_candidate_worker_{index}_summary"] = summary_path
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required status evidence is missing: " + ", ".join(missing))

    decision = read_json(paths["canonical_decision"])
    global_rows = read_tsv(paths["global_map"])
    reference_overlap_provenance = read_json(paths["power_screen_reference_overlap_provenance"])
    promotion_reassessment = read_json(paths["power_promotion_rule_reassessment"])
    overlap_audit = read_json(paths["power_overlap_integrity_audit"])
    family_feasibility = read_json(paths["power_family_feasibility_audit"])
    candidate_decision = read_json(paths["power_candidates_v3_provenance"])
    for relative, expected_hash in candidate_decision.get("input_sha256", {}).items():
        candidate_input = ROOT / relative
        if not candidate_input.is_file() or sha256(candidate_input) != expected_hash:
            raise ValueError(f"Candidate decision reassessment input changed: {relative}")
    for relative, expected_hash in candidate_decision.get("output_sha256", {}).items():
        candidate_output = ROOT / relative
        if not candidate_output.is_file() or sha256(candidate_output) != expected_hash:
            raise ValueError(f"Candidate decision reassessment output changed: {relative}")
    candidate_rows_v3 = read_tsv(paths["power_candidates_v3"])
    dashti_decision_v3 = next(row["decision"] for row in candidate_rows_v3
                              if row["candidate"].startswith("Dashti 2019"))
    overlap_provenance = read_json(paths["power_overlap_external_provenance"])
    if (promotion_reassessment.get("decision_status") != "FULL_SENSITIVITY_TECHNICAL_ADVANCEMENT_JUSTIFIED" or
            promotion_reassessment.get("association_results_consulted") is not False or
            overlap_audit.get("status") != "PASS_OVERLAP_AUDIT_FAMILY_FAILS_FROZEN_QC" or
            overlap_audit.get("family_decision", {}).get("bivariate_lava_started") is not False or
            family_feasibility.get("status") != "PASS_SLEEP_ONLY_REPLACEMENTS_CANNOT_PASS_FROZEN_FAMILY_QC" or
            family_feasibility.get("association_results_consulted") is not False or
            candidate_decision.get("status") != "PASS_TECHNICAL_ADVANCEMENT_FAMILY_QC_BLOCKED" or
            candidate_decision.get("technical_advancement_candidates") != 1 or
            candidate_decision.get("full_family_qc_passing_candidates") != 0 or
            candidate_decision.get("bivariate_lava_started") is not False or
            len(candidate_rows_v3) != 6 or
            "TECHNICAL ADVANCEMENT CRITERIA MET" not in dashti_decision_v3 or
            "Bivariate LAVA remains blocked" not in dashti_decision_v3 or
            family_feasibility.get("best_case_lower_bounds", {}).get("one_longsleep_replacement_zero_not_run", {}).get("excess_over_frozen_limit") != 1556 or
            family_feasibility.get("best_case_lower_bounds", {}).get("both_sleep_traits_zero_not_run", {}).get("excess_over_frozen_limit") != 885 or
            overlap_provenance.get("status") != "PASS_LDSC_INTERCEPTS_AND_MATRICES_CREATED"):
        raise ValueError("Continuous-duration sensitivity decision or LDSC overlap audit is not in its expected state")
    complete_sleep_audit = read_tsv(paths["power_screen_complete_trait_audit"])
    if (reference_overlap_provenance.get("status") != "PASS_HASH_VERIFIED_EXACT_ID_OVERLAP" or
            {row["trait"] for row in complete_sleep_audit} != {"insomnia", "longsleep"}):
        raise ValueError("Canonical sleep source/reference overlap audit is incomplete")
    novelty_provenance = read_json(paths["novelty_crosswalk_provenance"])
    grover_provenance = read_json(paths["grover_sleep_pd_ldsc_context_provenance"])
    grover_rows = read_tsv(paths["grover_sleep_pd_ldsc_context"])
    grover_source_context = grover_provenance.get("source", {})
    for relative, expected_hash in grover_provenance.get("input_hashes", {}).items():
        source_path = ROOT / relative
        if not source_path.is_file() or sha256(source_path) != expected_hash:
            raise ValueError(f"Grover context input changed: {relative}")
    for relative, expected_hash in grover_provenance.get("source_file_hashes", {}).items():
        source_path = ROOT / relative
        if not source_path.is_file() or sha256(source_path) != expected_hash:
            raise ValueError(f"Grover article source XML changed: {relative}")
    grover_supplementary = grover_source_context.get("linked_supplementary_file", {})
    if (grover_provenance.get("status") != "PASS_PRIOR_NEGATIVE_LDSC_CONTEXT" or
            grover_provenance.get("builder_sha256") != sha256(paths["grover_sleep_pd_ldsc_context_builder"]) or
            grover_provenance.get("output", {}).get("sha256") != sha256(paths["grover_sleep_pd_ldsc_context"]) or
            len(grover_rows) != 1 or
            grover_rows[0].get("sleep_trait") != "longsleep" or
            grover_rows[0].get("brain_disorder") != "parkinson" or
            grover_rows[0].get("published_rg") != "NA" or
            grover_rows[0].get("numeric_pair_estimate_reviewed") != "False" or
            grover_source_context.get("numeric_rg_se_p_in_article_text") is not False or
            grover_source_context.get("linked_supplementary_workbook_reviewed") is not False or
            grover_supplementary.get("filename") != "Table_1.xlsx" or
            grover_supplementary.get("size_bytes") != 722516 or
            grover_supplementary.get("binary_retrieved") is not False or
            grover_supplementary.get("numeric_pair_estimate_reviewed") is not False):
        raise ValueError("Prior negative long-sleep–Parkinson LDSC context is missing or misrepresented")
    targeted_literature_rows = read_tsv(paths["targeted_literature_update_table"])
    if not any(row.get("source_id") == "grover_2022_sleep_pd_ldsc" for row in targeted_literature_rows):
        raise ValueError("Targeted literature update omits Grover et al. negative LDSC context")
    sleepchart_provenance = read_json(paths["sleepchart_longsleep_rg_context_provenance"])
    sleepchart_rows = read_tsv(paths["sleepchart_longsleep_rg_context"])
    for relative, expected_hash in sleepchart_provenance.get("input_hashes", {}).items():
        source_path = ROOT / relative
        if not source_path.is_file() or sha256(source_path) != expected_hash:
            raise ValueError(f"SleepChart long-sleep context input changed: {relative}")
    if (sleepchart_provenance.get("status") != "PASS_PRIOR_CONTEXT_ROUNDED_ESTIMATES_NOT_INDEPENDENT_REPLICATION" or
            sleepchart_provenance.get("builder_sha256") != sha256(paths["sleepchart_longsleep_rg_context_builder"]) or
            sleepchart_provenance.get("output", {}).get("sha256") != sha256(paths["sleepchart_longsleep_rg_context"]) or
            len(sleepchart_rows) != 4 or
            {row.get("brain_disorder") for row in sleepchart_rows} != {"adhd", "mdd", "scz", "bipolar"} or
            any(row.get("direction_concordant") != "True" or
                row.get("independent_replication") != "False" or
                row.get("power_replacement_eligible") != "False" for row in sleepchart_rows) or
            sleepchart_provenance.get("scope", {}).get("novelty_or_fdr_changed") is not False or
            sleepchart_provenance.get("source", {}).get("sleep_gwas", {}).get("cases") != 25049 or
            sleepchart_provenance.get("source", {}).get("sleep_gwas", {}).get("controls") != 300420):
        raise ValueError("SleepChart rounded prior-context estimates are missing or overinterpreted")
    if not any(row.get("source_id") == "sleepchart_2026_long_sleep_rg" for row in targeted_literature_rows):
        raise ValueError("Targeted literature update omits SleepChart 2026 long-sleep context")
    targeted_literature_provenance = read_json(paths["targeted_literature_update_provenance"])
    for relative, expected_hash in targeted_literature_provenance.get("outputs", {}).items():
        output_path = ROOT / relative
        if not output_path.is_file() or sha256(output_path) != expected_hash:
            raise ValueError(f"Targeted literature output changed: {relative}")
    annotation_gate_provenance = read_json(paths["partial_annotation_gate_provenance"])
    annotation_gate_rows = read_tsv(paths["partial_annotation_gate"])
    annotation_provenance = read_json(paths["partial_nearest_gene_provenance"])
    annotation_rows = read_tsv(paths["partial_nearest_gene_context"])
    annotation_raw = read_json(paths["partial_nearest_gene_raw_responses"])
    if (annotation_gate_provenance.get("status") != "PASS_NOT_TIERED_PREREQUISITES_BLOCKED" or
            annotation_gate_provenance.get("output", {}).get("sha256") != sha256(paths["partial_annotation_gate"]) or
            annotation_gate_provenance.get("builder", {}).get("sha256") != sha256(paths["partial_annotation_gate_builder"]) or
            len(annotation_gate_rows) != 19 or
            any(row.get("evidence_tier") != "NOT_TIERED" or row.get("final_shared_locus") != "False"
                for row in annotation_gate_rows) or
            annotation_provenance.get("status") != "PASS_DESCRIPTIVE_POSITIONAL_ANNOTATION_NOT_GENE_PRIORITIZATION" or
            annotation_provenance.get("output", {}).get("sha256") != sha256(paths["partial_nearest_gene_context"]) or
            annotation_provenance.get("builder", {}).get("sha256") != sha256(paths["partial_nearest_gene_builder"]) or
            annotation_provenance.get("raw_response_snapshot", {}).get("sha256") != sha256(paths["partial_nearest_gene_raw_responses"]) or
            len(annotation_rows) != 21 or len(annotation_raw.get("records", [])) != 21 or
            annotation_provenance.get("scope", {}).get("evidence_tier_changed") is not False or
            annotation_provenance.get("scope", {}).get("causal_gene_claims") is not False or
            annotation_provenance.get("scope", {}).get("fine_mapping_performed") is not False or
            annotation_provenance.get("scope", {}).get("colocalization_performed") is not False or
            any(row.get("evidence_tier") != "NOT_TIERED" for row in annotation_rows)):
        raise ValueError("Partial candidate annotation gate/context is missing, stale, or overinterpreted")
    annotation_inputs = annotation_provenance.get("sources", {})
    for item in ("gate", "gate_provenance", "candidate_variants", "candidate_sensitivity_provenance"):
        record = annotation_inputs.get(item, {})
        source_path = ROOT / record.get("path", "")
        if not source_path.is_file() or sha256(source_path) != record.get("sha256"):
            raise ValueError(f"Partial candidate gene-context source changed: {item}")
    region_overlap_provenance = read_json(paths["partial_region_overlap_provenance"])
    region_overlap_rows = read_tsv(paths["partial_region_overlap"])
    if (region_overlap_provenance.get("status") != "PASS_DESCRIPTIVE_NOT_SHARED_LOCUS_EVIDENCE" or
            region_overlap_provenance.get("builder", {}).get("sha256") != sha256(paths["partial_region_overlap_builder"]) or
            region_overlap_provenance.get("output", {}).get("sha256") != sha256(paths["partial_region_overlap"]) or
            len(region_overlap_rows) != 3 or
            region_overlap_provenance.get("scope", {}).get("cross_pair_region_comparisons_same_chromosome") != 9 or
            region_overlap_provenance.get("scope", {}).get("exact_shared_lead_variant_ids") != ["rs1051168"] or
            region_overlap_provenance.get("scope", {}).get("shared_locus_claims") is not False or
            any(row.get("evidence_tier") != "NOT_TIERED" for row in region_overlap_rows)):
        raise ValueError("Cross-pair candidate-region overlap diagnostic is missing or overinterpreted")
    lead_stability_provenance = read_json(paths["partial_locus_lead_stability_provenance"])
    lead_stability_rows = read_tsv(paths["partial_locus_lead_stability"])
    if (lead_stability_provenance.get("status") != "PASS_DIAGNOSTIC_ONLY_NO_PROMOTION" or
            lead_stability_provenance.get("builder", {}).get("sha256") != sha256(paths["partial_locus_lead_stability_builder"]) or
            lead_stability_provenance.get("output", {}).get("sha256") != sha256(paths["partial_locus_lead_stability"]) or
            len(lead_stability_rows) != 25 or
            lead_stability_provenance.get("scope", {}).get("stable_all_three") != 21 or
            lead_stability_provenance.get("scope", {}).get("range_filter_only") != 4 or
            lead_stability_provenance.get("scope", {}).get("full_five_track_family") is not False or
            lead_stability_provenance.get("scope", {}).get("genotype_level_ld_validation") is not False or
            lead_stability_provenance.get("scope", {}).get("tier_promotion") is not False or
            {row["lead_variant"] for row in lead_stability_rows if row["stability_class"] == "RANGE_FILTER_ONLY"} !=
            {"rs6430538", "rs199535", "rs199534", "rs4938023"}):
        raise ValueError("Partial PLACO lead-stability diagnostic is missing, stale, or overinterpreted")
    candidate_pvalue_provenance = read_json(paths["partial_candidate_pvalue_validation_provenance"])
    candidate_pvalue_rows = read_tsv(paths["partial_candidate_pvalue_validation"])
    candidate_pvalue_summary_rows = read_tsv(paths["partial_candidate_pvalue_validation_summary"])
    candidate_pvalue_summary = candidate_pvalue_provenance.get("pair_summary", {})
    if (candidate_pvalue_provenance.get("status") != "PASS_PINNED_SOURCE_CANDIDATE_PVALUES_RECOMPUTED" or
            candidate_pvalue_provenance.get("output", {}).get("sha256") != sha256(paths["partial_candidate_pvalue_validation"]) or
            candidate_pvalue_provenance.get("summary_output", {}).get("sha256") != sha256(paths["partial_candidate_pvalue_validation_summary"]) or
            candidate_pvalue_provenance.get("sources", {}).get("builder", {}).get("sha256") != sha256(paths["partial_candidate_pvalue_validation_builder"]) or
            candidate_pvalue_provenance.get("sources", {}).get("replay_helper", {}).get("sha256") != sha256(paths["partial_candidate_pvalue_validation_r_helper"]) or
            len(candidate_pvalue_rows) != 2297 or
            len(candidate_pvalue_summary_rows) != 4 or
            any(row.get("validation_status") != "PASS" for row in candidate_pvalue_rows) or
            any(row.get("matching_p_values") != row.get("candidate_variants") or
                float(row.get("maximum_relative_error", "1")) >= 1e-8
                for row in candidate_pvalue_summary_rows) or
            sum(item.get("variants", 0) for item in candidate_pvalue_summary.values()) != 2297 or
            any(item.get("passed") != item.get("variants") or item.get("max_relative_error", 1) >= 1e-8
                for item in candidate_pvalue_summary.values()) or
            candidate_pvalue_provenance.get("scope", {}).get("full_five_track_family") is not False or
            candidate_pvalue_provenance.get("scope", {}).get("candidate_locus_promotion") is not False):
        raise ValueError("Partial PLACO candidate p-value replay is missing, stale, or overinterpreted")
    replication_rows = read_tsv(paths["replication"])
    pgc_sensitivity_rows = read_tsv(paths["pgc_mdd2025_rg_sensitivity"])
    pgc_sensitivity_provenance = read_json(paths["pgc_mdd2025_rg_sensitivity_provenance"])
    if (len(pgc_sensitivity_rows) != 1 or
            pgc_sensitivity_provenance.get("status") != "PASS_AGGREGATE_ONLY_EXTERNAL_SENSITIVITY_RECORD" or
            pgc_sensitivity_provenance.get("output_sha256") != sha256(paths["pgc_mdd2025_rg_sensitivity"]) or
            pgc_sensitivity_rows[0].get("status") != "SENSITIVITY_ONLY_NOT_INDEPENDENT_REPLICATION" or
            pgc_sensitivity_rows[0].get("rg") != "0.4771" or
            pgc_sensitivity_rows[0].get("se") != "0.0228" or
            pgc_sensitivity_rows[0].get("p") != "4.6535e-97"):
        raise ValueError("PGC MDD2025 insomnia sensitivity record is missing, stale, or mislabeled")
    placo_rows = read_tsv(paths["placo_master"])
    cross_layer = read_tsv(paths["cross_layer"])
    placo_qc = read_json(paths["placo_qc"])
    round_audit = read_json(ROUND_AUDIT)
    roundoff = round_audit["lava_roundoff"]
    probe = read_tsv(paths["pair_gate_probe"])

    stages = [
        {"stage": "locked_global_rg", "status": "COMPLETE_VALID",
         "summary": "72 locked Brain6 pairs; significance inherited from the original 396-test family; all 35 significant pairs map to the dated parent-atlas literature audit.",
         "metrics": {"rows": len(global_rows), "inherited_significant": sum(
             row["significance_under_original_396_family"].lower() == "true" for row in global_rows),
                     "novelty_audit_search_date": novelty_provenance["summary"]["parent_audit_search_as_of_date"],
                     "targeted_literature_context_sources": len(targeted_literature_rows),
                     "grover_longsleep_parkinson_prior_ldsc_status": grover_rows[0]["source_pair_status"],
                     "grover_longsleep_parkinson_numeric_rg_available": False,
                     "grover_literature_context_independent": False,
                     "sleepchart_longsleep_rg_context_rows": len(sleepchart_rows),
                     "sleepchart_longsleep_rg_directions_concordant": sum(
                         row["direction_concordant"] == "True" for row in sleepchart_rows),
                     "sleepchart_longsleep_rg_independent_replication": False,
                     "sleepchart_longsleep_replacement_eligible": False,
                     "novelty_direct_rg_comparators": novelty_provenance["summary"]["direct_rg_count"],
                     "novelty_apparently_novel_pairs": novelty_provenance["summary"]["apparent_novelty_count"],
                     "novelty_classification_counts": novelty_provenance["summary"]["classification_counts"],
                     "published_rg_context_pairs": 21,
                     "published_rg_context_direction_concordant": 21,
                     "published_rg_context_comparison_independent": False},
         "evidence": ["global_map", "global_source", "novelty_crosswalk",
                      "novelty_crosswalk_provenance", "novelty_crosswalk_builder",
                      "novelty_crosswalk_tests", "parent_novelty_audit",
                      "parent_novelty_queries", "published_rg_context",
                      "published_rg_context_provenance", "published_rg_context_builder",
                      "published_rg_context_tests", "targeted_literature_update",
                      "targeted_literature_update_table", "targeted_literature_update_provenance",
                      "grover_sleep_pd_ldsc_context", "grover_sleep_pd_ldsc_context_provenance",
                      "grover_sleep_pd_ldsc_context_builder", "grover_sleep_pd_article_xml",
                      "sleepchart_longsleep_rg_context", "sleepchart_longsleep_rg_context_provenance",
                      "sleepchart_longsleep_rg_context_builder", "sleepchart_longsleep_rg_context_tests",
                      "supplement_index", "hostile_review"]},
        {"stage": "canonical_lava_v3", "status": "COMPLETE_FAILED_QC",
         "summary": "Immutable completed v3 family; no local-rg promotion.",
         "metrics": {"decision": decision["overall_status"], **decision["status_counts"],
                     "not_run_limit": decision["maximum_allowed_untested_cells"]},
         "evidence": ["canonical_decision", "canonical_aggregate", "canonical_causes", "canonical_cause_summary"]},
        {"stage": "lava_roundoff_recovery", "status": "PARTIAL",
         "summary": "All roundoff receipts verify, but pair-context univariate aggregation failed; no family decision.",
         "metrics": {"verified_loci": roundoff["verified_loci"], "planned_loci": roundoff["planned_loci"],
                     "receipt_errors": roundoff["receipt_errors"],
                     "status_counts": roundoff["pair_locus_status_counts"],
                     "family_decision": "NOT_CREATED_AGGREGATION_FAILED"},
         "evidence": ["roundoff_manifest", "roundoff_audit", "roundoff_comparison", "pair_context_audit", "pair_gate_probe", "pair_gate_probe_provenance"]},
        {"stage": "protected_track_b", "status": "BLOCKED_EXTERNAL",
         "summary": "Historical failed-partial artifact is non-importable; protected-pair rerun is prohibited.",
         "blocker": "No terminal-valid, source-complete, authorized immutable result package for insomnia–ADHD PLACO+.",
         "evidence": ["track_b_archive_audit"]},
        {"stage": "replication", "status": "PARTIAL",
         "summary": "One archived directional insomnia–FinnGen ADHD result and one aggregate PGC MDD2025 insomnia global-rg sensitivity are recorded. The PGC source excludes UK Biobank but nine cohort entries are also in locked Howard 2019 MDD; participant overlap remains unresolved, so this is not independent replication. Public cohort descriptions distinguish the UK Biobank insomnia, iPSYCH/deCODE/PGC ADHD discovery, and FinnGen register-based ADHD outcome, but ten PGC cohort identities and person-level overlap remain unverified. The PGC3 SCZ supplement names no UK Biobank cohort, but long sleep–SCZ participant overlap remains unknown. Four FinnGen R9 endpoints have manifest/object metadata only. A source-specific FinnGen R4 PD stratum is documented as a smaller independent-outcome fallback, but its standalone object is unverified and both R4/R9 require the official access workflow.",
         "metrics": {"rows": len(replication_rows), "completed_exact_or_directional": sum(
             row.get("replication_class", "").upper() in {"DIRECTIONAL_REPLICATION", "DIRECTIONAL_SUPPORT", "DIRECT_REPLICATION", "EXACT_REPLICATION"}
             for row in replication_rows), "finngen_r9_metadata_candidates": 4,
                     "pgc_mdd2025_insomnia_source_sensitivity_rg": float(pgc_sensitivity_rows[0]["rg"]),
                     "pgc_mdd2025_insomnia_source_sensitivity_status": pgc_sensitivity_rows[0]["status"],
                     "finngen_r9_gwas_files_acquired": 0,
                     "finngen_r4_pd_lower_power_fallbacks": 1,
                     "finngen_r4_pd_gwas_files_acquired": 0},
         "evidence": ["replication", "pgc_mdd2025_rg_sensitivity", "pgc_mdd2025_rg_sensitivity_provenance",
                      "pgc_mdd2025_rg_sensitivity_runner", "pgc_mdd2025_rg_sensitivity_materializer",
                      "pgc_mdd2025_rg_sensitivity_validator", "replication_source_feasibility",
                      "mdd2025_no_ukbb_cohort_sidecar", "mdd2025_no_ukbb_cohort_overlap_audit", "finngen_r9_manifest",
                      "scz_cohort_overlap_audit", "scz_cohort_overlap_provenance",
                      "adhd_replication_overlap_audit", "adhd_replication_overlap_provenance",
                      "finngen_r9_candidate_table", "finngen_r9_candidate_provenance",
                      "finngen_r9_feasibility_report", "finngen_r9_validator",
                      "finngen_r9_validator_tests", "finngen_pd_r4_fallback_candidate",
                      "finngen_pd_r4_fallback_provenance", "validation_makefile", "track_b_archive_audit"]},
        {"stage": "placo_and_locus_processing", "status": "PARTIAL",
         "summary": "Four pairs pass pair-level QC; protected fifth track is absent and LD independence is unresolved. Factor-normalized internal sensitivity matches original candidate leads/intervals, but genotype-level LD remains unvalidated.",
         "metrics": {"published_pairs": len(placo_rows), "pair_qc_status": placo_qc.get("status"),
                     "partial_candidate_intervals": 19, "invalid_ld_edges_preserved": 12317,
                     "factor_reconstruction_edges": 12317, "factor_sensitivity_edge_replacements": 157403,
                     "factor_sensitivity_roundoff_boundary_values": 1922,
                     "factor_sensitivity_leads_match_original": 21,
                     "factor_sensitivity_intervals_match_original": 19,
                     "promoted_from_diagnostic": 0},
         "evidence": ["placo_master", "placo_qc", "placo_candidates", "placo_valid_ld",
                      "placo_ld_factor_diagnostic", "placo_ld_block_psd_diagnostic",
                      "placo_factor_sensitivity_provenance", "placo_factor_sensitivity_variants",
                      "placo_factor_sensitivity_loci", "placo_factor_sensitivity_edges",
                      "placo_factor_sensitivity_summary", "placo_factor_sensitivity_builder",
                      "placo_factor_sensitivity_validator", "placo_factor_sensitivity_tests"]},
        {"stage": "fine_mapping_trait_coloc_qtl_regulatory_celltype", "status": "NOT_APPLICABLE",
         "summary": "Upstream complete-family eligibility gates are not met; no downstream loci are promoted.",
         "evidence": ["canonical_decision", "placo_qc", "readiness_report"]},
        {"stage": "cross_layer_integration", "status": "PARTIAL",
         "summary": "Five selected pairs remain genetic-only; no downstream molecular layer is admitted.",
         "metrics": {"pairs": len(cross_layer)}, "evidence": ["cross_layer"]},
        {"stage": "pair_independent_gate_probe", "status": "PARTIAL",
         "summary": "Two targeted local-rg cells were computed reproducibly; neither passes the full-family correction.",
         "metrics": {"cells": len(probe), "q_values": [float(row["q_bh_full_12475"]) for row in probe],
                     "promotion": "NONE"},
         "evidence": ["pair_gate_probe", "pair_gate_probe_provenance", "probe_script"]},
        {"stage": "sleep_power_optimized_screen", "status": runtime_status,
         "summary": ("The outcome-blinded continuous-duration screen meets the predeclared technical improvement examples. Candidate-specific LDSC overlap matrices pass source-bound integrity checks. The substituted seven-trait family fails the frozen 5% untested-cell ceiling; feasibility bounds show even perfect replacements of both sleep traits leave 885 excess NOT_RUN cells. Bivariate LAVA and downstream association analyses remain unstarted. Canonical insomnia and long-sleep source variants have exact hash-bound SNP-ID overlap counts against the sealed reference. Worker liveness is recorded separately and is not inferred from durable outputs."
                     if aggregate_summary is not None else
                     "Two Brain6 sleep traits audited; long sleep has greatest low-local-h2 burden. The continuous-duration candidate has four live workers verified in the latest runtime snapshot; no durable result TSVs or summaries exist yet."
                     if runtime_status == "RUNNING" else
                     "Two Brain6 sleep traits audited; the continuous-duration candidate is configured as four disjoint workers, but current liveness is incomplete or unverified."),
         "metrics": {"sleep_traits_audited": 2, "candidates_screened": len(read_tsv(paths["power_candidates"])),
                     "separate_literature_addendum_rows": 1,
                     "candidate_reference_overlap_fraction": 0.9604258571912654,
                     "canonical_sleep_reference_overlap_counts": {
                         row["trait"]: int(row["frozen_reference_overlap_count"])
                         for row in complete_sleep_audit},
                     "canonical_sleep_reference_overlap_status": reference_overlap_provenance["status"],
                     "candidate_screen_workers_configured": len(worker_config_data),
                     "candidate_worker_locus_assignments": [len(item["locus_ids"]) for item in worker_config_data],
                     "durable_worker_result_files": len(durable_worker_results),
                     "durable_worker_summary_files": len(durable_worker_summaries),
                     "running_workers": runtime_audit.get("active_workers", 0) if runtime_is_fresh else 0,
                     "runtime_audit_age_seconds": runtime_age_seconds,
                     "worker_liveness": "VERIFIED_RUNNING" if runtime_status == "RUNNING" else "NOT_VERIFIED",
                     "screen_completion_status": aggregate_summary.get("status") if aggregate_summary else "INCOMPLETE_OR_UNVERIFIED",
                     "screen_rows": aggregate_summary.get("rows") if aggregate_summary else None,
                     "screen_tested": aggregate_summary.get("tested") if aggregate_summary else None,
                     "screen_not_run": aggregate_summary.get("not_run") if aggregate_summary else None,
                     "screen_failed": aggregate_summary.get("failed") if aggregate_summary else None,
                     "screen_reason_counts": aggregate_summary.get("reason_counts") if aggregate_summary else None,
                     "screen_strict_gate_pass": aggregate_summary.get("strict_gate_pass") if aggregate_summary else None,
                     "screen_canonical_strict_gate_pass": aggregate_summary.get("canonical_longsleep_strict_gate_pass") if aggregate_summary else None,
                     "screen_additional_strict_gate_pass": aggregate_summary.get("additional_strict_gate_pass") if aggregate_summary else None,
                     "screen_minimum_additional_for_promotion": aggregate_summary.get("minimum_additional_for_promotion") if aggregate_summary else None,
                     "screen_eligible_under_historical_strict_gate_rule": aggregate_summary.get("eligible_for_full_sensitivity_screen_gate") if aggregate_summary else None,
                     "technical_advancement_trigger_met": promotion_reassessment.get("decision_status") == "FULL_SENSITIVITY_TECHNICAL_ADVANCEMENT_JUSTIFIED",
                     "technical_improvement_criteria_met": promotion_reassessment.get("qualifying_criteria", {}),
                     "substitution_family_not_run": overlap_audit["family_decision"]["family_not_run"],
                     "substitution_family_cells": overlap_audit["family_decision"]["family_cells"],
                     "substitution_family_not_run_fraction": overlap_audit["family_decision"]["family_not_run_fraction"],
                     "substitution_family_maximum_not_run_cells": overlap_audit["family_decision"]["maximum_not_run_cells"],
                     "substitution_family_passes_frozen_5_percent_qc": overlap_audit["family_decision"]["passes_frozen_5_percent_qc"],
                     "candidate_specific_overlap_matrices": len(overlap_provenance.get("pairs", [])),
                     "candidate_overlap_integrity_status": overlap_audit.get("status"),
                     "candidate_bivariate_lava_started": overlap_audit["family_decision"]["bivariate_lava_started"],
                     "sleep_replacement_family_feasibility_status": family_feasibility["status"],
                     "candidate_decision_status": candidate_decision["status"],
                     "candidate_technical_advancement_count": candidate_decision["technical_advancement_candidates"],
                     "candidate_full_family_qc_passing_count": candidate_decision["full_family_qc_passing_candidates"],
                     "best_case_one_longsleep_replacement_not_run": family_feasibility["best_case_lower_bounds"]["one_longsleep_replacement_zero_not_run"]["not_run"],
                     "best_case_one_longsleep_replacement_excess": family_feasibility["best_case_lower_bounds"]["one_longsleep_replacement_zero_not_run"]["excess_over_frozen_limit"],
                     "best_case_both_sleep_replacements_not_run": family_feasibility["best_case_lower_bounds"]["both_sleep_traits_zero_not_run"]["retained_five_disorder_traits_not_run"],
                     "best_case_both_sleep_replacements_excess": family_feasibility["best_case_lower_bounds"]["both_sleep_traits_zero_not_run"]["excess_over_frozen_limit"],
                     "screen_completion_comparison": aggregate_summary.get("completion_comparison") if aggregate_summary else None,
                     "sensitivity_families_launched": 0, "promotion": "NONE"},
         "evidence": ["power_screen", "power_screen_complete_trait_audit",
                      "power_screen_historical_trait_audit", "power_screen_canonical_metrics_provenance",
                      "power_screen_canonical_metrics_builder",
                      "power_screen_reference_overlap", "power_screen_reference_overlap_provenance",
                      "power_screen_reference_overlap_builder", "power_candidates", "power_candidates_historical",
                      "power_candidates_v2_provenance", "power_candidates_v2_builder", "power_screen_report",
                      "power_screen_provenance", "power_screen_builder", "power_screen_validator",
                      "power_screen_literature_addendum", "power_screen_literature_table",
                      "power_screen_literature_provenance",
                      "power_candidate_source_config", "power_candidate_normalization_receipt",
                      "power_candidate_materialization", "power_candidate_reference_coverage",
                      "power_candidate_lava_worker", "power_candidate_lava_config",
                      "power_candidate_worker_1_config", "power_candidate_worker_2_config",
                      "power_candidate_worker_3_config", "power_candidate_worker_4_config",
                      "power_screen_runtime_audit_builder",
                      "all_locked_sleep_power_audit", "all_locked_sleep_power_audit_provenance",
                      "all_locked_sleep_power_audit_builder", "all_locked_sleep_power_audit_tests",
                      "canonical_family_trait_power_audit", "canonical_family_trait_power_audit_provenance",
                      "canonical_family_trait_power_audit_builder",
                      "power_promotion_rule_reassessment", "power_overlap_integrity_audit",
                      "power_overlap_readiness_report", "power_overlap_external_provenance",
                      "power_overlap_validator", "power_overlap_validator_tests"] +
                     ["power_family_feasibility_audit", "power_family_feasibility_builder",
                      "power_family_feasibility_tests"] +
                     ["power_candidates_v3", "power_candidates_v3_provenance",
                      "power_candidates_v3_builder", "power_candidates_v3_tests",
                      "power_candidate_decision_addendum"] +
                     (["power_candidate_aggregate_tsv", "power_candidate_aggregate_summary",
                       "power_candidate_aggregate_provenance", "power_candidate_aggregate_builder",
                       "power_screen_completion_table", "power_screen_completion_report",
                       "power_screen_completion_provenance", "power_screen_completion_builder",
                       "power_candidate_worker_1_result", "power_candidate_worker_1_summary",
                       "power_candidate_worker_2_result", "power_candidate_worker_2_summary",
                       "power_candidate_worker_3_result", "power_candidate_worker_3_summary",
                       "power_candidate_worker_4_result", "power_candidate_worker_4_summary"]
                      if aggregate_summary is not None else []) +
                     (["power_screen_runtime_audit"] if runtime_audit_path else [])},
        {"stage": "manuscript_reporting_and_validation", "status": "PARTIAL",
         "summary": "The report remains provisional; external and scientific blockers remain.",
         "metrics": {"report_status": "PROVISIONAL; validation recorded separately",
                     "baseline_manifest_matches_receipts": round_audit["lava"]["manifest_matches_receipts"]},
         "evidence": ["readiness_report"]},
        {"stage": "partial_candidate_gene_context", "status": "COMPLETE_DESCRIPTIVE_DIAGNOSTIC",
         "summary": "Pinned-source p-value replay, nearest-gene context, cross-pair interval overlap, and three-method PLACO lead stability are recorded after the hash-bound NOT_TIERED gate; LD is not genotype-level validated, and method-sensitive leads are not validated independent loci or gene prioritization.",
         "metrics": {"candidate_regions_gated": len(annotation_gate_rows),
                     "regions_not_tiered": sum(row["evidence_tier"] == "NOT_TIERED" for row in annotation_gate_rows),
                     "lead_signals_annotated": len(annotation_rows),
                     "cross_pair_same_chromosome_comparisons": region_overlap_provenance["scope"]["cross_pair_region_comparisons_same_chromosome"],
                     "cross_pair_overlapping_region_pairs": len(region_overlap_rows),
                     "exact_shared_lead_variant_ids": region_overlap_provenance["scope"]["exact_shared_lead_variant_ids"],
                     "lead_signals_stable_across_ld_methods": lead_stability_provenance["scope"]["stable_all_three"],
                     "range_filtered_method_specific_leads": lead_stability_provenance["scope"]["range_filter_only"],
                     "method_specific_lead_ids": sorted(row["lead_variant"] for row in lead_stability_rows
                                                         if row["stability_class"] == "RANGE_FILTER_ONLY"),
                     "candidate_placo_pvalues_recomputed": len(candidate_pvalue_rows),
                     "candidate_placo_pvalues_matching": sum(row["validation_status"] == "PASS" for row in candidate_pvalue_rows),
                     "candidate_placo_pvalue_max_relative_error": max(item["max_relative_error"] for item in candidate_pvalue_summary.values()),
                     "assembly": annotation_provenance["scope"]["assembly"],
                     "causal_gene_claims": False, "fine_mapping_performed": False,
                     "colocalization_performed": False, "qtl_or_tissue_prioritization_performed": False},
         "evidence": ["partial_annotation_tier_policy", "partial_annotation_gate",
                      "partial_annotation_gate_provenance", "partial_annotation_gate_builder",
                      "partial_nearest_gene_context", "partial_nearest_gene_provenance",
                      "partial_nearest_gene_raw_responses", "partial_nearest_gene_builder",
                      "partial_placo_factor_sensitivity_provenance", "partial_region_overlap",
                      "partial_region_overlap_provenance", "partial_region_overlap_builder",
                      "partial_locus_lead_stability", "partial_locus_lead_stability_provenance",
                      "partial_locus_lead_stability_builder", "partial_candidate_pvalue_validation",
                      "partial_candidate_pvalue_validation_provenance", "partial_candidate_pvalue_validation_builder",
                      "partial_candidate_pvalue_validation_summary",
                      "partial_candidate_pvalue_validation_r_helper"]},
    ]
    status_text = git("status", "--porcelain")
    inventory = {
        "schema_version": 1,
        "inventory_id": "brain6-only-project-state-20260925",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "repository": {"root": str(ROOT), "branch": git("branch", "--show-current"),
                       "head": git("rev-parse", "HEAD"), "working_tree_dirty": bool(status_text)},
        "status_values": ["COMPLETE_VALID", "COMPLETE_FAILED_QC", "COMPLETE_DESCRIPTIVE_DIAGNOSTIC", "RUNNING", "PARTIAL",
                          "READY_TO_RUN", "BLOCKED_EXTERNAL", "NOT_APPLICABLE"],
        "stages": stages,
        "evidence_files": {name: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
                           for name, path in paths.items()},
        "interpretation": "Current evidence inventory, not a declaration of project completion or promotion of diagnostic results.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "stages": len(stages),
                      "output_sha256": sha256(output)}, indent=2))


if __name__ == "__main__":
    main()
