#!/usr/bin/env python3
"""Validate Brain6 global extraction and follow-up lock without recomputing inference."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

from build_brain6_novelty_crosswalk import validate as validate_brain6_novelty_crosswalk
from build_brain6_published_rg_context import (
    run as validate_brain6_published_rg_context,
)
from build_lava_canonical_v3_cause_summaries import (
    summarize as summarize_lava_v3_causes,
)
from build_profile_leave_one_out import FIELDS as PROFILE_LOO_FIELDS
from build_profile_leave_one_out import render_table as render_profile_loo
from build_profile_spearman import render_table
from validate_placo_factor_normalized_sensitivity_v1 import (
    validate as validate_placo_factor_sensitivity,
)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def profile_loo_tables_equivalent(actual: list[dict[str, str]], expected: list[dict[str, str]],
                                 absolute_tolerance: float = 5e-12) -> bool:
    """Compare frozen descriptive values at their serialized 12-digit precision."""
    if len(actual) != len(expected):
        return False
    numeric_fields = {
        "full_profile_r", "leave_one_out_r", "delta_from_full",
        "leave_one_out_min_r", "leave_one_out_max_r", "leave_one_out_range_width",
    }
    for actual_row, expected_row in zip(actual, expected):
        if actual_row.keys() != expected_row.keys():
            return False
        for field in actual_row:
            left, right = actual_row[field], expected_row[field]
            if field not in numeric_fields:
                if left != right:
                    return False
                continue
            try:
                left_value, right_value = float(left), float(right)
            except ValueError:
                return False
            if (not math.isfinite(left_value) or not math.isfinite(right_value) or
                    abs(left_value - right_value) > absolute_tolerance):
                return False
    return True


def validate_pgc_mdd2025_cohort_overlap(root: Path = ROOT) -> dict[str, object]:
    """Verify the public no-UKB MDD2025 cohort sidecar and its PGC29 comparison."""
    sidecar = root / "brain6/qc/replication_sources/pgc-mdd2025_no23andMe-noUKBB_eur_v3.49.24.11.txt"
    audit_path = root / "brain6/qc/replication_sources/pgc_mdd2025_noUKBB_overlap_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    raw = sidecar.read_bytes()
    expected = {
        "cof3": "COFAMS", "grdg": "DGN", "nes1": "NESDA/NTR",
        "qi3c": "QIMR I317", "qi6c": "QIMR I610", "qio2": "QIMR COEX",
        "rad3": "RADIANT-UK", "shp0": "SHIP-LEGEND", "shpt": "SHIP-TREND-0",
    }
    rows = read_tsv(sidecar)
    row_by_id = {row["cohort"]: row for row in rows}
    matched = audit.get("matched_legacy_pgc29_cohorts", [])
    observed_matches = {row.get("cohort_id"): row for row in matched}
    checks = {
        "schema": audit.get("schema_version") == 1,
        "audit_id": audit.get("audit_id") == "pgc_mdd2025_no23andme_noUKBB_cohort_overlap_v1",
        "source_url": audit.get("source_url") == "https://ndownloader.figshare.com/files/51486998",
        "file_id": audit.get("figshare_file_id") == 51486998,
        "license": audit.get("license") == "CC0",
        "filename": audit.get("filename") == sidecar.name,
        "bytes": audit.get("bytes") == len(raw) == 7981,
        "published_md5": audit.get("published_md5") == "3367c2fd90cba302b1e59f06e1448927",
        "observed_md5": audit.get("observed_md5") == hashlib.md5(raw, usedforsecurity=False).hexdigest(),
        "sha256": audit.get("sha256") == sha256(sidecar),
        "row_count": audit.get("rows") == len(rows) == 74,
        "decision": audit.get("decision") == (
            "DOCUMENTED_COHORT_LEVEL_OVERLAP_WITH_LOCKED_HOWARD_2019_PGC_COMPONENT; "
            "NOT_INDEPENDENT_REPLICATION"),
        "match_ids": set(observed_matches) == set(expected),
        "match_count": audit.get("matched_cohort_count") == len(expected) == len(matched),
    }
    if not all(checks.values()):
        raise ValueError("PGC MDD2025 overlap audit metadata or source checksum is invalid: " +
                         ", ".join(name for name, ok in checks.items() if not ok))
    for cohort_id, expected_name in expected.items():
        recorded = observed_matches[cohort_id]
        source_row = row_by_id.get(cohort_id)
        if (source_row is None or source_row["study_name"] != expected_name or
                recorded.get("mdd2025_study_name") != expected_name or
                recorded.get("ancestry") != source_row["ancestry"] or
                recorded.get("phenotype") != source_row["phenotype"] or
                recorded.get("N_cases") != int(source_row["N_cases"]) or
                recorded.get("N_controls") != int(source_row["N_controls"])):
            raise ValueError(f"PGC MDD2025 cohort comparison differs from exact sidecar row {cohort_id}")
    if audit.get("limitations", []).count(
            "Cohort metadata demonstrates shared cohort membership, not participant-level identity or the exact overlap count.") != 1:
        raise ValueError("PGC MDD2025 overlap audit must preserve participant-level uncertainty")
    return {"status": "PASS_COHORT_OVERLAP_NOT_INDEPENDENT_REPLICATION",
            "sidecar_rows": len(rows), "matched_pgc29_cohorts": len(expected),
            "sidecar_sha256": sha256(sidecar), "audit_sha256": sha256(audit_path)}


def validate_evidence_tier_policy() -> str:
    """Enforce the frozen pre-annotation tier policy and its source bindings."""
    policy_path = OUT / "config/shared_locus_evidence_tiers_v1.json"
    checksum_path = OUT / "config/shared_locus_evidence_tiers_v1.json.sha256"
    policy_hash = sha256(policy_path)
    assert checksum_path.read_text(encoding="utf-8").split()[0] == policy_hash
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    assert policy["rule_id"] == "brain6_shared_locus_evidence_tiers_v1"
    assert policy["status"] == "FROZEN_PRE_ANNOTATION"
    assert policy["frozen_before_functional_annotation"] is True
    assert policy["functional_annotation_reviewed_before_freeze"] is False
    assert set(policy["tiers"]) == {"TIER_1", "TIER_2", "TIER_3", "NOT_TIERED"}
    assert policy["output_schema"] == [
        "pair", "disorder", "sleep_trait", "chr", "start", "end", "lead_variant",
        "PLACO_P", "local_rg_LAVA", "local_rg_HDLL", "replication_level",
        "evidence_tier", "notes",
    ]
    for record in policy["source_locks"].values():
        path = (ROOT / record["path"]).resolve()
        path.relative_to(ROOT.resolve())
        assert path.is_file() and sha256(path) == record["sha256"]
    # No final tier table may appear before complete PLACO/LAVA families exist.
    assert policy["current_assignment_status"] == "NONE_ASSIGNED_UPSTREAM_FAMILIES_INCOMPLETE"
    assert not (OUT / "results/loci/shared_loci.tsv").exists()
    return policy_hash


def validate_canonical_lava_decision() -> dict[str, object]:
    """Validate a finalized v3 decision when its production run is present."""
    run_id = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
    run_dir = ROOT / "work/lava-canonical-v3-production" / run_id
    decision_path = run_dir / "canonical_family_decision.json"
    if not decision_path.is_file():
        return {"status": "NOT_FINALIZED", "run_id": run_id}

    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    latest = json.loads((run_dir / "latest_audit.json").read_text(encoding="utf-8"))
    identity = json.loads((run_dir / "run_identity.json").read_text(encoding="utf-8"))
    family_path = OUT / "config/lava_family_canonical_v3.json"
    execution_path = OUT / "config/lava_execution_canonical_v3.json"
    family = json.loads(family_path.read_text(encoding="utf-8"))
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    aggregate_path = run_dir / "results/canonical_family_results.tsv"
    status_counts = decision["status_counts"]
    untested = status_counts.get("NOT_RUN", 0) + status_counts.get("FAILED", 0)
    expected_rows = family["canonical_univariate"]["n_tests"]

    assert decision["analysis_id"] == family["analysis_id"] == latest["analysis_id"]
    assert decision["run_id"] == run_id == latest["run_id"] == identity["run_id"]
    assert decision["run_identity_sha256"] == identity["run_identity_sha256"] == run_id
    assert decision["family_lock_sha256"] == sha256(family_path) == latest["run_identity"]["family_lock_sha256"]
    assert decision["execution_lock_sha256"] == sha256(execution_path) == latest["run_identity"]["execution_lock_sha256"]
    assert decision["canonical_manifest_sha256"] == identity["canonical_manifest_sha256"]
    assert decision["canonical_script_sha256"] == identity["script_sha256"]
    assert execution["worker_count"] == execution["execution_policy"]["worker_count"] == latest["workers"] == 4
    assert decision["planned_loci"] == family["canonical_univariate"]["loci"] == 2495
    assert decision["planned_cells"] == decision["verified_cells"] == expected_rows == 17465
    assert decision["aggregate_rows_verified"] == decision["aggregate"]["rows"] == expected_rows
    assert Path(decision["aggregate"]["path"]).resolve() == aggregate_path.resolve()
    assert decision["aggregate"]["sha256"] == sha256(aggregate_path)
    assert latest["aggregate"] == decision["aggregate"]
    assert latest["audit"]["valid_cells"] == expected_rows
    assert latest["audit"]["intended_cells"] == expected_rows
    assert latest["audit"]["statuses"] == status_counts
    assert latest["audit"]["untested_fraction"] == decision["untested_fraction"]
    assert untested == decision["untested_cells_including_missing"]
    assert decision["maximum_allowed_untested_cells"] == int(
        expected_rows * family["canonical_univariate"]["maximum_untested_fraction"]
    )
    assert decision["qc_pass"] is False
    assert decision["overall_status"] == "FAILED_QC_NOT_PROMOTED"
    assert decision["promotion_permitted"] is False
    assert decision["pairwise_stage_authorized"] is False
    return {
        "status": decision["overall_status"],
        "run_id": run_id,
        "verified_cells": expected_rows,
        "tested": status_counts["TESTED"],
        "not_run": status_counts["NOT_RUN"],
        "failed": status_counts["FAILED"],
        "decision_sha256": sha256(decision_path),
        "aggregate_sha256": decision["aggregate"]["sha256"],
    }


def validate_canonical_lava_not_run_causes() -> dict[str, object]:
    """Verify the reason table covers every cell in the frozen NOT_RUN set."""
    table_path = OUT / "results/lava/canonical_v3_not_run_cells.tsv"
    provenance_path = OUT / "results/lava/canonical_v3_not_run_cells.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    run_id = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
    decision_path = ROOT / "work/lava-canonical-v3-production" / run_id / "canonical_family_decision.json"
    aggregate_path = ROOT / "work/lava-canonical-v3-production" / run_id / "results/canonical_family_results.tsv"
    builder_path = ROOT / "brain6/scripts/build_lava_canonical_v3_not_run_causes.py"
    rows = read_tsv(table_path)
    reasons = {"LOW_LOCAL_H2_UNDERPOWERED": 3564,
               "FEWER_THAN_MIN_K_SHARED_REFERENCE_VARIANTS": 154,
               "FEWER_THAN_MIN_K": 2}
    observed = {}
    for row in rows:
        observed[row["reason"]] = observed.get(row["reason"], 0) + 1
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    assert len(rows) == 3720 and observed == reasons
    assert len({(r["trait_id"], r["locus_id"]) for r in rows}) == len(rows)
    assert all(r["run_id"] == run_id and r["status"] == "NOT_RUN" for r in rows)
    assert all("not evidence for a null association" in r["interpretation"].lower() for r in rows)
    assert provenance["status"] == "COMPLETE_DIAGNOSTIC_ONLY"
    assert provenance["run_id"] == run_id and provenance["not_run_cells"] == len(rows)
    assert provenance["cause_counts"] == reasons
    assert provenance["canonical_decision_sha256"] == sha256(decision_path)
    assert provenance["canonical_aggregate_sha256"] == decision["aggregate"]["sha256"] == sha256(aggregate_path)
    assert provenance["canonical_latest_audit_sha256"] == sha256(decision_path.parent / "latest_audit.json")
    assert provenance["builder_script_sha256"] == sha256(builder_path)
    assert provenance["output_sha256"] == sha256(table_path)
    assert provenance["output_path"] == str(table_path)
    return {"status": provenance["status"], "not_run_cells": len(rows),
            "cause_counts": observed, "output_sha256": provenance["output_sha256"]}


def validate_canonical_lava_cause_summaries() -> dict[str, object]:
    """Validate immutable trait/locus summaries against the complete aggregate."""
    base = OUT / "results/lava"
    provenance_path = base / "canonical_v3_cause_breakdown_v2.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    run_id = "d730debf45266d298401564f3260bdecb14739d1c1f1835a5aebd615c83fa60b"
    run_dir = ROOT / "work/lava-canonical-v3-production" / run_id
    decision_path = run_dir / "canonical_family_decision.json"
    aggregate_path = run_dir / "results/canonical_family_results.tsv"
    latest_path = run_dir / "latest_audit.json"
    cell_path = base / "canonical_v3_not_run_cells.tsv"
    cell_provenance_path = base / "canonical_v3_not_run_cells.provenance.json"
    trait_path = base / "canonical_v3_not_run_by_trait_v1.tsv"
    locus_path = base / "canonical_v3_status_by_locus_v1.tsv"
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    latest = json.loads(latest_path.read_text(encoding="utf-8"))
    trait_expected, locus_expected, expected_summary = summarize_lava_v3_causes(read_tsv(aggregate_path))
    trait_rows, locus_rows = read_tsv(trait_path), read_tsv(locus_path)
    expected_sources = {
        "canonical_family_decision": sha256(decision_path),
        "canonical_aggregate": sha256(aggregate_path),
        "canonical_latest_audit": sha256(latest_path),
        "cell_level_not_run_table": sha256(cell_path),
        "cell_level_not_run_provenance": sha256(cell_provenance_path),
    }
    assert provenance["status"] == "COMPLETE_DIAGNOSTIC_ONLY"
    assert provenance["run_id"] == run_id and provenance["source_hashes"] == expected_sources
    assert provenance["summary"] == expected_summary
    assert provenance["trait_table"] == {
        "path": str(trait_path.relative_to(ROOT)), "rows": 7, "sha256": sha256(trait_path)}
    assert provenance["locus_table"] == {
        "path": str(locus_path.relative_to(ROOT)), "rows": 2495, "sha256": sha256(locus_path)}
    as_tsv_rows = lambda rows: [{key: str(value) for key, value in row.items()} for row in rows]
    assert trait_rows == as_tsv_rows(trait_expected) and locus_rows == as_tsv_rows(locus_expected)
    assert latest["audit"]["statuses"] == decision["status_counts"]
    assert decision["overall_status"] == "FAILED_QC_NOT_PROMOTED"
    assert decision["maximum_allowed_untested_cells"] == 873
    assert all(int(row["planned_cells"]) == 2495 for row in trait_rows)
    assert all(int(row["planned_cells"]) == 7 for row in locus_rows)
    return {"status": provenance["status"], "trait_rows": len(trait_rows),
            "locus_rows": len(locus_rows), "summary": expected_summary,
            "provenance_sha256": sha256(provenance_path)}


def validate_cross_layer_evidence() -> dict[str, object]:
    """Check the five-pair descriptive synthesis against its frozen source rules."""
    table_path = OUT / "results/evidence/brain6_cross_layer_evidence.tsv"
    provenance_path = OUT / "results/evidence/brain6_cross_layer_evidence.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    rules_path = OUT / "config/cross_layer_evidence_rules_v1.json"
    rules_hash_path = rules_path.with_suffix(rules_path.suffix + ".sha256")
    rules_sha = sha256(rules_path)
    assert rules_hash_path.read_text(encoding="utf-8").split()[0] == rules_sha
    rules = json.loads(rules_path.read_text(encoding="utf-8"))
    assert rules["status"] == "FROZEN_BEFORE_CROSS_LAYER_TABLE_BUILD"
    rows = read_tsv(table_path)
    expected_pairs = {
        "insomnia__adhd", "insomnia__mdd", "longsleep__scz",
        "longsleep__bipolar", "longsleep__parkinson",
    }
    assert {row["pair_id"] for row in rows} == expected_pairs and len(rows) == 5
    assert all(row["overall_evidence_category"] == "GENETIC_ONLY" for row in rows)
    assert all(row["canonical_v3_local_rg_family_status"] == "FAILED_QC_NOT_PROMOTED" for row in rows)
    assert all(row["canonical_v3_pairwise_status"] == "NOT_AVAILABLE_PAIRWISE_STAGE_PROHIBITED" for row in rows)
    assert all(row["placo_full_family_status"] == "INCOMPLETE_PROTECTED_TRACK_B_ABSENT" for row in rows)
    adhd = next(row for row in rows if row["pair_id"] == "insomnia__adhd")
    assert adhd["placo_pair_status"] == "PROTECTED_TRACK_B_ABSENT_NOT_RERUN"
    assert adhd["replication_class"] == "DIRECTIONAL_REPLICATION"
    assert all(row["fine_mapping_status"].startswith("NOT_RUN_") for row in rows)
    assert provenance["status"] == "DESCRIPTIVE_PARTIAL_CROSS_LAYER_INTEGRATION"
    assert provenance["rows"] == len(rows) and provenance["evidence_categories"] == {"GENETIC_ONLY": 5}
    assert provenance["output_sha256"] == sha256(table_path)
    assert provenance["output_path"] == str(table_path.relative_to(ROOT))
    assert provenance["builder_script_sha256"] == sha256(ROOT / "brain6/scripts/build_cross_layer_evidence.py")
    assert provenance["inputs"][str(rules_path.relative_to(ROOT))] == rules_sha
    for relpath, expected_hash in provenance["inputs"].items():
        source_path = ROOT / relpath
        assert source_path.is_file() and sha256(source_path) == expected_hash
    return {"status": provenance["status"], "pairs": len(rows),
            "category_counts": provenance["evidence_categories"], "output_sha256": provenance["output_sha256"]}


def validate_partial_placo_figure() -> dict[str, object]:
    """Validate the supplementary PLACO overview without promoting a partial family."""
    figure_dir = OUT / "results/placo/figures"
    figure_id = "figS_partial_placo_candidate_intervals_v1"
    provenance_path = figure_dir / f"{figure_id}.provenance.json"
    outputs = [figure_dir / f"{figure_id}.png", figure_dir / f"{figure_id}.pdf"]
    if not provenance_path.exists() and not any(path.exists() for path in outputs):
        return {"status": "NOT_GENERATED"}
    assert provenance_path.is_file() and all(path.is_file() and path.stat().st_size for path in outputs)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    source_table = OUT / "results/loci/placo_candidate_loci_partial.tsv"
    source_provenance = OUT / "results/loci/placo_candidate_loci_partial.provenance.json"
    family_path = OUT / "config/placo_family_v3/family_lock.json"
    builder_path = ROOT / "brain6/scripts/plot_partial_placo_candidate_loci.py"
    expected_inputs = {
        str(source_table.relative_to(ROOT)): sha256(source_table),
        str(source_provenance.relative_to(ROOT)): sha256(source_provenance),
        str(family_path.relative_to(ROOT)): sha256(family_path),
        str(builder_path.relative_to(ROOT)): sha256(builder_path),
    }
    assert provenance["status"] == "PASS_PARTIAL_FAMILY"
    assert provenance["figure_id"] == figure_id
    assert provenance["inputs"] == expected_inputs
    assert provenance["builder_sha256"] == sha256(builder_path)
    assert provenance["source_table_sha256"] == sha256(source_table)
    assert provenance["source_locus_provenance_sha256"] == sha256(source_provenance)
    assert provenance["family_lock_sha256"] == sha256(family_path)
    assert provenance["candidate_locus_count"] == len(read_tsv(source_table)) == 19
    assert provenance["pairs"] == [
        "insomnia__mdd", "longsleep__scz", "longsleep__bipolar", "longsleep__parkinson",
    ]
    assert "track b" in provenance["scope"].casefold()
    assert "not a complete-family" in provenance["interpretation"].casefold()
    recorded = {Path(record["path"]).name: record for record in provenance["outputs"]}
    assert set(recorded) == {path.name for path in outputs}
    for path in outputs:
        record = recorded[path.name]
        assert record["bytes"] == path.stat().st_size and record["sha256"] == sha256(path)
    return {"status": provenance["status"], "candidate_loci": provenance["candidate_locus_count"],
            "outputs": len(outputs)}


def validate_partial_placo_directions() -> dict[str, object]:
    """Check the four-pair post-selection direction description and provenance."""
    table = OUT / "results/loci/placo_candidate_directions_partial.tsv"
    provenance_path = table.with_suffix(".provenance.json")
    builder = ROOT / "brain6/scripts/build_partial_placo_direction_table.py"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    rows = read_tsv(table)
    pairs = {
        "insomnia__mdd": 124,
        "longsleep__bipolar": 27,
        "longsleep__parkinson": 15,
        "longsleep__scz": 2131,
    }
    assert provenance["status"] == "PASS_PARTIAL_FAMILY_DESCRIPTIVE"
    assert provenance["scope"] == {"pair_ids": list(pairs), "n_pairs": 4,
                                    "n_candidates": 2297, "complete_five_track_family": False}
    assert provenance["sources"]["builder_script"] == {
        "path": str(builder.relative_to(ROOT)), "sha256": sha256(builder)}
    for key in ("candidate_table", "placo_master", "placo_family_lock", "pair_lock", "pair_harmonizer"):
        record = provenance["sources"][key]
        path = ROOT / record["path"]
        assert path.is_file() and sha256(path) == record["sha256"]
    for relpath, expected_hash in provenance["sources"]["pair_variant_files"].items():
        source = ROOT / relpath
        assert source.is_file() and sha256(source) == expected_hash
    assert provenance["output"] == {"path": str(table.relative_to(ROOT)),
                                    "sha256": sha256(table), "rows": len(rows)}
    assert len(rows) == 2297
    observed_pair_counts = {pair: sum(row["pair_id"] == pair for row in rows) for pair in pairs}
    assert observed_pair_counts == pairs
    assert {row["direction_relation"] for row in rows} == {"CONCORDANT", "OPPOSING"}
    assert provenance["direction_counts"]["all"] == {"CONCORDANT": 2296, "OPPOSING": 1}
    opposing = [row for row in rows if row["direction_relation"] == "OPPOSING"]
    assert len(opposing) == 1 and opposing[0]["pair_id"] == "longsleep__bipolar"
    assert opposing[0]["SNP"] == "rs4790841"
    assert "not signed LD" in provenance["interpretation"]
    return {"status": provenance["status"], "candidate_rows": len(rows),
            "direction_counts": provenance["direction_counts"]["all"],
            "output_sha256": provenance["output"]["sha256"]}


def validate_partial_placo_signed_ld() -> dict[str, object]:
    """Verify partial signed LD against its exact inputs and prior r² clumps."""
    table = OUT / "results/loci/placo_candidate_signed_ld_partial.tsv"
    provenance_path = table.with_suffix(".provenance.json")
    builder = ROOT / "brain6/scripts/build_partial_placo_signed_ld.py"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    rows = read_tsv(table)
    assert provenance["status"] == "PASS_SIGNED_LD_WITH_REFERENCE_RANGE_EXCEPTIONS"
    assert provenance["scope"]["complete_five_track_family"] is False
    assert provenance["scope"]["candidate_rows"] == len(rows) == 2297
    assert provenance["scope"]["signed_ld_rows"] == 2269
    assert provenance["scope"]["out_of_range_rows"] == 4
    assert provenance["scope"]["unmatched_rows"] == 24
    assert provenance["sources"]["builder_script"] == {
        "path": str(builder.relative_to(ROOT)), "sha256": sha256(builder)}
    for name, record in provenance["sources"].items():
        if name in {"pair_variant_outputs", "reference_files"}:
            continue
        path = Path(record["path"])
        if not path.is_absolute():
            path = ROOT / path
        assert path.is_file() and sha256(path) == record["sha256"], name
    for pair, expected_hash in provenance["sources"]["pair_variant_outputs"].items():
        record = next(row for row in read_tsv(OUT / "results/placo/placo_master.tsv")
                      if row["pair_id"] == pair)
        source = ROOT / record["output_path"]
        assert source.is_file() and record["output_sha256"] == expected_hash == sha256(source)
    for source_name, expected_hash in provenance["sources"]["reference_files"].items():
        source = Path(source_name)
        assert source.is_file() and sha256(source) == expected_hash
    assert provenance["output"] == {"path": str(table.relative_to(ROOT)),
                                    "sha256": sha256(table), "rows": len(rows)}
    assert len({(row["pair_id"], row["SNP"]) for row in rows}) == len(rows)
    signed = [row for row in rows if row["ld_qc_status"] == "PASS_VALID_SIGNED_R"]
    out_of_range = [row for row in rows if row["ld_qc_status"] == "OUT_OF_RANGE_RAW_LAVA_R"]
    unmatched = [row for row in rows if row["ld_qc_status"] == "NO_EXACT_REFERENCE_OR_LEAD"]
    candidate_map = {(row["pair_id"], row["SNP"]): row
                     for row in read_tsv(OUT / "results/loci/placo_candidate_variants_partial.tsv")}
    assert len(signed) == 2269 and len(out_of_range) == 4 and len(unmatched) == 24
    for row in signed:
        r = float(row["signed_r_ref_A1"])
        r2 = float(row["r2_ref"])
        assert -1 <= r <= 1 and abs(r * r - r2) < 1e-10
        assert row["r2_matches_existing_clump"] == "TRUE"
        candidate = candidate_map[(row["pair_id"], row["SNP"])]
        assert abs(float(candidate["r2_to_lead"]) - r2) < 1e-6
    for row in out_of_range:
        assert abs(float(row["lava_r_raw"])) > 1
        assert row["signed_r_ref_A1"] == row["r2_ref"] == "NA"
        candidate = candidate_map[(row["pair_id"], row["SNP"])]
        assert abs(float(row["lava_r_raw"]) ** 2 - float(candidate["r2_to_lead"])) < 1e-6
        assert row["r2_matches_existing_clump"] == "RAW_R_SQUARED_MATCHES_CLUMP_ONLY"
    assert all(row["lava_r_raw"] == "NA" for row in unmatched)
    assert "Four raw LAVA values are outside" in provenance["interpretation"]
    assert "does not align trait Z signs" in provenance["interpretation"]
    return {"status": provenance["status"], "candidate_rows": len(rows),
            "signed_ld_rows": len(signed), "out_of_range_rows": len(out_of_range),
            "unmatched_rows": len(unmatched),
            "output_sha256": provenance["output"]["sha256"]}


def validate_partial_placo_valid_ld_rebuild() -> dict[str, object]:
    """Verify range-filtered partial clumping without replacing prior candidate outputs."""
    directory = OUT / "results/loci"
    variant_path = directory / "placo_candidate_variants_valid_ld_partial_v2.tsv"
    locus_path = directory / "placo_candidate_loci_valid_ld_partial_v2.tsv"
    exception_path = directory / "placo_candidate_ld_range_exceptions_partial_v2.tsv"
    summary_path = OUT / "results/supplement/table_S22_partial_placo_valid_ld_rebuild.tsv"
    provenance_path = directory / "placo_candidate_valid_ld_partial_v2.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["analysis_id"] == "brain6-placo-partial-candidate-valid-ld-v2"
    assert provenance["status"] == "PASS_PARTIAL_WITH_LD_RANGE_EXCEPTIONS"
    assert provenance["scope"]["complete_five_track_family"] is False
    assert provenance["scope"]["candidate_variants"] == 2297
    assert provenance["scope"]["reference_exact_matches"] == 2273
    assert provenance["scope"]["out_of_range_ld_edges"] == 12317
    assert provenance["scope"]["pair_candidate_exception_rows"] == 12321
    assert provenance["scope"]["candidate_intervals"] == 19
    builder = ROOT / "brain6/scripts/rebuild_partial_placo_clumps_valid_ld.py"
    extractor = ROOT / "brain6/scripts/extract_placo_windowed_ld.R"
    assert provenance["sources"]["builder"] == {
        "path": str(builder.relative_to(ROOT)), "sha256": sha256(builder)}
    assert provenance["sources"]["extractor"]["sha256"] == sha256(extractor)
    for name, record in provenance["sources"].items():
        if name in {"pair_outputs", "reference_files", "preserved_partial_outputs"}:
            continue
        source = Path(record["path"])
        if not source.is_absolute():
            source = ROOT / source
        assert source.is_file() and sha256(source) == record["sha256"], name
    for name, record in provenance["sources"]["preserved_partial_outputs"].items():
        source = ROOT / record["path"]
        assert source.is_file() and sha256(source) == record["sha256"], name
    for pair, expected_hash in provenance["sources"]["pair_outputs"].items():
        master = next(row for row in read_tsv(OUT / "results/placo/placo_master.tsv")
                      if row["pair_id"] == pair)
        source = ROOT / master["output_path"]
        assert source.is_file() and master["output_sha256"] == expected_hash == sha256(source)
    for path_string, expected_hash in provenance["sources"]["reference_files"].items():
        source = Path(path_string)
        assert source.is_file() and sha256(source) == expected_hash
    output_paths = {item["path"]: item for item in provenance["outputs"]}
    for path in (variant_path, locus_path, exception_path, summary_path):
        record = output_paths[str(path.relative_to(ROOT))]
        assert path.is_file() and path.stat().st_size == record["bytes"]
        assert sha256(path) == record["sha256"]
    variants, loci, exceptions, summary = (read_tsv(path) for path in
                                            (variant_path, locus_path, exception_path, summary_path))
    prior_variants = read_tsv(directory / "placo_candidate_variants_partial.tsv")
    prior_identity = {(r["pair_id"], r["SNP"], r["CHR"], r["BP"], r["P_PLACO"])
                      for r in prior_variants}
    new_identity = {(r["pair_id"], r["SNP"], r["CHR"], r["BP"], r["P_PLACO"])
                    for r in variants}
    assert len(variants) == 2297 and new_identity == prior_identity
    assert len(loci) == 19 and sum(int(row["n_lead_signals"]) for row in loci) == 25
    assert sum(row["candidate_status"] == "LEAD" for row in variants) == 25
    assert sum(row["reference_status"] != "EXACT_MATCH" for row in variants) == 24
    invalid_edges = set()
    for row in exceptions:
        assert abs(float(row["lava_r_raw"])) > 1
        assert row["disposition"] == "EXCLUDED_FROM_R2_CLUMPING_RAW_R_OUTSIDE_CORRELATION_RANGE"
        assert abs(int(row["BP2"]) - int(row["BP1"])) <= 500000
        invalid_edges.add((row["pair_id"], *sorted((row["SNP1"], row["SNP2"]))))
    assert len(exceptions) == 12321
    assert len(summary) == 4
    assert sum(int(row["candidate_rows"]) for row in summary) == 2297
    assert sum(int(row["exact_reference_matches"]) for row in summary) == 2273
    assert sum(int(row["unmatched_candidates"]) for row in summary) == 24
    assert sum(int(row["original_lead_variants"]) for row in summary) == 21
    assert sum(int(row["rebuild_valid_ld_lead_variants"]) for row in summary) == 25
    assert sum(int(row["out_of_range_pair_edges"]) for row in summary) == 12321
    for row in variants:
        if row["candidate_status"] != "LD_CLUMPED":
            continue
        edge = (row["pair_id"], *sorted((row["SNP"], row["lead_SNP"])))
        assert edge not in invalid_edges
    old_leads = sum(row["candidate_status"] == "LEAD" for row in prior_variants)
    assert old_leads == 21
    return {"status": provenance["status"], "candidate_rows": len(variants),
            "intervals": len(loci), "old_leads": old_leads,
            "valid_ld_leads": 25, "range_exception_edges": len(exceptions),
            "summary_table_sha256": output_paths[str(summary_path.relative_to(ROOT))]["sha256"],
            "output_sha256": output_paths[str(variant_path.relative_to(ROOT))]["sha256"]}


def validate_partial_placo_candidate_pvalues() -> dict[str, object]:
    directory = OUT / "results/loci"
    result_path = directory / "placo_candidate_pvalue_validation_v1.tsv"
    provenance_path = directory / "placo_candidate_pvalue_validation_v1.provenance.json"
    summary_path = OUT / "results/supplement/table_S39_partial_placo_pvalue_replay.tsv"
    builder = ROOT / "brain6/scripts/validate_partial_placo_pvalues.py"
    helper = ROOT / "brain6/scripts/recompute_partial_placo_candidate_pvalues.R"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    rows, summary = read_tsv(result_path), read_tsv(summary_path)
    assert provenance["status"] == "PASS_PINNED_SOURCE_CANDIDATE_PVALUES_RECOMPUTED"
    assert provenance["output"]["sha256"] == sha256(result_path)
    assert provenance["summary_output"]["sha256"] == sha256(summary_path)
    assert provenance["sources"]["builder"]["sha256"] == sha256(builder)
    assert provenance["sources"]["replay_helper"]["sha256"] == sha256(helper)
    assert provenance["scope"]["full_five_track_family"] is False
    assert provenance["scope"]["candidate_locus_promotion"] is False
    assert len(rows) == 2297 and len(summary) == 4
    assert all(row["validation_status"] == "PASS" for row in rows)
    assert sum(int(row["candidate_variants"]) for row in summary) == 2297
    assert sum(int(row["matching_p_values"]) for row in summary) == 2297
    assert max(float(row["maximum_relative_error"]) for row in summary) < 1e-8
    master = read_tsv(OUT / "results/placo/placo_master.tsv")
    by_output = {row["output_path"]: row for row in master}
    for output_path, expected_hash in provenance["sources"]["pair_outputs"].items():
        source = ROOT / output_path
        assert source.is_file() and sha256(source) == expected_hash
        assert by_output[output_path]["output_sha256"] == expected_hash
    return {"status": provenance["status"], "candidate_variants": len(rows),
            "matching_p_values": sum(row["validation_status"] == "PASS" for row in rows),
            "maximum_relative_error": max(float(row["relative_error"]) for row in rows)}


def validate_partial_placo_signed_effect_ld() -> dict[str, object]:
    """Verify reference-allele effect orientation and signed LD for rebuilt clumps."""
    directory = OUT / "results/loci"
    table = directory / "placo_candidate_signed_effect_ld_partial_v1.tsv"
    summary = OUT / "results/supplement/table_S23_partial_placo_signed_effect_ld.tsv"
    provenance_path = directory / "placo_candidate_signed_effect_ld_partial_v1.provenance.json"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    rows, summary_rows = read_tsv(table), read_tsv(summary)
    builder = ROOT / "brain6/scripts/build_partial_placo_signed_effect_ld.py"
    extractor = ROOT / "brain6/scripts/extract_placo_candidate_lead_ld.R"
    assert provenance["analysis_id"] == "brain6-placo-partial-signed-effect-ld-v1"
    assert provenance["status"] == "PASS_PARTIAL_FAMILY_DESCRIPTIVE_WITH_LD_RANGE_EXCEPTIONS"
    assert provenance["scope"]["complete_five_track_family"] is False
    assert provenance["scope"]["candidate_rows"] == len(rows) == 2297
    assert provenance["sources"]["builder_script"] == {
        "path": str(builder.relative_to(ROOT)), "sha256": sha256(builder)}
    assert provenance["sources"]["extractor_script"] == {
        "path": str(extractor.relative_to(ROOT)), "sha256": sha256(extractor)}
    clump_table = directory / "placo_candidate_variants_valid_ld_partial_v2.tsv"
    clump_provenance = directory / "placo_candidate_valid_ld_partial_v2.provenance.json"
    for key, path in (("candidate_table", clump_table), ("candidate_provenance", clump_provenance)):
        record = provenance["sources"][key]
        assert record["path"] == str(path.relative_to(ROOT)) and record["sha256"] == sha256(path)
    for pair, source in provenance["sources"]["joined_pair_inputs"].items():
        pair_path = Path(source["pair_tsv_path"])
        receipt_path = Path(source["pair_receipt_path"])
        assert pair_path.is_file() and sha256(pair_path) == source["pair_tsv_sha256"]
        assert receipt_path.is_file() and sha256(receipt_path) == source["pair_receipt_sha256"]
        assert pair in {"insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz"}
    ref_record = provenance["sources"]["reference_provenance"]
    ref_prov = Path(ref_record["path"])
    assert ref_prov.is_file() and sha256(ref_prov) == ref_record["sha256"]
    for source_path, expected_hash in provenance["sources"]["reference_files"].items():
        source = Path(source_path)
        assert source.is_file() and sha256(source) == expected_hash
    output_records = {x["path"]: x for x in provenance["outputs"]}
    for path in (table, summary):
        record = output_records[str(path.relative_to(ROOT))]
        assert path.is_file() and sha256(path) == record["sha256"]
        assert record["rows"] == len(read_tsv(path))
    assert len({(r["pair_id"], r["SNP"]) for r in rows}) == len(rows)
    clump_rows = {(r["pair_id"], r["SNP"]): r for r in read_tsv(clump_table)}
    direction_rows = {(r["pair_id"], r["SNP"]): r
                      for r in read_tsv(directory / "placo_candidate_directions_partial.tsv")}
    lead_rows = {(r["pair_id"], r["SNP"]): r for r in rows
                 if r["candidate_status"] == "LEAD" and r["SNP"] == r["lead_SNP"]}
    signed = [r for r in rows if r["ld_qc_status"] == "PASS_VALID_SIGNED_R"]
    invalid = [r for r in rows if r["ld_qc_status"] == "OUT_OF_RANGE_RAW_LAVA_R"]
    unmatched = [r for r in rows if r["ld_qc_status"] == "NO_EXACT_REFERENCE_OR_LEAD"]
    assert len(signed) == 2273 and len(invalid) == 0 and len(unmatched) == 24
    for row in signed:
        r, r2 = float(row["signed_r_reference_A1"]), float(row["r2_reference"])
        assert -1 <= r <= 1 and abs(r*r-r2) < 1e-10
        candidate = clump_rows[(row["pair_id"], row["SNP"])]
        assert row["lead_SNP"] == candidate["lead_SNP"]
        assert abs(float(candidate["r2_to_lead"]) - r2) < 1e-6
        assert row["effect_orientation_status"] == "PASS_REFERENCE_A1_ALIGNED"
        direction = direction_rows[(row["pair_id"], row["SNP"])]
        assert row["direction_relation"] == direction["direction_relation"]
        lead = lead_rows[(row["pair_id"], row["lead_SNP"])]
        assert (int(row["CHR"]), row["lead_reference_A1"], row["lead_reference_A2"]) == (
            int(lead["CHR"]), lead["reference_A1"], lead["reference_A2"])
        assert row["sleep_effect_vs_lead_ld"] == row["disorder_effect_vs_lead_ld"] == "CONSISTENT_WITH_LD_SIGN"
    for row in invalid:
        assert abs(float(row["lava_r_raw"])) > 1
        assert row["signed_r_reference_A1"] == row["r2_reference"] == "NA"
        assert row["sleep_effect_vs_lead_ld"] == row["disorder_effect_vs_lead_ld"] == "NOT_ASSESSED"
    assert {r["pair_id"] for r in summary_rows} == {
        "insomnia__mdd", "longsleep__bipolar", "longsleep__parkinson", "longsleep__scz"}
    assert sum(int(r["candidate_rows"]) for r in summary_rows) == len(rows)
    assert sum(int(r["valid_signed_ld_rows"]) for r in summary_rows) == 2273
    assert sum(int(r["both_effects_consistent_with_ld_sign"]) for r in summary_rows) == 2273
    assert {r["direction_relation"] for r in signed} == {"CONCORDANT", "OPPOSING"}
    assert sum(r["direction_relation"] == "CONCORDANT" for r in signed) == 2272
    assert sum(r["direction_relation"] == "OPPOSING" for r in signed) == 1
    assert "not an" in provenance["interpretation"].casefold()
    return {"status": provenance["status"], "candidate_rows": len(rows),
            "valid_signed_ld_rows": len(signed), "out_of_range_rows": len(invalid),
            "output_sha256": output_records[str(table.relative_to(ROOT))]["sha256"]}


def validate_placo_ld_range_exception_replay() -> dict[str, object]:
    """Verify native LAVA replay of each preserved raw out-of-range LD edge."""
    qc = OUT / "qc"
    provenance_path = qc / "placo_ld_range_exception_replay_v1.provenance.json"
    output_path = qc / "placo_ld_range_exception_replay_v1.tsv"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert provenance["analysis_id"] == "brain6-placo-ld-range-exception-replay-v1"
    assert provenance["scope"]["replay_status"] == "PASS"
    assert provenance["scope"]["unique_edges"] == 12317
    assert provenance["scope"]["replayed_edges"] == 12317
    assert provenance["scope"]["maximum_absolute_replay_delta"] <= 1e-12
    assert provenance["scope"]["minimum_submatrix_eigenvalue"] < 0
    assert provenance["scope"]["reference_matrix_diagonal_error"] <= 1e-6
    assert provenance["scope"]["reference_matrix_symmetry_error"] <= 1e-8
    for name in ("exception_table", "prior_extraction_provenance", "reference_provenance"):
        record = provenance["inputs"][name]
        path = Path(record["path"])
        if not path.is_absolute():
            path = ROOT / path
        assert path.is_file() and sha256(path) == record["sha256"]
    for path_string, record in provenance["inputs"]["reference_files"].items():
        path = Path(path_string)
        assert path.is_file() and sha256(path) == record["sha256"]
    for record in provenance["sources"].values():
        path = ROOT / record["path"]
        assert path.is_file() and sha256(path) == record["sha256"]
    assert output_path.is_file() and sha256(output_path) == provenance["output"]["sha256"]
    rows = read_tsv(output_path)
    assert [int(row["chromosome"]) for row in rows] == provenance["scope"]["chromosomes"]
    assert sum(int(row["replay_matches"]) for row in rows) == 12317
    assert max(float(row["max_abs_replay_delta"]) for row in rows) <= 1e-12
    assert max(float(row["diagonal_max_abs_deviation"]) for row in rows) <= 1e-6
    assert max(float(row["symmetry_max_abs"]) for row in rows) <= 1e-8
    assert min(float(row["minimum_eigenvalue"]) for row in rows) < 0
    return {"status": "PASS_NATIVE_REPLAY_WITH_NON_PSD_SUBMATRICES",
            "replayed_edges": sum(int(row["replay_matches"]) for row in rows),
            "output_sha256": sha256(output_path)}


def validate_placo_ld_factor_normalization() -> dict[str, object]:
    """Validate the diagnostic-only factor reconstruction without promoting LD."""
    qc = OUT / "qc"
    provenance_path = qc / "placo_ld_factor_normalization_v3.provenance.json"
    output_path = qc / "placo_ld_factor_normalization_v3.tsv"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    summary = provenance["summary"]
    assert provenance["analysis_id"] == "brain6-placo-ld-factor-normalization-v1"
    assert summary["unique_edges"] == 12317 and summary["all_within_primary_blocks"] is True
    assert summary["promoted_edges"] == 0 and summary["scaled_abs_r_gt_1"] == 0
    assert summary["max_factor_vs_lava_abs_delta"] <= 1e-12
    assert 0.999 < summary["factor_diagonal_min"] < 1 < summary["factor_diagonal_max"]
    for field in ("exception_table", "native_replay", "replay_provenance", "audit_script"):
        item = provenance["sources"][field]
        path = ROOT / item["path"]
        assert path.is_file() and sha256(path) == item["sha256"]
    for path_string, item in provenance["reference_files"].items():
        path = Path(path_string)
        assert path.is_file() and path.stat().st_size == item["bytes"] and sha256(path) == item["sha256"]
    assert output_path.is_file() and sha256(output_path) == provenance["output"]["sha256"]
    rows = read_tsv(output_path)
    assert len(rows) == 12317
    assert len({(row["CHR"], row["SNP1"], row["SNP2"]) for row in rows}) == len(rows)
    assert all(row["promoted_for_analysis"] == "NO" for row in rows)
    assert all(abs(float(row["lava_r_raw"])) > 1 for row in rows)
    assert all(abs(float(row["factor_product_r"]) - float(row["lava_r_raw"])) <= 1e-12 for row in rows)
    assert all(abs(float(row["unit_diagonal_scaled_r_diagnostic_only"])) <= 1 for row in rows)
    return {"status": "PASS_DIAGNOSTIC_ONLY_NO_PROMOTION", "unique_edges": len(rows),
            "max_factor_delta": summary["max_factor_vs_lava_abs_delta"], "promoted_edges": 0,
            "output_sha256": sha256(output_path)}


def validate_placo_ld_primary_block_psd() -> dict[str, object]:
    """Validate block-level PSD evidence as a non-promoted reference diagnostic."""
    qc = OUT / "qc"
    provenance_path = qc / "placo_ld_primary_block_psd_v2.provenance.json"
    output_path = qc / "placo_ld_primary_block_psd_v2.tsv"
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    summary = provenance["summary"]
    assert provenance["analysis_id"] == "brain6-placo-ld-primary-block-psd-v1"
    assert summary["blocks"] == 19 and summary["promoted_blocks"] == 0
    assert summary["blocks_with_positive_stored_eigenvalues"] == 19
    assert summary["blocks_with_lava_style_negative_min_eigenvalue"] == 19
    assert summary["blocks_normalized_psd_with_numeric_tolerance"] == 19
    for field in ("factor_edge_output", "factor_edge_provenance", "exception_table", "factor_parser", "audit_script"):
        item = provenance["sources"][field]
        path = ROOT / item["path"]
        assert path.is_file() and sha256(path) == item["sha256"]
    for path_string, item in provenance["reference_files"].items():
        path = Path(path_string)
        assert path.is_file() and path.stat().st_size == item["bytes"] and sha256(path) == item["sha256"]
    assert output_path.is_file() and sha256(output_path) == provenance["output"]["sha256"]
    rows = read_tsv(output_path)
    assert len(rows) == 19 and len({(row["CHR"], row["primary_block"]) for row in rows}) == 19
    assert all(int(row["n_snps"]) == 250 and row["promoted_for_analysis"] == "NO" for row in rows)
    assert all(float(row["stored_eigenvalue_min"]) > 0 for row in rows)
    assert all(float(row["lava_unit_diagonal_min_eigenvalue"]) < 0 for row in rows)
    assert all(row["unit_scaled_psd_with_numeric_tolerance"] == "True" for row in rows)
    return {"status": "PASS_DIAGNOSTIC_ONLY_NO_PROMOTION", "blocks": len(rows),
            "lava_style_non_psd_blocks": sum(float(row["lava_unit_diagonal_min_eigenvalue"]) < 0 for row in rows),
            "normalized_psd_blocks": sum(row["unit_scaled_psd_with_numeric_tolerance"] == "True" for row in rows),
            "promoted_blocks": 0, "output_sha256": sha256(output_path)}


def main() -> None:
    pgc_mdd2025_overlap = validate_pgc_mdd2025_cohort_overlap()
    evidence_tier_policy_sha256 = validate_evidence_tier_policy()
    lava_decision = validate_canonical_lava_decision()
    lava_not_run_causes = validate_canonical_lava_not_run_causes()
    lava_cause_summaries = validate_canonical_lava_cause_summaries()
    cross_layer_evidence = validate_cross_layer_evidence()
    partial_placo_figure = validate_partial_placo_figure()
    partial_placo_directions = validate_partial_placo_directions()
    partial_placo_signed_ld = validate_partial_placo_signed_ld()
    partial_placo_valid_ld = validate_partial_placo_valid_ld_rebuild()
    partial_placo_signed_effect_ld = validate_partial_placo_signed_effect_ld()
    placo_ld_exception_replay = validate_placo_ld_range_exception_replay()
    placo_ld_factor_normalization = validate_placo_ld_factor_normalization()
    placo_ld_primary_block_psd = validate_placo_ld_primary_block_psd()
    placo_candidate_pvalues = validate_partial_placo_candidate_pvalues()
    placo_factor_normalized_sensitivity = validate_placo_factor_sensitivity(ROOT)
    brain6_novelty_crosswalk = validate_brain6_novelty_crosswalk()
    brain6_published_rg_context = validate_brain6_published_rg_context(validate=True)
    atlas = read_tsv(ROOT / "results/atlas/trait_pairs.tsv")
    panel = read_tsv(ROOT / "config/analysis_panel.tsv")
    global_rows = read_tsv(OUT / "results/global/brain6_72_locked.tsv")
    source = read_tsv(OUT / "results/global/brain6_72_figure_source.tsv")
    assert len(atlas) == 396 and len(global_rows) == 72 and len(source) == 72
    original = {(r["sleep_trait"], r["non_sleep_trait"]): r for r in atlas}
    assert len(original) == 396
    seen = set()
    n_sig = 0
    for row in global_rows:
        key = row["sleep_trait"], row["brain_disorder"]
        assert key not in seen and key in original
        seen.add(key)
        parent = original[key]
        assert abs(float(row["original_BH_FDR_q"]) - float(parent["global_rg_fdr_all_396"])) < 1e-15
        assert row["significance_under_original_396_family"].lower() == parent["global_rg_primary_significant"].lower()
        assert abs(float(row["rg"]) - float(parent["global_rg"])) < 1e-12
        n_sig += row["significance_under_original_396_family"].lower() == "true"
    expected = {(t["trait_id"], d) for t in panel if t["domain"] == "sleep"
                for d in ("adhd", "mdd", "scz", "bipolar", "alz", "parkinson")}
    assert seen == expected and n_sig == 35
    assert sum(r["brain_disorder"] == "alz" and r["significance_under_original_396_family"].lower() == "true"
               for r in global_rows) == 0
    assert [r for r in global_rows] == [r for r in source]

    profile_config = OUT / "config/brain6_locked_family.yaml"
    profile_map = OUT / "results/global/brain6_72_locked.tsv"
    profile_output = OUT / "results/global/disorder_profile_similarity_spearman.tsv"
    profile_provenance_path = OUT / "results/global/disorder_profile_similarity_spearman.provenance.json"
    profile_provenance = json.loads(profile_provenance_path.read_text(encoding="utf-8"))
    assert profile_output.read_text(encoding="utf-8") == render_table(profile_config, profile_map)
    assert profile_provenance["status"] == "PASS"
    assert profile_provenance["method"] == "Spearman rank correlation across the locked sleep-trait rg profile"
    assert profile_provenance["interpretation"].startswith("descriptive sensitivity only")
    assert profile_provenance["config_sha256"] == sha256(profile_config)
    assert profile_provenance["global_map_sha256"] == sha256(profile_map)
    assert profile_provenance["output_sha256"] == sha256(profile_output)

    profile_loo_plan_path = OUT / "config/brain6_profile_leave_one_out_v1.json"
    profile_loo_plan_checksum = profile_loo_plan_path.with_suffix(profile_loo_plan_path.suffix + ".sha256")
    profile_loo_plan = json.loads(profile_loo_plan_path.read_text(encoding="utf-8"))
    assert sha256(profile_loo_plan_path) == profile_loo_plan_checksum.read_text(encoding="utf-8").split()[0]
    assert profile_loo_plan["status"] == "FROZEN_BEFORE_COMPUTATION"
    assert profile_loo_plan["method"]["inferential_p_values"] is False
    assert profile_loo_plan["method"]["expected_sleep_trait_count"] == 12
    profile_loo_output = OUT / "results/global/disorder_profile_leave_one_sleep_trait_out.tsv"
    profile_loo_provenance_path = profile_loo_output.with_suffix(".provenance.json")
    profile_loo_provenance = json.loads(profile_loo_provenance_path.read_text(encoding="utf-8"))
    profile_loo_rows = read_tsv(profile_loo_output)
    expected_profile_loo_pairs = {
        (row["claim_id"], row["disorder_a"], row["disorder_b"])
        for row in profile_loo_plan["profiles"]
    }
    expected_profile_loo_text = render_profile_loo(profile_loo_plan, profile_config, profile_map)[0]
    expected_profile_loo_rows = list(csv.DictReader(expected_profile_loo_text.splitlines(), delimiter="\t"))
    assert profile_loo_tables_equivalent(profile_loo_rows, expected_profile_loo_rows)
    assert profile_loo_provenance["status"] == "PASS"
    assert profile_loo_provenance["plan_sha256"] == sha256(profile_loo_plan_path)
    assert profile_loo_provenance["locked_family_config_sha256"] == sha256(profile_config)
    assert profile_loo_provenance["global_map_sha256"] == sha256(profile_map)
    assert profile_loo_provenance["script_sha256"] == sha256(ROOT / "brain6/scripts/build_profile_leave_one_out.py")
    assert profile_loo_provenance["output_sha256"] == sha256(profile_loo_output)
    assert profile_loo_provenance["row_count"] == len(profile_loo_rows) == len(expected_profile_loo_pairs) * 12
    assert set(profile_loo_rows[0]) == set(PROFILE_LOO_FIELDS)
    assert not {"p", "p_value", "se", "standard_error", "confidence_interval"} & set(profile_loo_rows[0])
    for claim_id, disorder_a, disorder_b in expected_profile_loo_pairs:
        pair_rows = [row for row in profile_loo_rows if row["claim_id"] == claim_id]
        assert len(pair_rows) == 12
        assert {(row["disorder_a"], row["disorder_b"]) for row in pair_rows} == {(disorder_a, disorder_b)}
        assert len({row["omitted_sleep_trait"] for row in pair_rows}) == 12
        loo_values = [float(row["leave_one_out_r"]) for row in pair_rows]
        assert all(-1 <= value <= 1 for value in loo_values)
        assert abs(min(loo_values) - float(pair_rows[0]["leave_one_out_min_r"])) < 1e-10
        assert abs(max(loo_values) - float(pair_rows[0]["leave_one_out_max_r"])) < 1e-10
        assert abs(max(loo_values) - min(loo_values) - float(pair_rows[0]["leave_one_out_range_width"])) < 1e-10

    tracks = read_tsv(OUT / "config/deep_tracks_v1.tsv")
    pair_lock = json.loads((ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json").read_text())
    lock_hash = (OUT / "config/deep_tracks_v1.tsv.sha256").read_text().split()[0]
    assert sha256(OUT / "config/deep_tracks_v1.tsv") == lock_hash
    assert len(tracks) == 6 and len(pair_lock["pairs"]) == 5
    assert pair_lock["selection_is_post_global_screen"] is True
    assert pair_lock["protected_legacy_track"] == "A/B/CONTROL immutable"
    eligible = {(r["sleep_trait"], r["disease_trait"]) for r in pair_lock["pairs"]}
    selected = {(r["primary_sleep_trait"], r["disorder"]) for r in tracks if r["eligibility"] == "POST_ATLAS_PRIMARY"}
    assert eligible == selected
    assert sum(r["eligibility"] == "NO_ELIGIBLE_PAIR" and r["disorder"] == "alz" for r in tracks) == 1

    figures = list((OUT / "results/global/figures").glob("fig1*.png"))
    figures += list((OUT / "results/global/figures").glob("fig1*.pdf"))
    assert len(figures) == 10 and all(p.stat().st_size > 0 for p in figures)
    figure_provenance_path = OUT / "results/global/figures/figure_build_provenance.json"
    figure_provenance = json.loads(figure_provenance_path.read_text(encoding="utf-8"))
    assert figure_provenance["status"] == "PASS"
    assert figure_provenance["script_sha256"] == sha256(ROOT / "brain6/scripts/build_locked_global.py")
    assert figure_provenance["style_helper_sha256"] == sha256(ROOT / "brain6/scripts/figure_style.py")
    expected_figure_sources = {
        "config/analysis_panel.tsv": sha256(ROOT / "config/analysis_panel.tsv"),
        "results/atlas/traits.tsv": sha256(ROOT / "results/atlas/traits.tsv"),
        "results/atlas/trait_pairs.tsv": sha256(ROOT / "results/atlas/trait_pairs.tsv"),
        "results/analysis/phase1_master_analysis.tsv": sha256(ROOT / "results/analysis/phase1_master_analysis.tsv"),
        "brain6/results/global/brain6_72_locked.tsv": sha256(OUT / "results/global/brain6_72_locked.tsv"),
        "brain6/results/global/brain6_72_figure_source.tsv": sha256(OUT / "results/global/brain6_72_figure_source.tsv"),
        "brain6/results/global/disorder_profile_similarity.tsv": sha256(OUT / "results/global/disorder_profile_similarity.tsv"),
    }
    assert figure_provenance["sources"] == expected_figure_sources
    expected_figure_outputs = {path.name for path in figures}
    assert set(figure_provenance["outputs"]) == expected_figure_outputs
    for filename, record in figure_provenance["outputs"].items():
        path = OUT / "results/global/figures" / filename
        assert record["bytes"] == path.stat().st_size and record["sha256"] == sha256(path)
    replication_figure = OUT / "results/replication/figures/fig_replication_insomnia_adhd.provenance.json"
    replication_receipt = json.loads(replication_figure.read_text(encoding="utf-8"))
    assert replication_receipt["status"] == "PASS"
    assert replication_receipt["pair"] == "insomnia__adhd"
    assert replication_receipt["same_sleep_gwas_used_for_both_estimates"] is True
    for record in replication_receipt["inputs"] + replication_receipt["outputs"]:
        path = ROOT / record["path"]
        assert path.is_file() and path.stat().st_size == record["bytes"]
        assert sha256(path) == record["sha256"]
    assert len(replication_receipt["outputs"]) == 2
    supplement = list((OUT / "results/supplement").glob("table_S*.tsv"))
    s1_table = OUT / "results/supplement/table_S1_gwas_metadata.tsv"
    s1_provenance_path = OUT / "results/supplement/table_S1_gwas_metadata.provenance.json"
    s1_provenance = json.loads(s1_provenance_path.read_text(encoding="utf-8"))
    availability = OUT / "manifests/gwas_external_availability.tsv"
    archive_audit = OUT / "qc/source_archive_audit.tsv"
    dense_qc = OUT / "qc/dense_source_qc_master.tsv"
    assert s1_provenance["status"] == "PASS"
    assert s1_provenance["locked_gwas_master_sha256"] == sha256(OUT / "manifests/gwas_master.tsv")
    assert s1_provenance["locked_dense_input_audit_sha256"] == sha256(OUT / "manifests/locked_dense_input_audit.tsv")
    assert s1_provenance["source_archive_audit_sha256"] == sha256(archive_audit)
    assert s1_provenance["dense_source_qc_master_sha256"] == sha256(dense_qc)
    assert s1_provenance["external_availability_sha256"] == sha256(availability)
    assert s1_provenance["output_table_sha256"] == sha256(s1_table)
    s1_rows = read_tsv(s1_table)
    assert len(s1_rows) == s1_provenance["data_rows"]
    for row in s1_rows:
        if row["trait"] in {"insomnia", "longsleep", "adhd", "mdd", "scz", "bipolar", "parkinson"}:
            assert row["source_archive_status"] in {"VERIFIED_EXACT_ARCHIVE", "VERIFIED_DECOMPRESSED_SOURCE"}
            assert row["raw_file_path"] == row["source_archive_path"]
            assert row["raw_file_SHA256"] == row["source_archive_container_SHA256"]
            assert row["harmonized_input_path"] != row["raw_file_path"]
    software_source = OUT / "manifests/software_versions.tsv"
    software_table = OUT / "results/supplement/table_S17_software_resources.tsv"
    software_provenance = json.loads(
        (OUT / "results/supplement/table_S17_software_resources.provenance.json").read_text(encoding="utf-8")
    )
    assert software_source.is_file() and software_table.read_bytes() == software_source.read_bytes()
    assert software_provenance["status"] == "PASS"
    assert software_provenance["source_sha256"] == sha256(software_source)
    assert software_provenance["output_sha256"] == sha256(software_table)
    assert software_provenance["data_rows"] == len(read_tsv(software_table))
    context_source = OUT / "manifests/published_context_literature.json"
    context_table = OUT / "results/supplement/table_S18_published_rg_context.tsv"
    context_provenance_path = OUT / "results/supplement/table_S18_published_rg_context.provenance.json"
    context_provenance = json.loads(context_provenance_path.read_text(encoding="utf-8"))
    context_rows = read_tsv(context_table)
    assert context_provenance["status"] == "PASS_CONTEXT_ONLY"
    assert context_provenance["source_manifest_sha256"] == sha256(context_source)
    assert context_provenance["discovery_table_sha256"] == sha256(OUT / "results/global/brain6_72_locked.tsv")
    assert context_provenance["output_table_sha256"] == sha256(context_table)
    assert context_provenance["statistical_reanalysis_performed"] is False
    assert context_provenance["replication_admitted"] is False
    assert len(context_rows) == 1
    context_row = context_rows[0]
    assert context_row["classification"] == "LEVEL_D_COMPARABLE_PHENOTYPE_CONTEXT"
    assert context_row["admission_status"] == "CONTEXT_ONLY_NOT_PHASE5_REPLICATION"
    assert abs(float(context_row["published_rg"]) - 0.30) < 1e-12
    assert abs(float(context_row["published_se"]) - 0.05) < 1e-12
    assert abs(float(context_row["published_p"]) - 1.58e-10) < 1e-20
    assert abs(float(context_row["brain6_rg"]) - 0.2889) < 1e-12
    assert abs(float(context_row["published_minus_brain6_rg_point_difference"]) - 0.0111) < 1e-12
    assert context_row["direction_concordant"] == "TRUE"
    assert len(supplement) == 13
    assert OUT / "results/supplement/table_S22_partial_placo_valid_ld_rebuild.tsv" in supplement
    assert OUT / "results/supplement/table_S23_partial_placo_signed_effect_ld.tsv" in supplement
    assert OUT / "results/supplement/table_S24_lava_roundoff_comparison.tsv" in supplement
    assert OUT / "results/supplement/table_S39_partial_placo_pvalue_replay.tsv" in supplement
    context_audit = OUT / "results/lava/lava_pair_context_univariates_v1.tsv"
    context_summary = OUT / "results/supplement/table_S27_lava_pair_context_univariates.tsv"
    context_provenance_path = OUT / "results/lava/lava_pair_context_univariates_v1.provenance.json"
    context_provenance = json.loads(context_provenance_path.read_text(encoding="utf-8"))
    context_summary_rows = read_tsv(context_summary)
    assert context_provenance["status"] == "DESCRIPTIVE_COMPLETE_NO_FAMILY_PROMOTION"
    assert context_provenance["receipt_verified_loci"] == {"baseline": 2495, "roundoff": 2495}
    assert context_provenance["intended_unique_trait_locus_tests"] == 17465
    assert context_provenance["pair_context_rows"] == 36819
    assert context_provenance["output_sha256"] == sha256(context_audit)
    assert context_provenance["summary_sha256"] == sha256(context_summary)
    assert len(read_tsv(context_audit)) == 36819
    assert len(context_summary_rows) == 14
    for run_name, repeated, conflicting, missing in (
        ("baseline", 2804, 2802, 2833), ("roundoff", 2806, 2804, 2832),
    ):
        subset = [r for r in context_summary_rows if r["run"] == run_name]
        assert sum(int(r["repeated_keys"]) for r in subset) == repeated
        assert sum(int(r["repeated_keys_disagree_gt_1e-12"]) for r in subset) == conflicting
        assert sum(int(r["missing_unique_trait_locus"]) for r in subset) == missing
    assert context_summary in supplement
    supplement_text = (OUT / "paper/supplement.md").read_text(encoding="utf-8")
    assert "Table S25" in supplement_text and "canonical_v3_not_run_by_trait_v1.tsv" in supplement_text
    assert "Table S26" in supplement_text and "canonical_v3_status_by_locus_v1.tsv" in supplement_text
    assert "Table S27" in supplement_text and "table_S27_lava_pair_context_univariates.tsv" in supplement_text
    lava_comparison = OUT / "results/lava/lava_baseline_roundoff_pair_locus_comparison_v1.tsv"
    lava_comparison_provenance_path = lava_comparison.with_suffix(".provenance.json")
    lava_comparison_provenance = json.loads(lava_comparison_provenance_path.read_text(encoding="utf-8"))
    lava_comparison_summary = OUT / "results/supplement/table_S24_lava_roundoff_comparison.tsv"
    lava_comparison_summary_rows = read_tsv(lava_comparison_summary)
    assert lava_comparison_provenance["status"] == "DESCRIPTIVE_COMPLETE_FAMILY_QC_UNRESOLVED"
    assert lava_comparison_provenance["planned_pair_locus_slots"] == 12475
    assert lava_comparison_provenance["comparison_table_sha256"] == sha256(lava_comparison)
    assert lava_comparison_provenance["summary_table_sha256"] == sha256(lava_comparison_summary)
    assert len(read_tsv(lava_comparison)) == 12475
    assert len(lava_comparison_summary_rows) == 5
    assert sum(int(row["loci"]) for row in lava_comparison_summary_rows) == 12475
    assert lava_comparison_provenance["both_tested_slots"] == 1
    assert lava_comparison_provenance["environmental_retry"]["status"] == "PASS_RETRY_ARTIFACT"
    assert len(read_tsv(OUT / "results/replication/replication_master.tsv")) == 6
    sensitivity_path = OUT / "results/sensitivity/headline_claim_sensitivity.tsv"
    sensitivity_rows = read_tsv(sensitivity_path)
    required_sensitivity = {
        "claim_id", "evidence_status", "evidence_sources", "alternative_gwas",
        "independent_replication", "sample_overlap_sensitivity",
        "phenotype_definition_sensitivity", "original_396_family_significance",
        "LAVA_UKB_v1_1", "HDL_L_local_rg", "local_h2_thresholds",
        "alternate_LD_reference", "MHC_exclusion", "fine_mapping_stability",
        "coloc_prior_sensitivity", "alternate_gene_mapping",
        "objective_vs_subjective_sleep", "ancestry_sensitivity",
        "profile_leave_one_sleep_trait_out",
        "overall_sensitivity_status",
    }
    assert required_sensitivity.issubset(sensitivity_rows[0])
    assert len({row["claim_id"] for row in sensitivity_rows}) == len(sensitivity_rows)
    placo_master_rows = read_tsv(OUT / "results/placo/placo_master.tsv")
    expected_sensitivity_ids = {
        "global_map_35_of_72", "alzheimer_no_inherited_signals", "insomnia_mdd_largest_global_rg",
        "profile_adhd_mdd", "profile_scz_bipolar", "insomnia_adhd_directional_replication",
        "longsleep_scz_published_context", "placo_partial_candidate_regions",
    } | {f"placo_{row['pair_id']}_headline_variants" for row in placo_master_rows
         if row["stage_status"].startswith("COMPLETE_QC_PASS")}
    assert {row["claim_id"] for row in sensitivity_rows} == expected_sensitivity_ids
    mhc_claim = next(row for row in sensitivity_rows if row["claim_id"] == "placo_partial_candidate_regions")
    assert mhc_claim["MHC_exclusion"].startswith("NO_CANDIDATES_IN_INTERVAL;")
    pair_mhc_claims = [row for row in sensitivity_rows
                       if row["claim_id"].startswith("placo_") and row["claim_id"].endswith("_headline_variants")]
    assert len(pair_mhc_claims) == 4
    assert all(row["MHC_exclusion"].startswith("NO_CANDIDATES_IN_INTERVAL;") for row in pair_mhc_claims)
    mhc_plan_path = OUT / "config/placo_candidate_mhc_exclusion_v1.json"
    mhc_path = OUT / "results/sensitivity/placo_candidate_mhc_exclusion.tsv"
    mhc_provenance_path = OUT / "results/sensitivity/placo_candidate_mhc_exclusion.provenance.json"
    mhc_plan = json.loads(mhc_plan_path.read_text(encoding="utf-8"))
    mhc_provenance = json.loads(mhc_provenance_path.read_text(encoding="utf-8"))
    assert mhc_provenance["status"] == "PASS_PARTIAL_FAMILY_DESCRIPTIVE_COUNTS"
    assert mhc_provenance["plan_sha256"] == sha256(mhc_plan_path)
    assert mhc_provenance["output_sha256"] == sha256(mhc_path)
    mhc_script_path = ROOT / mhc_provenance["script_path"]
    assert mhc_provenance["script_sha256"] == sha256(mhc_script_path)
    assert mhc_plan["coordinate_build"] == "GRCh37/hg19"
    assert mhc_plan["excluded_interval"]["start_inclusive"] == 25_000_000
    assert mhc_plan["excluded_interval"]["end_exclusive"] == 34_000_000
    s16_path = OUT / "results/supplement/table_S16_headline_claim_sensitivity.tsv"
    s16_provenance_path = s16_path.with_suffix(".provenance.json")
    s16_provenance = json.loads(s16_provenance_path.read_text(encoding="utf-8"))
    assert s16_provenance["status"] == "INTERIM_PARTIAL_SENSITIVITY_COVERAGE"
    assert s16_provenance["source_path"] == str(sensitivity_path.relative_to(ROOT))
    assert s16_provenance["source_sha256"] == sha256(sensitivity_path)
    assert s16_provenance["output_path"] == str(s16_path.relative_to(ROOT))
    assert s16_provenance["output_sha256"] == sha256(s16_path)
    assert s16_provenance["script_sha256"] == sha256(ROOT / s16_provenance["script_path"])
    assert s16_path.read_bytes() == sensitivity_path.read_bytes()
    assert s16_provenance["claim_rows"] == len(sensitivity_rows)
    context_sensitivity = next(row for row in sensitivity_rows if row["claim_id"] == "longsleep_scz_published_context")
    assert context_sensitivity["evidence_status"] == "SUPPORTED_PUBLISHED_CONTEXT_ONLY_LEVEL_D"
    assert context_sensitivity["independent_replication"].startswith("LEVEL_D_CONTEXT_ONLY;")
    assert "NOT_REPLICATION" in context_sensitivity["overall_sensitivity_status"]
    assert all(row["overall_sensitivity_status"] != "PASS" for row in sensitivity_rows)
    for claim_id, expected_range in {
        "profile_adhd_mdd": (0.960937241, 0.983377359),
        "profile_scz_bipolar": (0.830436784, 0.917377623),
    }.items():
        profile_claim = next(row for row in sensitivity_rows if row["claim_id"] == claim_id)
        prefix, range_text, _ = profile_claim["profile_leave_one_sleep_trait_out"].split("; ", 2)
        assert prefix == "12 omissions"
        observed_range = tuple(float(value) for value in range_text.removeprefix("Pearson r range ").split("–"))
        assert all(abs(observed - expected) < 0.001 for observed, expected in zip(observed_range, expected_range))
        assert "LOO_RANGES_REPORTED" in profile_claim["overall_sensitivity_status"]
    for row in sensitivity_rows:
        for source in row["evidence_sources"].split(";"):
            path_part, expected_hash = source.split("@sha256:", 1)
            source_path = ROOT / path_part
            assert source_path.is_file() and sha256(source_path) == expected_hash
    placo_locus_summary = {"candidate_locus_status": "NOT_GENERATED",
                           "candidate_pairs": 0, "candidate_variants": 0, "candidate_loci": 0}
    locus_prov_path = OUT / "results/loci/placo_candidate_loci_partial.provenance.json"
    if locus_prov_path.exists():
        provenance = json.loads(locus_prov_path.read_text(encoding="utf-8"))
        assert provenance["status"] in {"PASS_PARTIAL_FAMILY", "PASS"}
        rule_path = OUT / "config/shared_locus_rule_v1.json"
        family_path = OUT / "config/placo_family_v3/family_lock.json"
        assert sha256(rule_path) == provenance["rule_sha256"]
        assert (OUT / "config/shared_locus_rule_v1.json.sha256").read_text().split()[0] == sha256(rule_path)
        assert sha256(family_path) == provenance["placo_family_lock_sha256"]
        for output in provenance["outputs"]:
            path = ROOT / output["path"]
            assert path.is_file() and path.stat().st_size == output["bytes"]
            assert sha256(path) == output["sha256"]
        candidate_variants = read_tsv(OUT / "results/loci/placo_candidate_variants_partial.tsv")
        candidate_loci = read_tsv(OUT / "results/loci/placo_candidate_loci_partial.tsv")
        assert len(candidate_variants) == provenance["n_candidate_variants"]
        assert len(candidate_loci) == provenance["n_candidate_loci"]
        assert len({(r["pair_id"], r["SNP"]) for r in candidate_variants}) == len(candidate_variants)
        assert all(r["candidate_status"] in {"LEAD", "LD_CLUMPED", "NOT_CLUMPED_NO_EXACT_REFERENCE_MATCH"}
                   for r in candidate_variants)
        assert all(r["locus_status"] == "PLACO_ONLY_CANDIDATE_LOCUS" for r in candidate_loci)
        assert {r["pair_id"] for r in candidate_variants}.issubset(set(provenance["pairs_included"]))
        assert sum(int(r["n_candidate_variants"]) for r in candidate_loci) == sum(
            r["candidate_status"] in {"LEAD", "LD_CLUMPED"} for r in candidate_variants)
        placo_locus_summary = {"candidate_locus_status": provenance["status"],
                               "candidate_pairs": len(provenance["pairs_included"]),
                               "candidate_variants": len(candidate_variants),
                               "candidate_loci": len(candidate_loci)}
    placo_audit_path = OUT / "results/placo/placo_v3_pair_qc_validation.json"
    if placo_audit_path.exists():
        placo_audit = json.loads(placo_audit_path.read_text(encoding="utf-8"))
        assert placo_audit["status"] == "PASS_FOUR_PAIR_QC; FIVE_TRACK_FAMILY_INCOMPLETE"
        assert placo_audit["family_lock_sha256"] == sha256(OUT / "config/placo_family_v3/family_lock.json")
        assert placo_audit["q_values_independently_verified"] == sum(int(r["tested_rows"]) for r in placo_master_rows)
        assert placo_audit["full_family_complete"] is False
        assert placo_audit["protected_track_b_published"] is False
        assert set(placo_audit["published_output_sha256"]) == {r["pair_id"] for r in placo_master_rows}
        for row in placo_master_rows:
            assert placo_audit["published_output_sha256"][row["pair_id"]] == row["output_sha256"]
    report = {
        "validation": "PASS",
        "pgc_mdd2025_cohort_overlap": pgc_mdd2025_overlap,
        "global_rows": len(global_rows),
        "inherited_396_family_significant": n_sig,
        "alzheimer_global_significant": 0,
        "deep_track_rows": len(tracks),
        "deep_tracks_sha256": lock_hash,
        "figure_files_png_pdf": len(figures),
        "replication_figure_outputs_hash_verified": len(replication_receipt["outputs"]),
        "supplementary_tables_generated": len(supplement),
        "interim_sensitivity_claims": len(sensitivity_rows),
        "spearman_profile_sensitivity": "PASS",
        "leave_one_sleep_trait_out_profile_sensitivity": "PASS_DESCRIPTIVE_ONLY",
        **placo_locus_summary,
        "final_freeze_created": (OUT / "results/final_v1").exists(),
        "canonical_lava_v3_decision": lava_decision,
        "canonical_lava_v3_not_run_causes": lava_not_run_causes,
        "canonical_lava_v3_cause_summaries": lava_cause_summaries,
        "cross_layer_evidence": cross_layer_evidence,
        "partial_placo_candidate_figure": partial_placo_figure,
        "partial_placo_candidate_directions": partial_placo_directions,
        "partial_placo_signed_ld": partial_placo_signed_ld,
        "partial_placo_valid_ld_rebuild": partial_placo_valid_ld,
        "partial_placo_signed_effect_ld": partial_placo_signed_effect_ld,
        "placo_ld_exception_replay": placo_ld_exception_replay,
        "placo_ld_factor_normalization": placo_ld_factor_normalization,
        "placo_ld_primary_block_psd": placo_ld_primary_block_psd,
        "placo_candidate_pvalue_replay": placo_candidate_pvalues,
        "placo_factor_normalized_sensitivity": placo_factor_normalized_sensitivity,
        "brain6_novelty_crosswalk": brain6_novelty_crosswalk,
        "brain6_published_rg_context": brain6_published_rg_context,
        "evidence_tier_policy": "PASS_PRE_ANNOTATION_LOCK",
        "evidence_tier_policy_sha256": evidence_tier_policy_sha256,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
