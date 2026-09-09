#!/usr/bin/env python3
"""Shared immutable-contract checks for interpretation, robustness, and release."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


INTERPRETATION_SCRIPTS = (
    "scripts/downstream_contract.py", "scripts/liftover_chain.py",
    "scripts/pathway_sources.py", "scripts/catlas_policy.py",
    "scripts/74_interpretation_preflight.py", "scripts/75_prepare_interpretation_tasks.py",
    "scripts/76_run_interpretation_task.py", "scripts/76_record_interpretation_task.py",
    "scripts/77_collate_interpretation.py", "scripts/78_build_atlas_edges.py",
    "scripts/82_run_regulatory_task.py", "scripts/83_run_motif_task.py",
    "scripts/84_prepare_abc_overlap_cache.py", "scripts/85_run_abc_task.py",
    "scripts/86_prepare_pchic_overlap_cache.py", "scripts/87_run_pchic_task.py",
    "scripts/88_materialize_fuma_resources.py", "scripts/89_prepare_magma_annotation.py",
    "scripts/90_prepare_magma_gene_results.py", "scripts/91_run_fuma_scrna_task.py",
    "scripts/92_run_pathway_task.py", "scripts/93_setup_causal_runtime.sh",
    "scripts/94_prepare_causal_reference.py", "scripts/95_prepare_catlas_reference.py",
    "scripts/96_prepare_catlas_trait_cache.py", "scripts/97_run_catlas_task.py",
    "scripts/98_prepare_ldsc_seg_gtex_source.py", "scripts/98_prepare_ldsc_seg_reference.py",
    "scripts/98_prepare_ldsc_seg_trait.py", "scripts/98_run_ldsc_seg_task.py",
)
ROBUSTNESS_SCRIPTS = (
    "scripts/downstream_contract.py", "scripts/79_prepare_robustness_tasks.py",
    "scripts/80_run_robustness_task.py", "scripts/80_record_robustness_task.py",
    "scripts/81_collate_robustness.py", "scripts/53_validate_robustness.py",
)
RELEASE_SCRIPTS = (
    "scripts/downstream_contract.py", "scripts/52_validate_integrated_atlas.py",
    "scripts/53_validate_robustness.py", "scripts/54_build_release.py",
    "scripts/55_validate_release.py", "scripts/99_atlas_acceptance.py",
)


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty checksum-bound artifact: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty tabular artifact: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def script_hashes(root: Path, family: str) -> dict[str, str]:
    paths = {
        "interpretation": INTERPRETATION_SCRIPTS,
        "robustness": ROBUSTNESS_SCRIPTS,
        "release": RELEASE_SCRIPTS,
    }.get(family)
    if paths is None:
        fail(f"unknown downstream script family: {family}")
    return {relative: sha256(root / relative) for relative in paths}


def validate_interpretation_preflight(
    root: Path, policy_path: Path, downstream_path: Path, report_path: Path,
) -> dict[str, object]:
    try:
        report = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"interpretation preflight is unreadable: {exc}")
    readiness_path = root / "results/tables/interpretation_source_readiness.tsv"
    registry_path = root / "config/interpretation_source_registry.tsv"
    references_path = root / "config/interpretation_method_references.tsv"
    if (
        report.get("schema_version") != "atlas-v1.0-interpretation.1"
        or report.get("analysis_id") != "atlas-v1.0-interpretation"
        or report.get("policy_sha256") != sha256(policy_path)
        or report.get("downstream_policy_sha256") != sha256(downstream_path)
        or report.get("source_registry_sha256") != sha256(registry_path)
        or report.get("method_references_sha256") != sha256(references_path)
        or report.get("readiness_sha256") != sha256(readiness_path)
        or report.get("script_sha256") != script_hashes(root, "interpretation")
        or report.get("code_ready") is not True
        or report.get("production_ready") is not True
    ):
        fail("interpretation preflight is not production-ready or has drifted")
    upstream = report.get("upstream", {})
    for relative, record in upstream.items():
        path = root / relative
        if not isinstance(record, dict) or record.get("sha256") != sha256(path):
            fail(f"interpretation upstream differs from preflight: {relative}")
    return report


def validate_interpretation_manifest(
    root: Path, manifest_path: Path, lock_path: Path, policy_path: Path,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    fields, tasks = read_tsv(manifest_path)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"interpretation task lock is unreadable: {exc}")
    ids = [row["task_id"] for row in tasks]
    if (
        fields != policy["task_manifest_fields"]
        or lock.get("schema_version") != policy["schema_version"]
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("task_manifest_sha256") != sha256(manifest_path)
        or lock.get("task_ids_in_locked_order") != ids
        or lock.get("task_count") != len(tasks) or len(ids) != len(set(ids))
        or lock.get("interpretation_results_accessed_before_task_lock") is not False
        or lock.get("script_sha256") != script_hashes(root, "interpretation")
    ):
        fail("interpretation task family differs from its immutable pre-result lock")
    return tasks, lock


def validate_robustness_manifest(
    root: Path, manifest_path: Path, lock_path: Path, policy_path: Path,
    downstream_path: Path,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    fields, tasks = read_tsv(manifest_path)
    try:
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        fail(f"robustness task lock is unreadable: {exc}")
    ids = [row["task_id"] for row in tasks]
    if (
        fields != policy["robustness"]["task_manifest_fields"]
        or lock.get("schema_version") != policy["schema_version"]
        or lock.get("policy_sha256") != sha256(policy_path)
        or lock.get("downstream_policy_sha256") != sha256(downstream_path)
        or lock.get("task_manifest_sha256") != sha256(manifest_path)
        or lock.get("task_ids_in_locked_order") != ids
        or lock.get("task_count") != len(tasks) or len(ids) != len(set(ids))
        or lock.get("robustness_results_accessed_before_task_lock") is not False
        or lock.get("script_sha256") != script_hashes(root, "robustness")
    ):
        fail("robustness task family differs from its immutable pre-result lock")
    return tasks, lock
