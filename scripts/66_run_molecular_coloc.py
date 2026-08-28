#!/usr/bin/env python3
"""Run one checksum-locked trait-molecular SuSiE-RSS/coloc-SuSiE comparison."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


OUTPUTS = {
    "variants.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role", "dataset_type", "SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "MAF", "INFO", "prior_method", "normalized_prior_weight", "PIP", "credible_set_ids", "max_alpha_component", "model_converged", "rss_ld_s", "kriging_allele_switch_outlier"],
    "credible_sets.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role", "signal_id", "component_index", "lead_snp", "lead_pip", "credible_set_size", "credible_set_snps", "requested_coverage", "achieved_coverage", "min_abs_corr", "mean_abs_corr", "median_abs_corr", "cs_log10bf", "model_converged"],
    "colocalization.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset1_id", "dataset2_id", "molecular_feature_id", "tissue_cell_context", "coloc_method", "p1", "p2", "p12", "prior_role", "signal1", "signal2", "hit1", "hit2", "nsnps", "PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "PP_H4_over_PP_H3", "top_shared_variant", "top_shared_variant_PP_H4", "fine_mapping_qc", "analysis_status", "single_signal_fallback_justification", "claim_limit"],
    "shared_variant_posteriors.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "coloc_method", "p12", "signal1", "signal2", "SNP", "SNP_PP_H4"],
    "diagnostics.tsv": ["comparison_id", "pair_id", "locus_id", "comparison_type", "dataset_id", "dataset_role", "variant_count", "model_converged", "niter", "credible_set_count", "max_pip", "rss_ld_s", "kriging_allele_switch_outlier_count", "kriging_allele_switch_outliers", "diagnostic_status"],
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames or [], list(reader)


def probability(value: str, field: str, identity: str) -> None:
    if value in {"", "NA", "Inf"}:
        return
    try:
        parsed = float(value)
    except ValueError as exc:
        fail(f"invalid {field} for {identity}")
        raise AssertionError from exc
    if not math.isfinite(parsed) or not 0 <= parsed <= 1:
        fail(f"invalid {field} for {identity}")


def validate_output(path: Path, fields: list[str], comparison: str) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"molecular coloc engine omitted {path.name}")
    observed, rows = read_tsv(path)
    if observed != fields or any(row["comparison_id"] != comparison for row in rows):
        fail(f"unexpected or out-of-task rows in {path.name}")
    payload = path.read_text(encoding="utf-8", errors="ignore").upper()
    if any(marker in payload for marker in ("SYNTHETIC", "PLACEHOLDER", "FAKE_RESULT")):
        fail(f"non-real marker in {path.name}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("comparison_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--task-dir", default="results/molecular/tasks")
    parser.add_argument("--run-dir", default="results/molecular/coloc_runs")
    parser.add_argument("--policy", default="config/molecular_analysis_policy.json")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    task_path = root / args.task_dir / f"{args.comparison_id}.tsv"
    task_lock_path = root / args.task_dir / f"{args.comparison_id}.lock.json"
    policy_path = root / args.policy
    if not task_path.is_file() or not task_lock_path.is_file():
        fail("checksum-locked molecular task is absent")
    task_lock = json.loads(task_lock_path.read_text(encoding="utf-8"))
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    if (
        task_lock.get("comparison_id") != args.comparison_id
        or task_lock.get("task_sha256") != sha256(task_path)
        or task_lock.get("policy_sha256") != sha256(policy_path)
        or task_lock.get("trait_molecular_results_accessed_before_lock") is not False
    ):
        fail("molecular task differs from its pre-result lock")
    _, tasks = read_tsv(task_path)
    if len(tasks) != 1 or tasks[0]["queue_row_id"] != args.comparison_id:
        fail("molecular task identity is invalid")
    task = tasks[0]
    for path_field, hash_field in (
        ("summary1_path", "summary1_sha256"), ("summary2_path", "summary2_sha256"),
        ("ld_path", "ld_sha256"), ("ld_variant_order_path", "ld_variant_order_sha256"),
    ):
        path = Path(task[path_field])
        if not path.is_absolute(): path = root / path
        if not path.is_file() or sha256(path) != task[hash_field]:
            fail(f"locked molecular task input differs: {path_field}")
    engine = root / policy["generic_engine"]["path"]
    if not engine.is_file() or sha256(engine) != policy["generic_engine"]["sha256"]:
        fail("shared SuSiE/coloc engine differs from exact pin")
    final = root / args.run_dir / args.comparison_id
    if final.exists():
        fail("molecular coloc result already exists; refusing overwrite")
    if not args.execute:
        print("Task is locked. Re-run with --execute to perform the real trait-molecular analysis.")
        return 0
    work_root = root / "work/molecular_coloc_runs"
    work_root.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=args.comparison_id + ".", dir=work_root))
    try:
        output_paths = {name: staging / name for name in OUTPUTS}
        p12 = ",".join(format(value, ".12g") for value in policy["colocalization"]["p12_sensitivity"])
        command = [
            str(root / ".r-env/bin/Rscript"), str(engine), str(task_path), args.comparison_id,
            str(output_paths["variants.tsv"]), str(output_paths["credible_sets.tsv"]),
            str(output_paths["colocalization.tsv"]), str(output_paths["shared_variant_posteriors.tsv"]),
            str(output_paths["diagnostics.tsv"]), str(policy["fine_mapping"]["maximum_causal_signals"]),
            str(policy["fine_mapping"]["credible_set_coverage"]),
            str(policy["fine_mapping"]["minimum_absolute_correlation"]),
            str(policy["fine_mapping"]["maximum_iterations"]), str(policy["colocalization"]["p1"]),
            str(policy["colocalization"]["p2"]), p12,
        ]
        result = subprocess.run(command, cwd=root, capture_output=True, text=True)
        if result.returncode:
            fail(f"molecular SuSiE/coloc engine failed: {(result.stdout + result.stderr).strip()}")
        outputs = {name: validate_output(path, OUTPUTS[name], args.comparison_id) for name, path in output_paths.items()}
        diagnostics = outputs["diagnostics.tsv"]
        if len(diagnostics) != 2 or {row["dataset_id"] for row in diagnostics} != {task["dataset1_id"], task["dataset2_id"]}:
            fail("molecular diagnostics do not cover both locked datasets")
        for row in outputs["variants.tsv"]:
            probability(row["PIP"], "PIP", row["SNP"]); probability(row["normalized_prior_weight"], "normalized_prior_weight", row["SNP"])
        coloc_rows = outputs["colocalization.tsv"]
        if coloc_rows and {float(row["p12"]) for row in coloc_rows} != {float(value) for value in policy["colocalization"]["p12_sensitivity"]}:
            fail("molecular colocalization output does not cover exact prior grid")
        for row in coloc_rows:
            for field in ("PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "top_shared_variant_PP_H4"):
                probability(row[field], field, args.comparison_id)
        provenance = {
            "schema_version": "atlas-v1.0-molecular-coloc-run.1", "comparison_id": args.comparison_id,
            "task_sha256": sha256(task_path), "task_lock_sha256": sha256(task_lock_path),
            "policy_sha256": sha256(policy_path), "engine_sha256": sha256(engine),
            "stdout": result.stdout.strip(), "stderr": result.stderr.strip(),
            "outputs": {name: {"rows": len(outputs[name]), "sha256": sha256(path)} for name, path in output_paths.items()},
            "claim_limit": policy["claim_limit"],
        }
        (staging / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        final.parent.mkdir(parents=True, exist_ok=True); os.replace(staging, final)
        print(f"MOLECULAR_COLOC_RUN_OK comparison={args.comparison_id} variant_rows={len(outputs['variants.tsv'])} coloc_rows={len(coloc_rows)}")
        return 0
    finally:
        if staging.exists(): shutil.rmtree(staging)


if __name__ == "__main__":
    raise SystemExit(main())
