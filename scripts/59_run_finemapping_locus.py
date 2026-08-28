#!/usr/bin/env python3
"""Run one checksum-locked SuSiE-RSS/coloc-SuSiE trait-trait comparison."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path

import fine_mapping_contract

SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
OUTPUTS = fine_mapping_contract.ENGINE_OUTPUT_FIELDS


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
        fail(f"fine-mapping engine omitted {path.name}")
    observed, rows = read_tsv(path)
    if observed != fields:
        fail(f"unexpected {path.name} header")
    if any(row["comparison_id"] != comparison for row in rows):
        fail(f"out-of-task row in {path.name}")
    payload = path.read_text(encoding="utf-8", errors="ignore").upper()
    if any(marker in payload for marker in ("SYNTHETIC", "PLACEHOLDER", "FAKE_RESULT")):
        fail(f"non-real marker in {path.name}")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shared_locus_id")
    parser.add_argument("--root", default=".")
    parser.add_argument("--manifest", default="results/tables/fine_mapping_locus_manifest.tsv")
    parser.add_argument("--manifest-lock", default="results/tables/fine_mapping_locus_manifest.lock.json")
    parser.add_argument("--policy", default="config/fine_mapping_analysis_policy.json")
    parser.add_argument("--preflight", default="results/tables/fine_mapping_preflight.json")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not SAFE_ID.fullmatch(args.shared_locus_id):
        fail("unsafe shared-locus identifier")
    root = Path(args.root).resolve()
    manifest_path = root / args.manifest
    manifest_lock_path = root / args.manifest_lock
    policy_path, policy = fine_mapping_contract.load_policy(root, args.policy)
    preflight_path = root / args.preflight
    fine_mapping_contract.validate_preflight(root, policy_path, policy, preflight_path)
    manifest, manifest_lock = fine_mapping_contract.validate_manifest(
        root, manifest_path, manifest_lock_path, policy_path, policy, preflight_path,
    )
    selected = [row for row in manifest if row["shared_locus_id"] == args.shared_locus_id]
    if len(selected) != 1:
        fail("shared_locus_id does not identify exactly one locked locus")
    locus = selected[0]
    task_path = root / locus["task_path"]
    task_lock_path = task_path.with_suffix(".lock.json")
    if not task_path.is_file() or not task_lock_path.is_file():
        fail("checksum-locked fine-mapping task is absent")
    task_lock = json.loads(task_lock_path.read_text(encoding="utf-8"))
    if (
        task_lock.get("schema_version") != "atlas-v1.0-finemapping-task.2"
        or task_lock.get("comparison_id") != locus["comparison_id"]
        or task_lock.get("shared_locus_id") != locus["shared_locus_id"]
        or task_lock.get("task_sha256") != sha256(task_path)
        or task_lock.get("manifest_sha256") != sha256(manifest_path)
        or task_lock.get("manifest_lock_sha256") != sha256(manifest_lock_path)
        or task_lock.get("policy_sha256") != sha256(policy_path)
        or task_lock.get("preflight_sha256") != sha256(preflight_path)
        or task_lock.get("script_sha256") != fine_mapping_contract.script_hashes(root)
        or task_lock.get("results_accessed_before_lock") is not False
    ):
        fail("fine-mapping task differs from its pre-result lock")
    task_fields, task_rows = read_tsv(task_path)
    del task_fields
    if len(task_rows) != 1 or task_rows[0]["queue_row_id"] != locus["comparison_id"]:
        fail("fine-mapping task identity differs from the locus manifest")
    task = task_rows[0]
    for path_field, hash_field in (
        ("summary1_path", "summary1_sha256"), ("summary2_path", "summary2_sha256"),
        ("ld_path", "ld_sha256"), ("ld_variant_order_path", "ld_variant_order_sha256"),
    ):
        relative = fine_mapping_contract.safe_relative(task[path_field], path_field)
        path = root / relative
        if (
            not path.is_file() or sha256(path) != task[hash_field]
            or task_lock.get("outputs", {}).get(str(relative)) != task[hash_field]
        ):
            fail(f"locked task input differs: {path_field}")
    source_records = task_lock.get("source_inputs", {})
    if set(source_records) != {"sleep", "non_sleep"}:
        fail("fine-mapping task omits dense source-input provenance")
    for role, record in source_records.items():
        source_path = root / fine_mapping_contract.safe_relative(
            str(record.get("path", "")), f"{role} dense source",
        )
        if (
            not source_path.is_file() or source_path.stat().st_size != int(record.get("bytes", 0))
            or sha256(source_path) != record.get("sha256")
        ):
            fail(f"fine-mapping dense source input drifted: {role}")
    reference_provenance_path, _ = fine_mapping_contract.validate_lava_reference(
        root, policy, rehash_payloads=False,
    )
    bcor_path = root / fine_mapping_contract.safe_relative(
        str(task_lock.get("reference_bcor_path", "")), "reference bcor",
    )
    if (
        task_lock.get("reference_provenance_sha256") != sha256(reference_provenance_path)
        or task_lock.get("reference_content_verified_by_preflight") is not True
        or not bcor_path.is_file()
        or bcor_path.stat().st_size != int(task_lock.get("reference_bcor_bytes", 0))
    ):
        fail("fine-mapping task no longer binds the sealed signed-LD reference")
    engine = root / policy["shared_engine"]["path"]
    if not engine.is_file() or sha256(engine) != policy["shared_engine"]["sha256"]:
        fail("shared SuSiE/coloc engine differs from its exact pin")
    final = root / "results/fine_mapping/runs" / args.shared_locus_id
    staging = root / "work/fine_mapping_runs" / f"{args.shared_locus_id}.building"
    if final.exists() or staging.exists():
        fail("fine-mapping output/staging directory already exists; results are never overwritten")
    if not args.execute:
        print("Task is locked. Re-run with --execute to perform the real SuSiE-RSS/coloc analysis.")
        return 0
    staging.mkdir(parents=True)
    output_paths = {name: staging / name for name in OUTPUTS}
    p12 = ",".join(format(value, ".12g") for value in policy["colocalization"]["p12_sensitivity"])
    command = [
        str(root / ".r-env/bin/Rscript"), str(engine), str(task_path), task["queue_row_id"],
        str(output_paths["variants.tsv"]), str(output_paths["credible_sets.tsv"]),
        str(output_paths["colocalization.tsv"]), str(output_paths["shared_variant_posteriors.tsv"]),
        str(output_paths["diagnostics.tsv"]), str(policy["fine_mapping"]["maximum_causal_signals"]),
        str(policy["fine_mapping"]["credible_set_coverage"]),
        str(policy["fine_mapping"]["minimum_absolute_correlation"]),
        str(policy["fine_mapping"]["maximum_iterations"]), str(policy["colocalization"]["p1"]),
        str(policy["colocalization"]["p2"]), p12,
    ]
    result = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
    if result.returncode:
        fail(f"SuSiE/coloc engine failed: {(result.stdout + result.stderr).strip()}")
    outputs = {name: validate_output(path, OUTPUTS[name], task["queue_row_id"]) for name, path in output_paths.items()}
    diagnostics = outputs["diagnostics.tsv"]
    if len(diagnostics) != 2 or {row["dataset_id"] for row in diagnostics} != {task["dataset1_id"], task["dataset2_id"]}:
        fail("fine-mapping diagnostics do not cover both locked traits")
    for row in outputs["variants.tsv"]:
        probability(row["PIP"], "PIP", row["SNP"])
        probability(row["normalized_prior_weight"], "normalized_prior_weight", row["SNP"])
    coloc_rows = outputs["colocalization.tsv"]
    observed_p12 = {float(row["p12"]) for row in coloc_rows}
    if observed_p12 != {float(value) for value in policy["colocalization"]["p12_sensitivity"]}:
        fail("colocalization output does not cover the exact prior-sensitivity grid")
    for row in coloc_rows:
        for field in ("PP_H0", "PP_H1", "PP_H2", "PP_H3", "PP_H4", "top_shared_variant_PP_H4"):
            probability(row[field], field, task["queue_row_id"])
    provenance = {
        "schema_version": "atlas-v1.0-finemapping-run.2", "comparison_id": task["queue_row_id"],
        "task_sha256": sha256(task_path), "task_lock_sha256": sha256(task_lock_path),
        "manifest_sha256": sha256(manifest_path), "manifest_lock_sha256": sha256(manifest_lock_path),
        "policy_sha256": sha256(policy_path), "engine_sha256": sha256(engine),
        "preflight_sha256": sha256(preflight_path),
        "reference_provenance_sha256": sha256(reference_provenance_path),
        "script_sha256": fine_mapping_contract.script_hashes(root),
        "stdout": result.stdout.strip(), "stderr": result.stderr.strip(),
        "outputs": {name: {"rows": len(outputs[name]), "sha256": sha256(path)} for name, path in output_paths.items()},
        "claim_limit": policy["claim_limit"],
    }
    (staging / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    final.parent.mkdir(parents=True, exist_ok=True)
    os.replace(staging, final)
    print(
        f"FINEMAPPING_LOCUS_RUN_OK locus={args.shared_locus_id} "
        f"variant_rows={len(outputs['variants.tsv'])} coloc_rows={len(coloc_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
