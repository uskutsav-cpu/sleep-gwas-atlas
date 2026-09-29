#!/usr/bin/env python3
"""Retry one receipt-verified environmental LAVA locus failure in an isolated path.

The retry is a separate provenance artifact. It never replaces the original
roundoff locus receipt or contributes to that run's frozen family audit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
import sys

sys.path.insert(0, str(ROOT / "extensions/brain6"))
sys.path.insert(0, str(ROOT / "brain6/scripts"))
from brain6.io import read_json
from run_lava_family_roundoff_v1 import (
    EXPECTED_PAIRS,
    rows,
    run_locus_unit,
    sha256,
    validate_inputs,
    verify_unit,
)

DEFAULT_MANIFEST = ROOT / "brain6/manifests/lava_roundoff_v1_run.json"
DEFAULT_TEMP = Path("/private/tmp/brain6-lava-roundoff-tmp")


def environmental_failure_reason(pair_row: dict[str, str]) -> bool:
    return (pair_row.get("status") == "FAILED" and
            "No write permission for directory:" in pair_row.get("reason", ""))


def validate_retry_id(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}", value):
        raise ValueError("retry id must be 1–64 letters, digits, underscores, or hyphens")
    return value


def failed_attempt_paths(retry_dir: Path, locus_id: str) -> list[Path]:
    """Find only this retry's failed attempts under the coordinator's layout."""
    return sorted((retry_dir.parent / "failed_attempts").glob(f"{locus_id}-*"))


def immutable_json(path: Path, record: dict[str, Any]) -> None:
    data = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if path.exists():
        if path.read_bytes() != data:
            raise FileExistsError(f"Refusing to overwrite retry provenance: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.rename(temp, path)
    except BaseException:
        Path(temp).unlink(missing_ok=True)
        raise


def identity_hash(identity: dict[str, Any]) -> str:
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def retry(manifest_path: Path, locus_id: str, pair_id: str, retry_id: str,
          temp_dir: Path) -> Path:
    retry_id = validate_retry_id(retry_id)
    manifest = read_json(manifest_path)
    if manifest.get("run_label") != "roundoff-stability follow-up v1":
        raise ValueError("Manifest is not the protected roundoff v1 run")
    run_id = str(manifest["run_id"])
    run_root = Path(manifest["output_root"]).resolve()
    if run_root.name != run_id:
        raise ValueError("Roundoff output directory does not match its frozen run id")
    run_lock = run_root.parent / f".{run_id}.lock"
    if run_lock.exists():
        raise RuntimeError(f"Roundoff coordinator lock is present; do not run a concurrent retry: {run_lock}")

    family_path = ROOT / manifest["family_lock_path"]
    execution_path = ROOT / manifest["execution_lock_path"]
    input_root = Path(manifest["input_root"]).resolve()
    reference = ROOT / "ref/lava/ukb_v1.1"
    family, execution, pair_lock, loci, coverage_map = validate_inputs(
        family_path, execution_path, input_root, reference)
    if tuple(family["pairs"]) != EXPECTED_PAIRS:
        raise ValueError("Retry pair family differs from the frozen five-pair family")

    run_identity_path = run_root / "run_identity.json"
    run_identity = read_json(run_identity_path)
    if identity_hash(run_identity) != run_id:
        raise ValueError("Run identity does not reproduce the immutable roundoff run id")
    if sha256(family_path) != run_identity["family_lock_sha256"]:
        raise ValueError("Family lock differs from the roundoff run identity")
    if sha256(execution_path) != run_identity["execution_lock_sha256"]:
        raise ValueError("Execution lock differs from the roundoff run identity")
    if sha256(input_root / "provenance.json") != run_identity["materialized_inputs_sha256"]:
        raise ValueError("Materialized LAVA inputs differ from the roundoff run identity")
    materialized = read_json(input_root / "provenance.json")
    locus_input_records = [record for record in materialized["records"]
                           if str(record.get("locus_id")) == locus_id]
    locus_input_hashes = {}
    for record in locus_input_records:
        input_path = Path(record["path"])
        if input_path.stat().st_size != int(record["bytes"]) or sha256(input_path) != record["sha256"]:
            raise ValueError(f"Target-locus input changed since materialization: {input_path}")
        key = record.get("trait", "input_info")
        locus_input_hashes[key] = record["sha256"]
    if set(locus_input_hashes) != set(family["trait_ids"]) | {"input_info"}:
        raise ValueError("Target locus does not have the exact seven-trait-plus-input-info inventory")

    if locus_id not in {str(row["LOC"]) for row in loci} or pair_id not in EXPECTED_PAIRS:
        raise ValueError("Requested retry slot is outside the frozen LAVA family")
    original_unit = run_root / "loci" / locus_id
    original = verify_unit(original_unit, locus_id, set(EXPECTED_PAIRS), run_id)
    if original is None:
        raise ValueError(f"Original locus receipt is missing or incomplete: {original_unit}")
    original_receipt_path = original_unit / "receipt.json"
    original_receipt_sha = sha256(original_receipt_path)
    original_results_path = original_unit / "worker_output/pair_results.tsv"
    original_results_sha = sha256(original_results_path)
    original_rows = [r for r in rows(original_results_path) if r["pair_id"] == pair_id]
    if len(original_rows) != 1 or not environmental_failure_reason(original_rows[0]):
        raise ValueError("The original slot is not a receipt-verified R temp-directory failure")
    original_reason = original_rows[0]["reason"]

    source_config_path = original_unit / "worker_config.json"
    source_config = read_json(source_config_path)
    if sha256(source_config_path) != original["worker_config_sha256"]:
        raise ValueError("Original worker configuration is not receipt-verified")
    if source_config.get("runtime_validator_passed") is not True:
        raise ValueError("Original run did not pass the pinned runtime validation")
    worker = ROOT / "brain6/scripts/run_lava_family_locus_roundoff_v1.R"
    if sha256(worker) != run_identity["worker_sha256"]:
        raise ValueError("Roundoff worker source changed after the original run")
    runner = ROOT / "brain6/scripts/run_lava_family_roundoff_v1.py"
    if sha256(runner) != run_identity["orchestrator_script_sha256"]:
        raise ValueError("Roundoff coordinator source changed after the original run")
    validator = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
    if sha256(validator) != run_identity["runtime_validator_sha256"]:
        raise ValueError("Pinned runtime-validator source changed after the original run")
    rscript = Path(source_config["rscript"]).resolve()
    if not rscript.is_file() or sha256(rscript) != run_identity["rscript_binary_sha256"]:
        raise ValueError("Pinned Rscript binary differs from the original run")
    runtime = subprocess.run([str(rscript), str(validator)], cwd=ROOT,
                             capture_output=True, text=True, check=False)
    runtime_sha = hashlib.sha256(runtime.stdout.encode("utf-8")).hexdigest()
    if (runtime.returncode != 0 or
            runtime.stdout.count("TRACK_B_LAVA_RUNTIME\tstatus=PASS\t") != 1 or
            runtime_sha != run_identity["runtime_validation_output_sha256"]):
        raise RuntimeError("Pinned R/LAVA runtime no longer matches the original run preflight")

    temp_dir.mkdir(parents=True, exist_ok=True)
    probe = temp_dir / f".brain6-write-check-{os.getpid()}"
    try:
        probe.write_text("ok\n", encoding="utf-8")
    finally:
        probe.unlink(missing_ok=True)

    retry_dir = run_root / "targeted_retries" / retry_id
    unit_target = retry_dir / f"locus_{locus_id}"
    retry_record_path = retry_dir / "retry_provenance.json"
    if retry_record_path.is_file():
        record = read_json(retry_record_path)
        existing = verify_unit(unit_target, locus_id, set(EXPECTED_PAIRS), run_id)
        if record.get("original_receipt_sha256") != original_receipt_sha or not existing:
            raise ValueError("Existing retry artifact does not match the preserved original receipt")
        return retry_dir
    if retry_dir.exists():
        raise FileExistsError(f"Retry directory already exists without a valid provenance record: {retry_dir}")

    locus = next(row for row in loci if str(row["LOC"]) == locus_id)
    temp_values = {key: os.environ.get(key) for key in ("TMPDIR", "TMP", "TEMP")}
    os.environ.update({key: str(temp_dir) for key in ("TMPDIR", "TMP", "TEMP")})
    try:
        retry_dir.mkdir(parents=True)
        execution_policy = {
            "locus_processing": execution["locus_processing"],
            "univariate": execution["univariate"],
            "bivariate": execution["bivariate"],
        }
        exit_code = run_locus_unit(
            unit_target, locus, family, execution, pair_lock, input_root, family_path,
            execution_path, run_id, run_identity, worker,
            validator, rscript,
            execution_policy, coverage_map,
        )
    finally:
        for key, value in temp_values.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

    if exit_code != 0:
        failures = failed_attempt_paths(retry_dir, locus_id)
        immutable_json(retry_record_path, {
            "status": "RETRY_FAILED_PRESERVED",
            "retry_id": retry_id, "run_id": run_id, "locus_id": locus_id,
            "pair_id": pair_id, "original_reason": original_reason,
            "original_receipt_sha256": original_receipt_sha,
            "retry_worker_sha256": sha256(worker),
            "failed_attempt_directories": [str(path) for path in failures],
        })
        raise RuntimeError(f"Isolated retry failed with exit code {exit_code}; original receipt is unchanged")

    retry_receipt = verify_unit(unit_target, locus_id, set(EXPECTED_PAIRS), run_id)
    if retry_receipt is None:
        raise RuntimeError("Retry worker returned success without a verifiable locus receipt")
    if (sha256(original_receipt_path) != original_receipt_sha or
            verify_unit(original_unit, locus_id, set(EXPECTED_PAIRS), run_id) is None or
            sha256(original_results_path) != original_results_sha):
        raise RuntimeError("Original roundoff receipt changed during the isolated retry")
    retry_rows = [r for r in rows(unit_target / "worker_output/pair_results.tsv") if r["pair_id"] == pair_id]
    if len(retry_rows) != 1 or retry_rows[0]["status"] == "FAILED":
        raise RuntimeError("Retry did not resolve the targeted environmental failure")

    pair_cfg = next(item for item in source_config["pairs"] if item["pair_id"] == pair_id)
    overlap_path = Path(pair_cfg["sample_overlap_file"])
    input_info_path = Path(pair_cfg["input_info"])
    immutable_json(retry_record_path, {
        "status": "PASS_RETRY_ARTIFACT",
        "retry_id": retry_id, "run_id": run_id,
        "run_identity_sha256": run_id,
        "locus_id": locus_id, "pair_id": pair_id,
        "original_status": "FAILED_ENVIRONMENTAL_TEMP_DIRECTORY",
        "original_reason": original_reason,
        "original_receipt_path": str(original_receipt_path),
        "original_receipt_sha256": original_receipt_sha,
        "original_results_sha256": sha256(original_results_path),
        "original_worker_config_sha256": sha256(source_config_path),
        "family_lock_sha256": sha256(family_path),
        "execution_lock_sha256": sha256(execution_path),
        "materialized_inputs_manifest_sha256": sha256(input_root / "provenance.json"),
        "locus_input_hashes": locus_input_hashes,
        "locus_input_info_sha256": sha256(input_info_path),
        "pair_overlap_file_sha256": sha256(overlap_path),
        "reference_provenance_sha256": sha256(reference / "reference.provenance.json"),
        "worker_sha256": sha256(worker),
        "retry_coordinator_sha256": sha256(Path(__file__)),
        "runtime_validator_sha256": sha256(validator),
        "runtime_validation_output_sha256": runtime_sha,
        "rscript_binary_sha256": sha256(rscript),
        "retry_receipt_sha256": sha256(unit_target / "receipt.json"),
        "retry_pair_result": retry_rows[0],
        "temporary_directory": str(temp_dir),
        "receipt_boundary": "Separate targeted retry. Does not replace or merge the original locus receipt or family result.",
    })
    return retry_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--locus-id", default="740")
    parser.add_argument("--pair-id", default="longsleep__parkinson")
    parser.add_argument("--retry-id", required=True, help="Unique immutable output subdirectory name")
    parser.add_argument("--temp-dir", type=Path, default=DEFAULT_TEMP)
    args = parser.parse_args()
    print(retry(args.manifest.resolve(), args.locus_id, args.pair_id,
                args.retry_id, args.temp_dir.resolve()))


if __name__ == "__main__":
    main()
