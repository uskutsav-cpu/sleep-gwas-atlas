#!/usr/bin/env python3
"""Run the frozen Brain6 five-pair × 2,495-locus LAVA family resumably.

Each fresh R process handles one prespecified locus and all five pairs. Work
is checkpointed by locus, all failures remain family slots, and the final
local-rg table is corrected over the fixed family through Brain6's validator.
"""
from __future__ import annotations

import argparse
import atexit
import csv
import hashlib
import io
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "extensions" / "brain6"))
from brain6.artifacts import verify_artifact
from brain6.family24 import local_family
from brain6.io import read_json, sha256, write_json

EXPECTED_PAIRS = ("insomnia__adhd", "insomnia__mdd", "longsleep__scz",
                  "longsleep__bipolar", "longsleep__parkinson")
PAIR_FIELDS = ("pair_id", "locus_id", "status", "p", "local_rg", "reason")
VALID_STATUSES = {"TESTED", "UNIVARIATE_UNDERPOWERED", "NO_OVERLAP", "FAILED"}


def rows(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_bytes_immutable(path: Path, content: bytes) -> None:
    if path.exists():
        if path.read_bytes() != content:
            raise FileExistsError(f"Refusing to replace changed immutable run input: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.rename(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass
        raise


def remove_appledouble_sidecars(directory: Path) -> None:
    """Drop macOS resource-fork companions that external filesystems add."""
    for item in directory.iterdir():
        if item.name.startswith("._") and item.is_file():
            item.unlink()


def require_preflight_success(locus_id: str, worker_exit: int) -> None:
    if worker_exit != 0:
        raise RuntimeError(f"LAVA first-locus preflight failed for locus {locus_id} (worker exit {worker_exit})")


def tsv_bytes(fields: list[str], records: list[dict[str, Any]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    return stream.getvalue().encode("utf-8")


def load_loci(path: Path, expected_sha: str, expected_count: int) -> list[dict[str, str]]:
    if sha256(path) != expected_sha:
        raise ValueError("LAVA locus-definition checksum differs from the frozen family")
    with path.open(encoding="utf-8") as handle:
        lines = [line.split() for line in handle if line.strip()]
    if not lines or lines[0] != ["LOC", "CHR", "START", "STOP"]:
        raise ValueError("LAVA locus-definition header drifted")
    if any(len(row) != 4 for row in lines[1:]):
        raise ValueError("Malformed LAVA locus-definition row")
    values = [dict(zip(lines[0], row)) for row in lines[1:]]
    if len(values) != expected_count or not all({"LOC", "CHR", "START", "STOP"} <= set(r) for r in values):
        raise ValueError("LAVA locus file does not match the exact frozen locus family")
    ids = [str(row["LOC"]) for row in values]
    if len(ids) != len(set(ids)):
        raise ValueError("LAVA locus identifiers are duplicated")
    return values


def validate_inputs(family_path: Path, execution_path: Path, input_root: Path,
                    reference: Path) -> tuple[dict, dict, dict[str, dict], list[dict[str, str]], dict]:
    family = read_json(family_path)
    execution = read_json(execution_path)
    if family.get("analysis_id") != "brain6-lava-local-rg-v2" or family.get("scientific_status") != "NOT_RUN":
        raise ValueError("Active LAVA family lock is not the expected unrun v2 family")
    if sha256(family_path) != "18bc23cfda9d3752a6233ee24c75aa51c4d4dec1e375c469a2a9548c1da60105":
        raise ValueError("Active LAVA family lock checksum mismatch")
    if execution.get("analysis_id") != family["analysis_id"] or execution.get("family_lock_sha256") != sha256(family_path):
        raise ValueError("LAVA execution-policy lock does not bind the active family lock")
    if tuple(family["pairs"]) != EXPECTED_PAIRS:
        raise ValueError("The expected Brain6 pair order has changed")
    pair_lock_path = ROOT / "extensions/brain6/work/overnight-v03/brain6-pairs.lock.json"
    pair_lock = read_json(pair_lock_path)
    lock_body = {key: value for key, value in pair_lock.items() if key != "lock_sha256"}
    if hashlib.sha256(json.dumps(lock_body, sort_keys=True, separators=(",", ":")).encode()).hexdigest() != family["pair_lock_sha256"]:
        raise ValueError("Reviewed pair lock does not match the Brain6 family")
    pair_rows = {row["pair_id"]: row for row in pair_lock["pairs"] if row["pair_id"] in EXPECTED_PAIRS}
    if set(pair_rows) != set(EXPECTED_PAIRS):
        raise ValueError("Reviewed pair lock does not cover the exact LAVA family")

    ref_lock = read_json(reference / "reference.provenance.json")
    if ref_lock.get("verification") != "SHA-256 verified after official HTTPS acquisition":
        raise ValueError("Official LAVA reference is not sealed")
    overlap = family["overlap_gate"]["pairwise_sample_overlap_file_manifest"]
    if overlap.get("participant_level_overlap_verified") is not False or overlap.get("assume_zero_overlap") is not False:
        raise ValueError("Overlap interpretation changed; participant-level overlap must remain explicitly unverified")
    matrix_by_pair = {record["pair_id"]: record for record in overlap["matrices"]}
    if set(matrix_by_pair) != set(EXPECTED_PAIRS):
        raise ValueError("LAVA overlap matrices do not cover all five pairs")
    for pair_id, matrix in matrix_by_pair.items():
        matrix_path = ROOT / matrix["path"]
        if sha256(matrix_path) != matrix["sha256"]:
            raise ValueError(f"Frozen pairwise overlap matrix changed: {pair_id}")
        lock_pair = pair_rows[pair_id]
        if (lock_pair["sleep_trait"], lock_pair["disease_trait"]) != (matrix_path.stem.split("__")[0], matrix_path.stem.split("__")[1]):
            raise ValueError(f"Overlap matrix trait identity does not match pair lock: {pair_id}")

    materialized_path = input_root / "provenance.json"
    materialized = read_json(materialized_path)
    if materialized.get("schema_version") != "brain6-lava-locus-inputs.1" or materialized.get("family_lock_sha256") != sha256(family_path):
        raise ValueError("LAVA chromosome inputs are not bound to the active family")
    if materialized.get("materializer_script_sha256") != sha256(ROOT / "brain6/scripts/materialize_lava_inputs.py"):
        raise ValueError("LAVA input-materializer implementation differs from its provenance lock")
    if materialized.get("case_control_ledger_sha256") != sha256(ROOT / "brain6/manifests/gwas_master.tsv"):
        raise ValueError("Locked case/control sample sizes changed after LAVA input materialization")
    if materialized.get("official_reference_provenance_sha256") != sha256(reference / "reference.provenance.json"):
        raise ValueError("LAVA chromosome inputs reference a different UKB reference lock")
    loci_file = ROOT / family["locus_definition"]["path"]
    loci = load_loci(loci_file, family["locus_definition"]["sha256"], int(family["locus_definition"]["n_loci"]))
    expected_loci = {str(row["LOC"]) for row in loci}
    if int(materialized.get("n_loci", -1)) != len(expected_loci):
        raise ValueError("Materialized LAVA input locus count differs from the frozen family")
    file_records = {}
    for record in materialized.get("records", []):
        if record.get("kind") == "sumstats":
            p = Path(record["path"])
            if p.stat().st_size != int(record["bytes"]) or sha256(p) != record["sha256"]:
                raise ValueError(f"Materialized LAVA input changed: {p}")
            file_records[(str(record["locus_id"]), record["trait"])] = record
        elif record.get("kind") == "input_info":
            p = Path(record["path"])
            if p.stat().st_size != int(record["bytes"]) or sha256(p) != record["sha256"]:
                raise ValueError(f"Materialized LAVA input-info changed: {p}")
            file_records[(str(record["locus_id"]), "input_info")] = record
    if len(file_records) != len(expected_loci) * 8:
        raise ValueError("Materialized LAVA input inventory is not exactly 2,495 loci × (7 traits + input info)")
    expected_inventory = {(locus, trait) for locus in expected_loci for trait in family["trait_ids"]}
    expected_inventory |= {(locus, "input_info") for locus in expected_loci}
    if set(file_records) != expected_inventory:
        raise ValueError("Materialized LAVA input inventory has missing, duplicate, or unexpected trait/locus keys")
    coverage = materialized.get("pair_coverage", {})
    coverage_path = Path(coverage.get("path", ""))
    if coverage_path.stat().st_size != int(coverage.get("bytes", -1)) or sha256(coverage_path) != coverage.get("sha256"):
        raise ValueError("Pair-by-locus reference-coverage table differs from its manifest")
    coverage_map = {}
    for row in rows(coverage_path):
        key = (row["pair_id"], row["locus_id"])
        if key in coverage_map or key[0] not in EXPECTED_PAIRS or key[1] not in expected_loci:
            raise ValueError(f"Duplicate or unexpected pair-locus coverage row: {key}")
        shared = int(row["shared_reference_variants"])
        if shared < 0 or (row["status"] == "READY") != (shared >= 2):
            raise ValueError(f"Invalid pair-locus overlap accounting: {key}")
        coverage_map[key] = shared
    if set(coverage_map) != {(pair, locus) for pair in EXPECTED_PAIRS for locus in expected_loci}:
        raise ValueError("Pair-locus coverage table does not preserve all 12,475 frozen slots")
    return family, execution, pair_rows, loci, coverage_map


def expected_pair_config(pair_id: str, pair_lock: dict, locus_id: str, chrom: int,
                         input_root: Path, family: dict, coverage_map: dict) -> dict[str, Any]:
    row = pair_lock[pair_id]
    matrix = next(r for r in family["overlap_gate"]["pairwise_sample_overlap_file_manifest"]["matrices"]
                  if r["pair_id"] == pair_id)
    return {
        "pair_id": pair_id,
        "phenotypes": [row["sleep_trait"], row["disease_trait"]],
        "input_info": str(input_root / f"locus_{locus_id}" / "input_info.tsv"),
        "sample_overlap_file": str(ROOT / matrix["path"]),
        "reference_prefix": str((ROOT / "ref/lava/ukb_v1.1" / "lava-ukb-v1.1")) + f"_chr{chrom}",
        "shared_reference_variants": int(coverage_map[(pair_id, locus_id)]),
    }


def verify_unit(path: Path, locus_id: str, expected_pairs: set[str], run_id: str) -> dict[str, Any] | None:
    receipt_path = path / "receipt.json"
    if not receipt_path.is_file():
        return None
    receipt = read_json(receipt_path)
    if receipt.get("locus_id") != locus_id or receipt.get("run_id") != run_id:
        raise ValueError(f"Locus receipt identity mismatch: {locus_id}")
    if receipt.get("status") != "COMPLETE":
        raise ValueError(f"Locus receipt is not complete: {locus_id}")
    config_path = path / "worker_config.json"
    if not config_path.is_file() or sha256(config_path) != receipt.get("worker_config_sha256"):
        raise ValueError(f"Locus worker configuration changed: {locus_id}")
    outputs = receipt.get("outputs", [])
    expected_outputs = {"worker_output/pair_results.tsv", "worker_output/univariate.tsv",
                        "worker_output/status.json"}
    if {output.get("path") for output in outputs} != expected_outputs:
        raise ValueError(f"Locus checkpoint output inventory changed: {locus_id}")
    for output in outputs:
        item = path / output["path"]
        if not item.is_file() or item.stat().st_size != output["bytes"] or sha256(item) != output["sha256"]:
            raise ValueError(f"Corrupt LAVA locus checkpoint: {item}")
    result_rows = list(rows(path / "worker_output/pair_results.tsv"))
    if len(result_rows) != len(expected_pairs) or {r["pair_id"] for r in result_rows} != expected_pairs:
        raise ValueError(f"Locus checkpoint does not cover all selected pairs: {locus_id}")
    if any(r["status"] not in VALID_STATUSES for r in result_rows):
        raise ValueError(f"Unknown locus result status: {locus_id}")
    return receipt


def run_locus_unit(target: Path, locus: dict[str, str], family: dict, execution: dict,
                   pair_lock: dict, input_root: Path, family_path: Path,
                   execution_path: Path, run_id: str, run_identity: dict,
                   worker: Path, validator: Path, rscript: Path,
                   execution_policy: dict, coverage_map: dict) -> int:
    locus_id = str(locus["LOC"])
    chromosome = int(locus["CHR"])
    maximum_attempts = int(execution.get("maximum_process_attempts_per_locus", 1))
    if maximum_attempts < 1:
        raise ValueError("maximum_process_attempts_per_locus must be positive")
    failed_root = target.parent.parent / "failed_attempts"
    for attempt_number in range(1, maximum_attempts + 1):
        attempt = Path(tempfile.mkdtemp(prefix=f".{locus_id}.attempt-{attempt_number}-", dir=target.parent))
        config_path = attempt / "worker_config.json"
        out_dir = attempt / "worker_output"
        out_dir.mkdir()
        cfg = {
            "schema_version": 1,
            "analysis_id": family["analysis_id"],
            "family_lock_sha256": sha256(family_path),
            "execution_lock_sha256": sha256(execution_path),
            "locus_id": locus_id,
            "loci_file": str((ROOT / family["locus_definition"]["path"]).resolve()),
            "random_seed": int(family["execution"]["random_seed"]),
            "execution_policy": execution_policy,
            "rscript": str(rscript.resolve()),
            "runtime_validator": str(validator.resolve()),
            "runtime_validator_passed": True,
            "runtime_validation_output_sha256": run_identity["runtime_validation_output_sha256"],
            "output_dir": str(out_dir.resolve()),
            "pairs": [expected_pair_config(pair, pair_lock, locus_id, chromosome, input_root,
                                             family, coverage_map)
                      for pair in family["pairs"]],
        }
        write_json(config_path, cfg)
        proc = subprocess.run([str(rscript), str(worker), str(config_path)], cwd=ROOT,
                              capture_output=True, text=True, check=False)
        (attempt / "worker.log").write_text(proc.stdout + proc.stderr, encoding="utf-8")
        remove_appledouble_sidecars(out_dir)
        if proc.returncode == 0:
            pair_path = out_dir / "pair_results.tsv"
            result_rows = list(rows(pair_path)) if pair_path.is_file() else []
            if len(result_rows) != 5 or {r.get("pair_id") for r in result_rows} != set(family["pairs"]):
                proc = subprocess.CompletedProcess(proc.args, 1, proc.stdout,
                    "Worker output failed exact five-pair validation")
            elif any(r.get("status") not in VALID_STATUSES for r in result_rows):
                proc = subprocess.CompletedProcess(proc.args, 1, proc.stdout,
                    "Worker output contains an unknown result status")
        if proc.returncode == 0:
            outputs = []
            for item in sorted(out_dir.iterdir()):
                if item.is_file():
                    outputs.append({"path": str(item.relative_to(attempt)), "bytes": item.stat().st_size,
                                    "sha256": sha256(item)})
            expected_outputs = {"worker_output/pair_results.tsv", "worker_output/univariate.tsv",
                                "worker_output/status.json"}
            if {record["path"] for record in outputs} != expected_outputs:
                proc = subprocess.CompletedProcess(proc.args, 1, proc.stdout,
                    "Worker output inventory is incomplete or unexpected")
            else:
                receipt = {"schema_version": 1, "locus_id": locus_id,
                    "run_id": run_id, "run_identity": run_identity,
                    "worker_config_sha256": sha256(config_path), "outputs": outputs,
                    "status": "COMPLETE"}
                write_json(attempt / "receipt.json", receipt)
                os.rename(attempt, target)
                return 0
        failure = attempt / "failure.json"
        write_json(failure, {"locus_id": locus_id, "attempt": attempt_number,
            "returncode": proc.returncode,
            "stdout_tail": proc.stdout[-3000:], "stderr_tail": proc.stderr[-3000:],
            "worker_log_sha256": sha256(attempt / "worker.log")})
        failed_root.mkdir(exist_ok=True)
        os.rename(attempt, failed_root / f"{locus_id}-{attempt_number:02d}-{len(list(failed_root.iterdir())) + 1:04d}")
    return proc.returncode


def run(family_path: Path, execution_path: Path, input_root: Path, reference: Path,
        output_root: Path, rscript: Path, *, only_locus: str | None = None) -> Path:
    family, execution, pair_lock, loci, coverage_map = validate_inputs(
        family_path, execution_path, input_root, reference)
    if not rscript.is_file():
        raise FileNotFoundError(f"Pinned Rscript is missing: {rscript}")
    validator = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
    runtime = subprocess.run([str(rscript), str(validator)], cwd=ROOT,
                             capture_output=True, text=True, check=False)
    if runtime.returncode or runtime.stdout.count("TRACK_B_LAVA_RUNTIME\tstatus=PASS\t") != 1:
        raise RuntimeError("Pinned LAVA runtime failed preflight: " + (runtime.stdout + runtime.stderr)[-4000:])
    reference_check = subprocess.run([sys.executable, "scripts/lava_contract.py", "--verify-reference"],
        cwd=ROOT, capture_output=True, text=True, check=False)
    if reference_check.returncode or "LAVA_REFERENCE_VALIDATED" not in reference_check.stdout:
        raise RuntimeError("Official LAVA reference failed preflight: " +
                           (reference_check.stdout + reference_check.stderr)[-4000:])
    worker = ROOT / "brain6/scripts/run_lava_family_locus.R"
    materializer = ROOT / "brain6/scripts/materialize_lava_inputs.py"
    family_aggregator = ROOT / "extensions/brain6/brain6/family24.py"
    run_identity = {
        "analysis_id": family["analysis_id"],
        "family_lock_sha256": sha256(family_path),
        "execution_lock_sha256": sha256(execution_path),
        "materialized_inputs_sha256": sha256(input_root / "provenance.json"),
        "materializer_script_sha256": sha256(materializer),
        "orchestrator_script_sha256": sha256(Path(__file__)),
        "family_aggregator_source_sha256": sha256(family_aggregator),
        "worker_sha256": sha256(worker),
        "runtime_validator_sha256": sha256(validator),
        "runtime_validation_output_sha256": hashlib.sha256(runtime.stdout.encode("utf-8")).hexdigest(),
        "rscript": str(rscript.resolve()),
        "rscript_binary_sha256": sha256(rscript.resolve()),
        "r_version": "4.3.3",
        "lava_version": "0.1.5",
        "lava_commit": "e729a245f7b6923967a96804fbf5246eadf2d6c6",
    }
    run_id = hashlib.sha256(json.dumps(run_identity, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    run_root = output_root / run_id
    output_root.mkdir(parents=True, exist_ok=True)
    lock_path = output_root / f".{run_id}.lock"
    try:
        lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise RuntimeError(f"LAVA family run is locked; inspect the lock owner before resuming: {lock_path}") from error
    with os.fdopen(lock_fd, "w", encoding="utf-8") as handle:
        json.dump({"pid": os.getpid(), "run_id": run_id}, handle, sort_keys=True)
        handle.write("\n")
    def release_run_lock():
        lock_path.unlink(missing_ok=True)
    atexit.register(release_run_lock)
    units_root = run_root / "loci"
    units_root.mkdir(parents=True, exist_ok=True)
    identity_path = run_root / "run_identity.json"
    if identity_path.exists():
        if read_json(identity_path) != run_identity:
            raise ValueError("Existing LAVA run identity differs from its inputs")
    else:
        write_json(identity_path, run_identity)
    execution_policy = {
        "locus_processing": execution["locus_processing"],
        "univariate": execution["univariate"],
        "bivariate": execution["bivariate"],
    }
    expected_pairs = set(family["pairs"])
    selected_loci = loci if only_locus is None else [r for r in loci if str(r["LOC"]) == only_locus]
    if only_locus is not None and len(selected_loci) != 1:
        raise ValueError(f"Requested LAVA locus is absent/duplicated: {only_locus}")
    if only_locus is not None and only_locus != str(loci[0]["LOC"]):
        raise ValueError("The optional preflight may run only the first prespecified locus")

    for index, locus in enumerate(selected_loci, 1):
        locus_id = str(locus["LOC"])
        target = units_root / locus_id
        if verify_unit(target, locus_id, expected_pairs, run_id) is not None:
            continue
        worker_exit = run_locus_unit(target, locus, family, execution, pair_lock, input_root,
            family_path, execution_path, run_id, run_identity, worker, validator, rscript,
            execution_policy, coverage_map)
        if only_locus is not None:
            print(f"BRAIN6_LAVA_PREFLIGHT_LOCUS locus={locus_id} worker_exit={worker_exit}")
            require_preflight_success(locus_id, worker_exit)
            return run_root
        elif index % 25 == 0 or index == len(selected_loci):
            complete = sum((units_root / str(row["LOC"]) / "receipt.json").is_file() for row in selected_loci[:index])
            print(f"BRAIN6_LAVA_PROGRESS completed_loci={complete}/{len(selected_loci)}")

    planned_path = run_root / "planned_units.tsv"
    planned = [dict(pair_id=pair, locus_id=str(locus["LOC"]))
               for locus in loci for pair in family["pairs"]]
    write_bytes_immutable(planned_path, tsv_bytes(["pair_id", "locus_id"], planned))

    result_rows = []
    receipt_rows = []
    unique_univariate = {}
    for locus in loci:
        locus_id = str(locus["LOC"])
        checkpoint = units_root / locus_id
        receipt = verify_unit(checkpoint, locus_id, expected_pairs, run_id)
        by_pair = {row["pair_id"]: row for row in rows(checkpoint / "worker_output/pair_results.tsv")} if receipt else {}
        if receipt:
            receipt_rows.append({"locus_id": locus_id, "receipt_path": str(checkpoint / "receipt.json"),
                                 "receipt_sha256": sha256(checkpoint / "receipt.json")})
            for univariate in rows(checkpoint / "worker_output/univariate.tsv"):
                trait = univariate.get("phen")
                if trait not in family["trait_ids"]:
                    raise ValueError(f"Unexpected LAVA univariate trait at {locus_id}: {trait}")
                p_value = float(univariate["p"])
                if not math.isfinite(p_value) or not 0 <= p_value <= 1:
                    raise ValueError(f"Invalid LAVA univariate P at {trait}/{locus_id}")
                key = (trait, locus_id)
                if key in unique_univariate and not math.isclose(unique_univariate[key], p_value, rel_tol=0, abs_tol=1e-12):
                    raise ValueError(f"Pair-specific inputs disagree on the shared univariate test: {trait}/{locus_id}")
                unique_univariate[key] = p_value
        for pair in family["pairs"]:
            row = by_pair.get(pair)
            if row is None:
                row = {"pair_id": pair, "locus_id": locus_id, "status": "FAILED", "p": "NA",
                       "local_rg": "NA", "reason": "LOCUS_PROCESS_NOT_COMPLETED"}
            else:
                row["locus_id"] = locus_id
            result_rows.append({field: row.get(field, "") for field in PAIR_FIELDS})
    if len(result_rows) != int(family["bivariate_family"]["n_slots"]):
        raise ValueError("Result aggregation did not preserve the exact frozen pair-locus denominator")
    results_path = run_root / "results.tsv"
    write_bytes_immutable(results_path, tsv_bytes(list(PAIR_FIELDS), result_rows))
    univariate_records = []
    expected_univariate = int(family["univariate_gate"]["n_tests"])
    for locus in loci:
        locus_id = str(locus["LOC"])
        for trait in family["trait_ids"]:
            p_value = unique_univariate.get((trait, locus_id))
            univariate_records.append({
                "trait_id": trait, "locus_id": locus_id,
                "status": "TESTED" if p_value is not None else "NOT_RUN",
                "p": p_value if p_value is not None else "NA",
                "reason": "" if p_value is not None else "NO_VALID_PAIRWISE_LOCAL_UNIVARIATE_RESULT",
            })
    if len(univariate_records) != expected_univariate:
        raise ValueError("Unique trait-by-locus univariate family denominator drifted")
    univariate_path = run_root / "univariate_results.tsv"
    write_bytes_immutable(univariate_path, tsv_bytes(
        ["trait_id", "locus_id", "status", "p", "reason"], univariate_records))
    receipt_index = run_root / "worker_receipts.tsv"
    write_bytes_immutable(receipt_index, tsv_bytes(["locus_id", "receipt_path", "receipt_sha256"], receipt_rows))
    manifest_path = run_root / "local_family_manifest.json"
    analysis_manifest = {
        "reviewed": True,
        "synthetic": False,
        "analysis_id": family["analysis_id"],
        "family_lock_sha256": sha256(family_path),
        "execution_lock_sha256": sha256(execution_path),
        "run_identity_sha256": sha256(run_root / "run_identity.json"),
        "worker_receipts_path": str(receipt_index),
        "worker_receipts_sha256": sha256(receipt_index),
        "univariate_results_path": str(univariate_path),
        "univariate_results_sha256": sha256(univariate_path),
        "planned_units": str(planned_path),
        "results": str(results_path),
        "maximum_failure_rate": float(family["execution"]["maximum_locus_failure_fraction"]),
    }
    if not manifest_path.exists():
        write_json(manifest_path, analysis_manifest)
    elif read_json(manifest_path) != analysis_manifest:
        raise ValueError("Existing local-family manifest differs from the frozen run")
    final_name = f"brain6_local_rg_family_{run_id}"
    final_path = output_root / final_name
    if final_path.exists():
        verify_artifact(final_path)
        final = final_path
    else:
        final = local_family(manifest_path, output_root, final_name)
    status = read_json(final / "status.json")
    unique_univariate_missing = sum(row["status"] == "NOT_RUN" for row in univariate_records)
    missing_rate = unique_univariate_missing / expected_univariate
    univariate_cap = float(family["execution"]["maximum_univariate_untested_fraction"])
    overall_status = status["status"]
    if missing_rate > univariate_cap:
        overall_status = "FAILED_QC_NOT_CONSUMED"
    local_rows = list(rows(final / "local_family.tsv"))
    correction_rows = []
    for row in local_rows:
        tested = row["status"] == "TESTED"
        correction_rows.append({
            "pair_id": row["pair_id"], "locus_id": row["locus_id"],
            "status": row["status"], "observed_p": row["p"],
            "local_rg": row["local_rg"],
            "p_for_family_correction": row["p"] if tested else 1.0,
            "family_fdr": row["family_fdr"] if tested else 1.0,
        })
    correction_path = run_root / "family_correction.tsv"
    write_bytes_immutable(correction_path, tsv_bytes(
        ["pair_id", "locus_id", "status", "observed_p", "local_rg",
         "p_for_family_correction", "family_fdr"], correction_rows))
    decision = {
        "analysis_id": family["analysis_id"],
        "run_identity_sha256": sha256(run_root / "run_identity.json"),
        "local_family_artifact": str(final),
        "local_family_receipt_sha256": sha256(final / "receipt.json"),
        "local_family_status": status["status"],
        "overall_status": overall_status,
        "planned_bivariate_slots": len(correction_rows),
        "non_tested_slots_use_p_one": True,
        "family_correction_path": str(correction_path),
        "family_correction_sha256": sha256(correction_path),
        "unique_univariate_tests": expected_univariate,
        "unique_univariate_not_run": unique_univariate_missing,
        "unique_univariate_not_run_fraction": missing_rate,
        "maximum_univariate_untested_fraction": univariate_cap,
        "unique_univariate_gate_misses": sum(
            row["status"] == "TESTED" and float(row["p"]) >= float(family["univariate_gate"]["p_threshold_strictly_less_than"])
            for row in univariate_records),
        "participant_level_overlap_verified": False,
    }
    decision_path = run_root / "family_decision.json"
    decision_bytes = (json.dumps(decision, indent=2, sort_keys=True) + "\n").encode()
    write_bytes_immutable(decision_path, decision_bytes)
    print(f"BRAIN6_LAVA_FAMILY status={overall_status} planned={status['planned_units']} "
          f"execution_failures_or_missing={status['execution_failures_or_missing']} "
          f"unique_univariate_not_run={unique_univariate_missing}")
    release_run_lock()
    atexit.unregister(release_run_lock)
    return final


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family-lock", type=Path, default=ROOT / "brain6/config/lava_family_v2.json")
    parser.add_argument("--execution-lock", type=Path, default=ROOT / "brain6/config/lava_execution_v1.json")
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=ROOT / "ref/lava/ukb_v1.1")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--rscript", type=Path, default=ROOT / ".r-env/bin/Rscript")
    parser.add_argument("--only-locus", help="Run one first-class frozen slot for runtime preflight; omit for the full family")
    args = parser.parse_args()
    result = run(args.family_lock.resolve(), args.execution_lock.resolve(), args.input_root.resolve(),
                 args.reference.resolve(), args.output_root.resolve(), args.rscript.resolve(),
                 only_locus=args.only_locus)
    print(f"BRAIN6_LAVA_OUTPUT path={result}")


if __name__ == "__main__":
    main()
