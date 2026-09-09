#!/usr/bin/env python3
"""Run one Track B PLACO+ full-P pair through the additive terminal gate.

The frozen V1 materialization/validation engine is reused without modification.
An additive exact-science copy of the V1 R runner binds worker identity into
resumes.  A V2 task binds the additive terminal-gate lock; the bridge never passes
LAVA result rows to PLACO+, and it never performs LD clumping or locus publication.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager, nullcontext
import fcntl
import hashlib
import importlib.util
import json
import math
import os
import shutil
import stat
import subprocess
from pathlib import Path
from typing import Any, Callable, Iterator


ROOT = Path(__file__).resolve().parents[1]
ENGINE_PATH = ROOT / "scripts/125_materialize_track_b_placo_pair.py"
GATE_PATH = ROOT / "scripts/139_build_track_b_placo_terminal_gate_v2.py"
V1_RUNNER = Path("scripts/126_run_track_b_placo_pair.R")
V2_RUNNER = Path("scripts/141_run_track_b_placo_pair_v2.R")
NUISANCE_LOCK_FIELDS = [
    "run_fingerprint", "input_sha256", "task_sha256", "placo_source_sha256",
    "runner_sha256", "workers", "marginal_p_threshold", "variance_null_rows",
    "correlation_null_rows", "nuisance_rds_sha256",
]


class BridgeError(RuntimeError):
    """A fail-closed V2 bridge violation."""


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise BridgeError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


ENGINE = load_module("track_b_placo_v1_engine_for_v2_bridge", ENGINE_PATH)
GATE = load_module("track_b_placo_terminal_gate_for_v2_bridge", GATE_PATH)
ENGINE.RUNNER = V2_RUNNER

# This is deliberately one repository-wide lock for A, B, and CONTROL.  The
# file persists, but flock ownership never does, so process death is stale-safe.
EXECUTION_LOCK = ROOT / GATE.GATE_ROOT / ".execution.lock"


@contextmanager
def execution_lock(action: str, pair_id: str) -> Iterator[None]:
    path = EXECUTION_LOCK
    descriptor: int | None = None
    try:
        GATE.require_real_repository_path(path.parent, "PLACO V2 execution-lock parent")
        path.parent.mkdir(parents=True, exist_ok=True)
        GATE.require_real_repository_path(path.parent, "PLACO V2 execution-lock parent")
        flags = os.O_RDWR | os.O_CREAT
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path, flags, 0o600)
        opened = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            not stat.S_ISREG(opened.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (opened.st_dev, opened.st_ino) != (current.st_dev, current.st_ino)
            or opened.st_nlink != 1
            or current.st_nlink != 1
        ):
            raise BridgeError("PLACO V2 execution lock is not a private stable regular file")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise BridgeError(
                "BLOCKED_BY_COMPUTE another A/B/CONTROL PLACO V2 action holds the repository lock"
            ) from error
        locked = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            not stat.S_ISREG(locked.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (locked.st_dev, locked.st_ino) != (current.st_dev, current.st_ino)
            or locked.st_nlink != 1
            or current.st_nlink != 1
        ):
            raise BridgeError("PLACO V2 execution-lock path changed while acquiring flock")
        record = f"pid={os.getpid()} pair={pair_id} action={action}\n".encode()
        os.ftruncate(descriptor, 0)
        if os.write(descriptor, record) != len(record):
            raise BridgeError("could not write complete PLACO V2 execution-lock ownership record")
        os.fsync(descriptor)
        locked = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            (locked.st_dev, locked.st_ino) != (current.st_dev, current.st_ino)
            or locked.st_nlink != 1
            or current.st_nlink != 1
        ):
            raise BridgeError("PLACO V2 execution-lock path changed while recording ownership")
    except BridgeError:
        if descriptor is not None:
            os.close(descriptor)
        raise
    except (OSError, GATE.GateError) as error:
        if descriptor is not None:
            os.close(descriptor)
        raise BridgeError(f"could not acquire the repository-scoped PLACO V2 execution lock: {error}") from error
    try:
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def validate_v2_paths(context: dict[str, Any]) -> None:
    """Reject symlinks/special files anywhere in the mutable V2 namespace."""
    try:
        paths = ENGINE.materialized_paths(context)
        for name, path in paths.items():
            GATE.require_real_repository_path(path, f"PLACO V2 {name} path")
        base = paths["base"]
        if not base.exists():
            return
        for directory, children, files in os.walk(base, followlinks=False):
            directory_path = Path(directory)
            GATE.require_real_repository_path(directory_path, "PLACO V2 work tree")
            for name in [*children, *files]:
                path = directory_path / name
                observed = os.lstat(path)
                if stat.S_ISLNK(observed.st_mode):
                    raise BridgeError(f"PLACO V2 work tree contains a symbolic link: {path}")
                if name in children and not stat.S_ISDIR(observed.st_mode):
                    raise BridgeError(f"PLACO V2 work tree contains a non-directory ancestor: {path}")
                if name in files and not stat.S_ISREG(observed.st_mode):
                    raise BridgeError(f"PLACO V2 work tree contains a non-regular artifact: {path}")
    except (OSError, GATE.GateError) as error:
        raise BridgeError(f"unsafe PLACO V2 execution namespace: {error}") from error


def rollback_v2_publication(published: list[tuple[Path, int, int]]) -> None:
    """Remove only the exact inodes linked by this V2 publication attempt."""
    failure: BaseException | None = None
    for path, device, inode in reversed(published):
        try:
            GATE.unlink_if_identity(path, (device, inode))
        except BaseException as error:  # preserve every other rollback attempt
            if failure is None:
                failure = error
    published.clear()
    if failure is not None:
        raise BridgeError(f"could not safely roll back PLACO V2 publication: {failure}") from failure


def exclusive_family_v2(
    sources: list[Path],
    destinations: list[Path],
    *,
    before_publication: Callable[[], None],
    after_publication: Callable[[], None],
    committed: list[tuple[Path, int, int]],
) -> None:
    """Hard-link one no-replace family with V2 path and inode protections."""
    if len(sources) != len(destinations):
        raise BridgeError("internal PLACO V2 publication family length mismatch")
    source_identities: list[tuple[int, int, int, int]] = []
    for source in sources:
        try:
            GATE.require_real_repository_path(source, "PLACO V2 publication source")
            observed = os.lstat(source)
        except (OSError, GATE.GateError) as error:
            raise BridgeError(f"unsafe PLACO V2 publication source: {error}") from error
        if not stat.S_ISREG(observed.st_mode) or observed.st_size <= 0:
            raise BridgeError(f"PLACO V2 publication source is not a real non-empty file: {source}")
        # ctime necessarily changes when the same source is hard-linked to
        # multiple destinations (the provenance source is deliberately used
        # twice), so stability here binds inode, size, and content mtime.
        source_identities.append(
            (observed.st_dev, observed.st_ino, observed.st_size, observed.st_mtime_ns)
        )
    for destination in destinations:
        try:
            GATE.require_real_repository_path(destination.parent, "PLACO V2 publication parent")
            destination.parent.mkdir(parents=True, exist_ok=True)
            GATE.require_real_repository_path(destination.parent, "PLACO V2 publication parent")
            GATE.require_real_repository_path(destination, "PLACO V2 publication destination")
        except (OSError, GATE.GateError) as error:
            raise BridgeError(f"unsafe PLACO V2 publication destination: {error}") from error
        if destination.exists() or destination.is_symlink():
            raise BridgeError(f"immutable PLACO V2 publication destination already exists: {destination}")

    local: list[tuple[Path, int, int]] = []
    try:
        # This callback deliberately runs inside the V1 engine at the last
        # possible point, after its complete-ledger semantic validation and
        # before the first canonical hard link.
        before_publication()
        for source, destination, expected in zip(
            sources, destinations, source_identities, strict=True,
        ):
            GATE.require_real_repository_path(destination.parent, "PLACO V2 publication parent")
            if destination.exists() or destination.is_symlink():
                raise BridgeError(f"PLACO V2 publication destination appeared concurrently: {destination}")
            source_stat = os.lstat(source)
            observed_source = (
                source_stat.st_dev, source_stat.st_ino, source_stat.st_size, source_stat.st_mtime_ns,
            )
            if observed_source != expected or not stat.S_ISREG(source_stat.st_mode):
                raise BridgeError(f"PLACO V2 publication source changed before linking: {source}")
            try:
                os.link(source, destination, follow_symlinks=False)
            except FileExistsError as error:
                raise BridgeError(
                    f"PLACO V2 publication destination appeared concurrently: {destination}"
                ) from error
            destination_stat = os.lstat(destination)
            local.append((destination, source_stat.st_dev, source_stat.st_ino))
            if (
                not stat.S_ISREG(destination_stat.st_mode)
                or (destination_stat.st_dev, destination_stat.st_ino)
                != (source_stat.st_dev, source_stat.st_ino)
            ):
                raise BridgeError(f"PLACO V2 no-replace publication identity mismatch: {destination}")
            GATE.fsync_directory(destination.parent)
        after_publication()
    except BaseException:
        rollback_v2_publication(local)
        raise
    committed.extend(local)


def inherited_v1_upstream(row: dict[str, str]) -> None:
    """Validate the immutable V1 snapshot while replacing only its local proxy."""
    blockers = []
    if row.get("dense_input_gate") != "PASS_ALL_THREE_PAIRS":
        blockers.append("DENSE_INPUT")
    if not row.get("replication_gate", "").startswith("PASS_TERMINAL"):
        blockers.append("PRIMARY_REPLICATION")
    if row.get("local_analysis_gate") != "BLOCKED_LOCAL_ANALYSIS_NOT_TERMINAL":
        blockers.append("UNEXPECTED_V1_LOCAL_GATE_STATE")
    if row.get("software_gate") != "READY_PINNED_PLACO_SOURCE":
        blockers.append("PLACO_SOURCE")
    if row.get("claim_status") != "NO_SCIENTIFIC_RESULT":
        blockers.append("UNEXPECTED_V1_CLAIM_STATE")
    if blockers:
        raise SystemExit(
            "ERROR: frozen V1 PLACO+ snapshot differs from the additive bridge expectation: "
            + ",".join(blockers)
        )


def configure_context(pair_id: str) -> tuple[dict[str, Any], dict[str, Any]]:
    try:
        gate_lock = GATE.verify_gate()
    except GATE.GateError as error:
        raise BridgeError(str(error)) from error
    evidence = gate_lock.get("lava_scientific_evidence", {})
    if (
        evidence.get("status") != "FAILED_QC_NOT_CONSUMED"
        or evidence.get("scientific_validation_passed") is not False
        or evidence.get("passed_lava_results_used_by_placo") is not False
        or evidence.get("consumption") != "ORDERING_TERMINALITY_ATTESTATION_ONLY"
    ):
        raise BridgeError("V2 gate does not preserve FAILED_QC_NOT_CONSUMED semantics")
    if gate_lock.get("pair_scope") != ["A", "B", "CONTROL"] or pair_id not in ENGINE.PAIR_IDS:
        raise BridgeError("V2 bridge scope is not exactly A/B/CONTROL")

    original_assert = ENGINE.assert_upstream_ready
    original_input_lock = ENGINE.INPUT_LOCK
    original_runner = ENGINE.RUNNER
    ENGINE.assert_upstream_ready = inherited_v1_upstream
    ENGINE.INPUT_LOCK = GATE.V1_INPUT_LOCK
    ENGINE.RUNNER = V2_RUNNER
    try:
        base = ENGINE.production_context(ROOT, pair_id)
    except SystemExit as error:
        ENGINE.RUNNER = original_runner
        raise BridgeError(str(error)) from error
    except BaseException:
        ENGINE.RUNNER = original_runner
        raise
    finally:
        ENGINE.assert_upstream_ready = original_assert
        ENGINE.INPUT_LOCK = original_input_lock

    # Every materialized task and result provenance now binds the additive lock.
    ENGINE.INPUT_LOCK = GATE.TERMINAL_GATE_LOCK
    gate_hash = ENGINE.sha256(ROOT / GATE.TERMINAL_GATE_LOCK)
    payload = {
        "schema_version": "sleep-atlas-track-b-placo-execution-fingerprint.2",
        "analysis_id": base["policy"]["analysis_id"],
        "pair_id": pair_id,
        "traits": [base["trait1"], base["trait2"]],
        "family_role": base["family_role"],
        "policy_sha256": ENGINE.sha256(ROOT / ENGINE.POLICY),
        "v1_contract_lock_sha256": ENGINE.sha256(ROOT / ENGINE.CONTRACT_LOCK),
        "v1_input_gate_lock_sha256": ENGINE.sha256(ROOT / GATE.V1_INPUT_LOCK),
        "v2_terminal_gate_lock_sha256": gate_hash,
        "v1_materialization_validation_engine_sha256": ENGINE.sha256(ENGINE_PATH),
        "v2_materializer_bridge_sha256": ENGINE.sha256(Path(__file__).resolve()),
        "v1_runner_sha256": ENGINE.sha256(ROOT / V1_RUNNER),
        "v2_worker_bound_runner_sha256": ENGINE.sha256(ROOT / V2_RUNNER),
        "placo_source_sha256": base["policy"]["placo_plus"]["source_sha256"],
        "lava_scientific_evidence": "FAILED_QC_NOT_CONSUMED",
        "sources": [
            {key: value for key, value in source.items() if key != "absolute_path"}
            for source in base["sources"]
        ],
        "resource_envelope": gate_lock["resource_envelope"],
        "permitted_scope": "FULL_P_MATERIALIZATION_BENCHMARK_AND_SCAN_ONLY",
        "ld_locus_publication_allowed": False,
    }
    fingerprint = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    context = {
        **base,
        "input_lock": gate_lock,
        "fingerprint": fingerprint,
        "fingerprint_payload": payload,
    }
    validate_v2_paths(context)
    return context, gate_lock


def validate_benchmark(
    context: dict[str, Any], gate_lock: dict[str, Any], required_workers: int | None = None,
) -> dict[str, str]:
    paths = ENGINE.materialized_paths(context)
    try:
        ENGINE.verify_materialized(context)
        fields, rows = ENGINE.read_tsv(paths["benchmark"])
        task = ENGINE.one_tsv_row(paths["task"])
        nuisance_checkpoint = paths["checkpoint_dir"] / "nuisance.rds"
        nuisance_lock_path = paths["checkpoint_dir"] / "nuisance.sha256.tsv"
        nuisance_hash = ENGINE.sha256(nuisance_checkpoint)
        nuisance_lock_fields, nuisance_lock_rows = ENGINE.read_tsv(nuisance_lock_path)
    except (SystemExit, OSError) as error:
        raise BridgeError(f"missing or invalid exact-pair benchmark: {error}") from error
    expected_fields = ENGINE.BENCHMARK_RAW_FIELDS + [
        "peak_process_tree_rss_bytes", "rss_measurement_method", "wrapper_wall_seconds",
        "rss_sample_interval_seconds", "host_physical_memory_bytes", "host_free_bytes_at_benchmark",
        "peak_rss_fraction_of_physical_memory",
    ]
    if fields != expected_fields or len(rows) != 1:
        raise BridgeError("V2 PLACO+ benchmark schema or row count drifted")
    if nuisance_lock_fields != NUISANCE_LOCK_FIELDS or len(nuisance_lock_rows) != 1:
        raise BridgeError("global nuisance checkpoint lock schema or row count drifted")
    row = rows[0]
    nuisance_lock = nuisance_lock_rows[0]
    expected = {
        "analysis_id": task["analysis_id"],
        "pair_id": task["pair_id"],
        "run_fingerprint": task["run_fingerprint"],
        "input_rows": task["aligned_input_rows"],
        "input_sha256": task["aligned_input_sha256"],
        "task_sha256": ENGINE.sha256(paths["task"]),
        "policy_sha256": task["policy_sha256"],
        "contract_lock_sha256": task["contract_lock_sha256"],
        "input_gate_lock_sha256": ENGINE.sha256(ROOT / GATE.TERMINAL_GATE_LOCK),
        "materializer_sha256": ENGINE.sha256(ENGINE_PATH),
        "runner_sha256": ENGINE.sha256(ROOT / ENGINE.RUNNER),
        "placo_source_sha256": task["placo_source_sha256"],
        "nuisance_checkpoint_sha256": nuisance_hash,
        "global_nuisance_input_rows": task["aligned_input_rows"],
        "scientific_equivalence": (
            "GLOBAL_OFFICIAL_NUISANCE_ESTIMATED_ON_ALL_VALID_VARIANTS_"
            "THEN_IDENTICAL_SINGLE_VARIANT_PLACO_PLUS_CALLS_SAMPLED"
        ),
    }
    if any(row.get(key) != value for key, value in expected.items()):
        raise BridgeError("benchmark is not bound to the exact V2 task, software, and global nuisance fit")
    try:
        input_rows = int(row["input_rows"])
        variance_rows = int(row["nuisance_null_rows_variance"])
        correlation_rows = int(row["nuisance_null_rows_correlation"])
        failures = int(row["numerical_failures"])
        benchmark_n = int(row["benchmark_variants"])
        workers = int(row["workers"])
        peak = int(row["peak_process_tree_rss_bytes"])
        host_memory = int(row["host_physical_memory_bytes"])
        host_free = int(row["host_free_bytes_at_benchmark"])
        nuisance_elapsed = float(row["global_nuisance_elapsed_seconds"])
        testing_elapsed = float(row["single_variant_testing_elapsed_seconds"])
        runner_elapsed = float(row["measured_runner_elapsed_seconds"])
        projected = float(row["projected_full_family_seconds"])
        projected_hours = float(row["projected_full_family_hours"])
        rate = float(row["variants_per_second"])
        variance1 = float(row["VarZ1"])
        variance2 = float(row["VarZ2"])
        correlation = float(row["CorZ"])
        wrapper_wall = float(row["wrapper_wall_seconds"])
        sample_interval = float(row["rss_sample_interval_seconds"])
        peak_fraction = float(row["peak_rss_fraction_of_physical_memory"])
        nuisance_threshold = float(nuisance_lock["marginal_p_threshold"])
        locked_workers = int(nuisance_lock["workers"])
        locked_variance_rows = int(nuisance_lock["variance_null_rows"])
        locked_correlation_rows = int(nuisance_lock["correlation_null_rows"])
    except (TypeError, ValueError, OverflowError) as error:
        raise BridgeError("benchmark contains invalid execution measurements") from error
    nuisance_expected = {
        "run_fingerprint": task["run_fingerprint"],
        "input_sha256": task["aligned_input_sha256"],
        "task_sha256": ENGINE.sha256(paths["task"]),
        "placo_source_sha256": task["placo_source_sha256"],
        "runner_sha256": ENGINE.sha256(ROOT / ENGINE.RUNNER),
        "nuisance_rds_sha256": nuisance_hash,
    }
    if (
        any(nuisance_lock.get(key) != value for key, value in nuisance_expected.items())
        or not math.isclose(
            nuisance_threshold, float(task["marginal_p_threshold"]), rel_tol=0, abs_tol=0,
        )
        or locked_workers != workers
        or locked_variance_rows != variance_rows
        or locked_correlation_rows != correlation_rows
    ):
        raise BridgeError("global nuisance checkpoint lock differs from the exact V2 benchmark and task")
    envelope = gate_lock["resource_envelope"]
    if not 1 <= workers <= int(envelope["maximum_workers"]):
        raise BridgeError("benchmark worker count exceeds the V2 resource envelope")
    if required_workers is not None and workers != required_workers:
        raise BridgeError("full scan worker count must exactly match its admitted benchmark")
    method = row["rss_measurement_method"]
    if method not in {
        "AGGREGATE_PROCESS_TREE_PS_SAMPLED",
        "CHILD_RUSAGE_SESSION_MAXRSS_FALLBACK_PS_UNAVAILABLE",
    }:
        raise BridgeError("benchmark uses an unknown RSS measurement method")
    if workers > 1 and method != "AGGREGATE_PROCESS_TREE_PS_SAMPLED":
        raise BridgeError("multiworker scan lacks aggregate process-tree RSS evidence")
    reserve = int(envelope["ram_safety_reserve_bytes"])
    physical = ENGINE.physical_memory_bytes()
    if host_memory != physical:
        raise BridgeError("full scan host RAM differs from its matched benchmark host")
    if peak <= 0 or peak > host_memory or peak + reserve > physical:
        raise BridgeError(
            f"BLOCKED_BY_COMPUTE benchmark peak plus reserve exceeds RAM: "
            f"peak={peak} reserve={reserve} physical={physical}"
        )
    if (
        input_rows != int(task["aligned_input_rows"])
        or not 2 <= variance_rows <= input_rows
        or not 2 <= correlation_rows <= input_rows
        or not 1 <= benchmark_n <= input_rows
        or not 0 <= failures <= benchmark_n
        or row["global_nuisance_checkpoint_reused"] not in {"TRUE", "FALSE"}
        or not math.isfinite(nuisance_elapsed) or nuisance_elapsed < 0
        or not math.isfinite(testing_elapsed) or testing_elapsed <= 0
        or not math.isfinite(runner_elapsed) or runner_elapsed <= 0
        or not math.isfinite(rate) or rate <= 0
        or not math.isclose(rate, benchmark_n / testing_elapsed, rel_tol=1e-10, abs_tol=1e-12)
        or not math.isfinite(projected) or projected <= 0
        or not math.isclose(projected, input_rows / rate, rel_tol=1e-10, abs_tol=1e-12)
        or not math.isfinite(projected_hours) or projected_hours <= 0
        or not math.isclose(projected_hours, projected / 3600, rel_tol=1e-10, abs_tol=1e-12)
        or not math.isfinite(variance1) or variance1 <= 0
        or not math.isfinite(variance2) or variance2 <= 0
        or not math.isfinite(correlation) or not -1 < correlation < 1
        or not math.isfinite(wrapper_wall) or wrapper_wall <= 0
        or not math.isfinite(sample_interval) or not math.isclose(sample_interval, 0.05)
        or host_free < int(task["minimum_free_bytes"])
        or not math.isfinite(peak_fraction)
        or not math.isclose(
            peak_fraction, peak / host_memory, rel_tol=1e-10, abs_tol=1e-12,
        )
        or any(not row[key].strip() for key in ("r_version", "data_table_version", "runtime_platform"))
    ):
        raise BridgeError("benchmark lacks valid global-nuisance, timing, runtime, or RSS provenance")
    if shutil.disk_usage(ROOT).free < int(task["minimum_free_bytes"]):
        raise BridgeError("BLOCKED_BY_COMPUTE free storage is below the locked runner minimum")
    return row


def validate_stage_against_benchmark(
    paths: dict[str, Path], benchmark: dict[str, str], workers: int,
) -> None:
    try:
        summary = ENGINE.one_tsv_row(paths["run_summary"])
        observed_workers = int(summary["workers"])
        comparisons = {
            "nuisance_checkpoint_sha256": benchmark["nuisance_checkpoint_sha256"],
            "nuisance_null_rows_variance": benchmark["nuisance_null_rows_variance"],
            "nuisance_null_rows_correlation": benchmark["nuisance_null_rows_correlation"],
            "r_version": benchmark["r_version"],
            "data_table_version": benchmark["data_table_version"],
            "runtime_platform": benchmark["runtime_platform"],
        }
        if observed_workers != workers or any(summary.get(key) != value for key, value in comparisons.items()):
            raise BridgeError(
                "staged PLACO+ result does not match the admitted worker count, nuisance fit, and runtime"
            )
        for key in ("VarZ1", "VarZ2", "CorZ"):
            if not math.isclose(
                float(summary[key]), float(benchmark[key]), rel_tol=1e-12, abs_tol=1e-15,
            ):
                raise BridgeError("staged PLACO+ nuisance estimates differ from the admitted benchmark")
    except (KeyError, TypeError, ValueError, OverflowError, SystemExit) as error:
        if isinstance(error, BridgeError):
            raise
        raise BridgeError(f"invalid staged PLACO+ admission provenance: {error}") from error


def revalidate_publication_state(
    context: dict[str, Any],
    expected_gate: dict[str, Any],
    paths: dict[str, Path],
    workers: int,
) -> None:
    """Revalidate every mutable admission boundary at the actual link point."""
    validate_v2_paths(context)
    task = ENGINE.one_tsv_row(paths["task"])
    ENGINE.validate_task_against_context(task, context, paths["task"])
    benchmark = validate_benchmark(context, expected_gate, required_workers=workers)
    validate_stage_against_benchmark(paths, benchmark, workers)
    current_gate = GATE.verify_gate()
    if current_gate != expected_gate:
        raise BridgeError("PLACO V2 terminal gate changed during the full scan")
    # Keep this last: a passing or partial canonical LAVA family must be absent
    # at the closest possible instant before (and immediately after) linking.
    GATE.assert_no_canonical_lava_results()


def _run_full_scan_locked(
    context: dict[str, Any], gate_lock: dict[str, Any], workers: int,
) -> dict[str, Any]:
    validate_v2_paths(context)
    validate_benchmark(context, gate_lock, required_workers=workers)
    try:
        ENGINE.resource_preflight(context)
        paths = ENGINE.materialized_paths(context)
        ENGINE.verify_materialized(context)
    except SystemExit as error:
        raise BridgeError(str(error)) from error
    task_relative = ENGINE.relative(ROOT, paths["task"])
    command = [
        str(ROOT / ".r-env/bin/Rscript"), str(ROOT / ENGINE.RUNNER), task_relative,
        "--root", str(ROOT), "--execute", "--workers", str(workers), "--stage-only",
    ]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise BridgeError("V2 PLACO+ full-P runner failed:\n" + result.stdout + result.stderr)
    paths = ENGINE.materialized_paths(context)
    published: list[tuple[Path, int, int]] = []
    original_exclusive_family = ENGINE.exclusive_family

    def publication_boundary(sources: list[Path], destinations: list[Path]) -> None:
        exclusive_family_v2(
            sources,
            destinations,
            before_publication=lambda: revalidate_publication_state(
                context, gate_lock, paths, workers,
            ),
            after_publication=lambda: revalidate_publication_state(
                context, gate_lock, paths, workers,
            ),
            committed=published,
        )

    try:
        # The first pass rejects a bad stage before entering the V1 semantic
        # scan.  The injected publication boundary repeats it after that scan.
        revalidate_publication_state(context, gate_lock, paths, workers)
        ENGINE.exclusive_family = publication_boundary
        provenance = ENGINE.publish_run(ROOT, paths["task"], paths["run_summary"])
        if (
            provenance.get("reference_sha256")
            != "NOT_APPLICABLE_TO_PLACO_FULL_P_SCAN_LD_REQUIRED_FOR_LATER_CLUMPING"
        ):
            raise BridgeError("full-P result incorrectly claims LD-backed locus publication")
        if provenance.get("execution_resource_envelope", {}).get("observed_workers") != workers:
            raise BridgeError("published PLACO+ result worker count differs from its admitted benchmark")
        # Cover the V1 engine's final canonical-ledger identity check too.
        revalidate_publication_state(context, gate_lock, paths, workers)
    except (SystemExit, GATE.GateError) as error:
        try:
            rollback_v2_publication(published)
        except BridgeError as rollback_error:
            raise rollback_error from error
        raise BridgeError(str(error)) from error
    except BaseException as error:
        try:
            rollback_v2_publication(published)
        except BridgeError as rollback_error:
            raise rollback_error from error
        raise
    finally:
        ENGINE.exclusive_family = original_exclusive_family
    return provenance


def run_full_scan(context: dict[str, Any], gate_lock: dict[str, Any], workers: int) -> dict[str, Any]:
    with execution_lock("full-scan", context["pair_id"]):
        return _run_full_scan_locked(context, gate_lock, workers)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pair_id", choices=ENGINE.PAIR_IDS)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--preflight", action="store_true")
    action.add_argument("--materialize", action="store_true")
    action.add_argument("--verify-materialized", action="store_true")
    action.add_argument("--benchmark", action="store_true")
    action.add_argument("--verify-benchmark", action="store_true")
    action.add_argument("--run", action="store_true")
    parser.add_argument("--benchmark-variants", type=int, default=2000)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.materialize or args.benchmark or args.run:
        raise BridgeError(
            "direct PLACO V2 mutation is retired; use "
            "scripts/143_run_track_b_placo_sequential_v2.py so the inherited "
            "pair worker lease and live RAM evidence cannot be bypassed"
        )
    mutation = (
        "materialize" if args.materialize else
        "benchmark" if args.benchmark else
        "full-scan" if args.run else None
    )
    lock = execution_lock(mutation, args.pair_id) if mutation else nullcontext()
    with lock:
        context, gate_lock = configure_context(args.pair_id)
        if args.workers < 1 or args.workers > int(gate_lock["resource_envelope"]["maximum_workers"]):
            raise BridgeError("--workers must be between 1 and the locked maximum")
        paths = ENGINE.materialized_paths(context)
        if args.preflight:
            try:
                resources = ENGINE.resource_preflight(context)
            except SystemExit as error:
                raise BridgeError(str(error)) from error
            print(
                "TRACK_B_PLACO_V2_PREFLIGHT_PASS "
                f"pair={args.pair_id} fingerprint={context['fingerprint']} "
                f"ram={resources['physical_memory_bytes']} free={resources['free_bytes']} "
                "lava_evidence=FAILED_QC_NOT_CONSUMED"
            )
        elif args.materialize:
            try:
                validate_v2_paths(context)
                ENGINE.resource_preflight(context)
                provenance = ENGINE.materialize_pair(context)
                validate_v2_paths(context)
            except SystemExit as error:
                raise BridgeError(str(error)) from error
            print(
                f"TRACK_B_PLACO_V2_PAIR_MATERIALIZED pair={args.pair_id} "
                f"eligible={provenance['alignment_counts']['eligible_written']} "
                f"task={ENGINE.relative(ROOT, paths['task'])}"
            )
        elif args.verify_materialized:
            try:
                validate_v2_paths(context)
                provenance = ENGINE.verify_materialized(context)
            except SystemExit as error:
                raise BridgeError(str(error)) from error
            print(
                f"TRACK_B_PLACO_V2_PAIR_VERIFIED pair={args.pair_id} "
                f"eligible={provenance['alignment_counts']['eligible_written']}"
            )
        elif args.benchmark:
            try:
                validate_v2_paths(context)
                ENGINE.resource_preflight(context)
                row = ENGINE.benchmark_pair(context, args.benchmark_variants, args.workers)
                validate_v2_paths(context)
            except SystemExit as error:
                raise BridgeError(str(error)) from error
            validate_benchmark(context, gate_lock, required_workers=args.workers)
            print(
                f"TRACK_B_PLACO_V2_BENCHMARK_COMPLETE pair={args.pair_id} "
                f"workers={row['workers']} peak_rss={row['peak_process_tree_rss_bytes']} "
                f"projected_full_seconds={row['projected_full_family_seconds']}"
            )
        elif args.verify_benchmark:
            validate_v2_paths(context)
            row = validate_benchmark(context, gate_lock, required_workers=args.workers)
            print(
                f"TRACK_B_PLACO_V2_BENCHMARK_VERIFIED pair={args.pair_id} "
                f"workers={row['workers']} peak_rss={row['peak_process_tree_rss_bytes']}"
            )
        else:
            if not args.execute:
                raise BridgeError("explicit --execute is required for a full PLACO+ scan")
            provenance = _run_full_scan_locked(context, gate_lock, args.workers)
            print(
                f"TRACK_B_PLACO_V2_FULL_P_PUBLISHED pair={args.pair_id} "
                f"rows={provenance['output_rows']} status={provenance['terminal_result_state']} "
                "locus_publication=BLOCKED_PENDING_FULL_LD_AND_COLLATOR"
            )


if __name__ == "__main__":
    try:
        main()
    except BridgeError as error:
        raise SystemExit(f"ERROR: {error}") from error
