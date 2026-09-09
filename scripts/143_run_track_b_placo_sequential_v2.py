#!/usr/bin/env python3
"""Run the frozen Track B PLACO+ family sequentially within an 8-GiB envelope.

This additive coordinator deliberately reuses scripts 139--141 and the frozen
123--127 implementation.  It fixes the execution order to A, B, CONTROL and
the initial/admitted worker count to one.  A complete global nuisance fit is
therefore made on every valid aligned variant before the mathematically
independent single-variant PLACO+ calls are checkpointed in fixed shards.

The coordinator also owns the storage lifecycle.  It never removes an input,
canonical result, gate, source, or LD artifact.  After a published pair passes
the frozen engine's complete-ledger and BH validation, a compact immutable
reproducibility bundle is sealed.  Only then may an explicit whitelist of
regenerable work artifacts be unlinked.  Cleanup is itself restartable.
"""

from __future__ import annotations

import argparse
import base64
from contextlib import contextmanager
import csv
import errno
import fcntl
import hashlib
import importlib.util
import io
import json
import math
import os
import re
import resource
import signal
import shutil
import stat
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIDGE_PATH = ROOT / "scripts/140_materialize_track_b_placo_pair_v2.py"
ARCHIVE_STATE_GUARD_PATH = ROOT / "scripts/146_manage_lava_reference_archives.py"
COORDINATOR_PATH = Path(__file__).resolve()


class CoordinatorError(RuntimeError):
    """A fail-closed sequential execution or cleanup violation."""


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise CoordinatorError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


BRIDGE = load_module("track_b_placo_v2_bridge_for_sequential_coordinator", BRIDGE_PATH)
ENGINE = BRIDGE.ENGINE
GATE = BRIDGE.GATE
ARCHIVE_STATE_GUARD = load_module(
    "track_b_lava_reference_archive_guard_for_placo_coordinator",
    ARCHIVE_STATE_GUARD_PATH,
)

ARCHIVES_PRESENT = "ARCHIVES_PRESENT_FULLY_VERIFIED"
ARCHIVES_EVICTED = (
    "ARCHIVES_INTENTIONALLY_EVICTED_WITH_FINAL_RECEIPT_AND_ALL_44_"
    "EXTRACTED_PAYLOADS_REHASHED"
)
ACCEPTED_ARCHIVE_STATES = {ARCHIVES_PRESENT, ARCHIVES_EVICTED}

PAIR_ORDER = ("A", "B", "CONTROL")
WORKERS = 1
BENCHMARK_VARIANTS = 2000
RAM_LIMIT_BYTES = 8 * 1024**3
EXECUTION_CONTRACT = GATE.GATE_ROOT / "sequential_execution.contract.json"
RECORD_ROOT = GATE.GATE_ROOT / "sequential_execution_records"
WORKER_LEASE_ROOT = GATE.GATE_ROOT / "sequential_worker_leases"
RAM_AGGREGATE = GATE.GATE_ROOT / "RAM_BENCHMARK.tsv"
RAM_AGGREGATE_PROVENANCE = GATE.GATE_ROOT / "RAM_BENCHMARK.provenance.json"
CONTRACT_SCHEMA = "sleep-atlas-track-b-placo-sequential-execution.1"
REPRODUCIBILITY_SCHEMA = "sleep-atlas-track-b-placo-cleanup-bundle.1"
CLEANUP_COMPLETE_SCHEMA = "sleep-atlas-track-b-placo-cleanup-complete.1"
RAM_FIELDS = [
    "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
    "runtime_sec", "exit_status", "output_hash",
]
PRESERVED_WORK_ARTIFACTS = {
    "materialization.provenance.json",
    "task.tsv",
    "benchmark.tsv",
    "run.summary.tsv",
    "checkpoints/nuisance.rds",
    "checkpoints/nuisance.sha256.tsv",
}
OPTIONAL_PRESERVED_WORK_ARTIFACTS = {"benchmark.raw.tsv"}
FIXED_REGENERABLE_ARTIFACTS = {
    "aligned.tsv.gz",
    "materialization.provenance.json",
    "task.tsv",
    "alignment.sqlite",
    "alignment.sqlite-journal",
    "benchmark.raw.tsv",
    "benchmark.tsv",
    "staged.full.tsv.gz",
    "run.summary.tsv",
    "staged.provenance.json",
    "checkpoints/nuisance.rds",
    "checkpoints/nuisance.sha256.tsv",
}
SHARD_RE = re.compile(r"checkpoints/shard_[0-9]{6}\.rds(?:\.sha256\.tsv)?$")
KNOWN_TEMP_RE = re.compile(
    r"(?:"
    r"\.(?:aligned\.[0-9]+\.tmp\.gz|provenance\.[0-9]+\.tmp\.json|task\.[0-9]+\.tmp\.tsv)"
    r"|\.(?:benchmark\.tsv|staged\.full\.tsv\.gz|staged\.provenance\.json)\.[0-9]+\.tmp"
    r"|(?:benchmark\.raw\.tsv|staged\.full\.tsv\.gz|run\.summary\.tsv)\.[0-9]+\.tmp"
    r"|checkpoints/(?:nuisance\.rds|nuisance\.sha256\.tsv|shard_[0-9]{6}\.rds(?:\.sha256\.tsv)?)\.[0-9]+\.tmp"
    r"|\.staged\.full\.tsv\.gz\.[0-9]+\.validate\.sqlite(?:-journal|-wal|-shm)?"
    r")$"
)
KNOWN_WORK_DIRECTORIES = {"checkpoints"}
RSS_MEASUREMENT_METHODS = {
    "AGGREGATE_PROCESS_TREE_PS_SAMPLED",
    "CHILD_RUSAGE_SESSION_MAXRSS_FALLBACK_PS_UNAVAILABLE",
    "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE",
}
MAX_EMBEDDED_PROVENANCE_BYTES = 64 * 1024**2
FULL_SCAN_ANALYSIS = "PLACO_PLUS_FULL_P_SCAN_GLOBAL_NUISANCE_FROZEN_COMPLETE_PAIR_BH"
MATCHED_BENCHMARK_ANALYSIS = "PLACO_PLUS_MATCHED_WORKER_BENCHMARK_GLOBAL_NUISANCE_ALL_VALID_VARIANTS"
FULL_SCAN_REJECTION_PREFIX = "full_scan_ram_rejected_"
ACTIVE_WORKER_LEASE_FD: int | None = None
ACTIVE_WORKER_LEASE_CONTEXT: dict[str, Any] | None = None
WORKER_LEASE_ENV = "TRACK_B_PLACO_WORKER_LEASE_FD"
WORKER_LEASE_PAYLOAD_ENV = "TRACK_B_PLACO_WORKER_LEASE_PAYLOAD_B64"
WORKER_LEASE_DIGEST_ENV = "TRACK_B_PLACO_WORKER_LEASE_PAYLOAD_SHA256"
WORKER_LEASE_SCHEMA = "sleep-atlas-track-b-placo-worker-lease.1"
FAILURE_EVENT_SCHEMA = "sleep-atlas-track-b-placo-failure-event.1"
FAILURE_EVENT_RE = re.compile(r"^failure_event\.([0-9a-f]{64})\.json$")
FAILURE_EVENT_CORE_FIELDS = frozenset({
    "analysis", "analysis_id", "attempt_log", "coordinator_sha256", "details",
    "failure_class", "input_sha256", "n_snps", "output_hash",
    "owned_child_exit_status", "pair_id", "peak_process_tree_rss_bytes",
    "phase_prefix", "process_group_termination_status", "public_exit_status",
    "rss_measurement_method", "run_fingerprint", "runtime_seconds",
    "sequential_contract", "task", "terminal_compute_rejection",
    "terminal_gate_lock", "termination_reason",
})
PROCESS_GROUP_TERMINATION_STATUSES = frozenset({
    "CONFIRMED", "NOT_REQUIRED_PROCESS_ALREADY_EXITED", "NOT_STARTED", "UNCONFIRMED",
})


def coordinator_sha256() -> str:
    return hashlib.sha256(COORDINATOR_PATH.read_bytes()).hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def tsv_bytes(fields: list[str], rows: list[dict[str, object]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode()


def safe_repo_path(path: Path, label: str) -> Path:
    try:
        path.relative_to(ROOT)
        root_resolved = ROOT.resolve(strict=True)
    except (OSError, RuntimeError, ValueError) as error:
        raise CoordinatorError(f"{label} escapes the repository: {path}") from error
    current = ROOT
    for part in path.relative_to(ROOT).parts:
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        if stat.S_ISLNK(observed.st_mode):
            raise CoordinatorError(f"{label} contains a symbolic link: {current}")
    try:
        path.resolve(strict=False).relative_to(root_resolved)
    except (OSError, RuntimeError, ValueError) as error:
        raise CoordinatorError(f"{label} escapes the repository: {path}") from error
    return path


def stable_identity(path: Path, *, allow_empty: bool = False) -> dict[str, object]:
    safe_repo_path(path, "artifact")
    digest = hashlib.sha256()
    try:
        before_path = os.lstat(path)
        if (
            not stat.S_ISREG(before_path.st_mode)
            or before_path.st_size < 0
            or (before_path.st_size == 0 and not allow_empty)
        ):
            qualifier = "regular file" if allow_empty else "non-empty regular file"
            raise CoordinatorError(f"artifact is not a real {qualifier}: {path}")
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = os.lstat(path)
    except OSError as error:
        raise CoordinatorError(f"could not hash artifact {path}: {error}") from error
    identities = [
        (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns)
        for value in (before_path, before, after, current)
    ]
    if len(set(identities)) != 1 or not stat.S_ISREG(current.st_mode):
        raise CoordinatorError(f"artifact changed while hashing: {path}")
    return {
        "path": str(path.relative_to(ROOT)),
        "device": int(current.st_dev),
        "inode": int(current.st_ino),
        "bytes": int(current.st_size),
        "sha256": digest.hexdigest(),
    }


def identity_matches(
    path: Path, expected: dict[str, Any], *, include_inode: bool = True,
    allow_empty: bool = False,
) -> bool:
    observed = stable_identity(path, allow_empty=allow_empty)
    keys = ["path", "bytes", "sha256"]
    if include_inode:
        keys.extend(["device", "inode"])
    return all(observed.get(key) == expected.get(key) for key in keys)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def worker_lease_path(pair_id: str) -> Path:
    if pair_id not in PAIR_ORDER:
        raise CoordinatorError(f"unknown PLACO+ worker-lease pair: {pair_id}")
    return ROOT / WORKER_LEASE_ROOT / f"{pair_id}.lock"


@contextmanager
def worker_lease(context: dict[str, Any]):
    """Hold a pair-stable flock that the owned R process inherits across exec.

    The lease path deliberately does not include the run fingerprint: an orphan
    from an older fingerprint must also block a new run.  `pass_fds` keeps this
    lock alive if the Python coordinator itself is SIGKILLed, so a restart can
    never classify a live R temporary as stale or launch a second pair worker.
    """
    global ACTIVE_WORKER_LEASE_FD, ACTIVE_WORKER_LEASE_CONTEXT
    if ACTIVE_WORKER_LEASE_FD is not None or ACTIVE_WORKER_LEASE_CONTEXT is not None:
        raise CoordinatorError("nested PLACO+ worker leases are forbidden")
    pair_id = str(context.get("pair_id", ""))
    run_fingerprint = str(context.get("fingerprint", ""))
    if pair_id not in PAIR_ORDER or not re.fullmatch(r"[0-9a-f]{64}", run_fingerprint):
        raise CoordinatorError("PLACO+ worker lease requires an exact pair and run fingerprint")
    path = worker_lease_path(pair_id)
    safe_repo_path(path.parent, "PLACO+ worker lease parent")
    path.parent.mkdir(parents=True, exist_ok=True)
    safe_repo_path(path.parent, "PLACO+ worker lease parent")
    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(path, flags, 0o600)
    except OSError as error:
        raise CoordinatorError(f"could not open PLACO+ worker lease: {path}: {error}") from error
    try:
        opened = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            not stat.S_ISREG(opened.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (opened.st_dev, opened.st_ino) != (current.st_dev, current.st_ino)
            or current.st_nlink != 1
        ):
            raise CoordinatorError("PLACO+ worker lease is not a private regular file")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno in {errno.EACCES, errno.EAGAIN}:
                raise CoordinatorError(
                    f"active or orphaned PLACO+ worker still owns pair {context['pair_id']}"
                ) from error
            raise CoordinatorError(f"could not lock PLACO+ worker lease: {error}") from error
        locked = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            not stat.S_ISREG(locked.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (locked.st_dev, locked.st_ino) != (current.st_dev, current.st_ino)
            or current.st_nlink != 1
        ):
            raise CoordinatorError("PLACO+ worker lease changed while acquiring its flock")
        payload = {
            "schema_version": WORKER_LEASE_SCHEMA,
            "pair_id": pair_id,
            "run_fingerprint": run_fingerprint,
            "device": int(locked.st_dev),
            "inode": int(locked.st_ino),
            "authorization_nonce": os.urandom(32).hex(),
            "coordinator_sha256": coordinator_sha256(),
        }
        encoded = canonical_json_bytes(payload) + b"\n"
        try:
            os.ftruncate(descriptor, 0)
            os.lseek(descriptor, 0, os.SEEK_SET)
            offset = 0
            while offset < len(encoded):
                written = os.write(descriptor, encoded[offset:])
                if written <= 0:
                    raise OSError("short write while authorizing worker lease")
                offset += written
            os.fsync(descriptor)
            fsync_directory(path.parent)
        except OSError as error:
            raise CoordinatorError(f"could not authorize PLACO+ worker lease: {error}") from error
        after_write = os.fstat(descriptor)
        after_path = os.lstat(path)
        if (
            not stat.S_ISREG(after_write.st_mode)
            or not stat.S_ISREG(after_path.st_mode)
            or (after_write.st_dev, after_write.st_ino) != (after_path.st_dev, after_path.st_ino)
            or after_path.st_nlink != 1
            or after_write.st_size != len(encoded)
        ):
            raise CoordinatorError("PLACO+ worker lease changed while authorizing its payload")
        digest = hashlib.sha256(encoded).hexdigest()
        ACTIVE_WORKER_LEASE_FD = descriptor
        ACTIVE_WORKER_LEASE_CONTEXT = {
            "pair_id": pair_id,
            "run_fingerprint": run_fingerprint,
            "path": path,
            "device": int(after_write.st_dev),
            "inode": int(after_write.st_ino),
            "encoded": encoded,
            "payload_b64": base64.b64encode(encoded).decode("ascii"),
            "payload_sha256": digest,
        }
        yield descriptor
    finally:
        if ACTIVE_WORKER_LEASE_FD == descriptor:
            ACTIVE_WORKER_LEASE_FD = None
            ACTIVE_WORKER_LEASE_CONTEXT = None
        # Never issue LOCK_UN: the R child shares this open-file description.
        # Closing only our descriptor releases the lease iff no child survived.
        os.close(descriptor)


def require_active_worker_lease(context: dict[str, Any]) -> dict[str, Any]:
    active = ACTIVE_WORKER_LEASE_CONTEXT
    descriptor = ACTIVE_WORKER_LEASE_FD
    if active is None or descriptor is None:
        raise CoordinatorError("PLACO+ action attempted without its durable worker lease")
    if (
        str(context.get("pair_id", "")) != active["pair_id"]
        or str(context.get("fingerprint", "")) != active["run_fingerprint"]
    ):
        raise CoordinatorError("active PLACO+ worker lease belongs to a different pair or run fingerprint")
    path = Path(active["path"])
    try:
        opened = os.fstat(descriptor)
        current = os.lstat(path)
        content = os.pread(descriptor, opened.st_size, 0)
        after_opened = os.fstat(descriptor)
        after_current = os.lstat(path)
    except OSError as error:
        raise CoordinatorError(f"could not revalidate active PLACO+ worker lease: {error}") from error
    if (
        not stat.S_ISREG(opened.st_mode)
        or not stat.S_ISREG(current.st_mode)
        or not stat.S_ISREG(after_opened.st_mode)
        or not stat.S_ISREG(after_current.st_mode)
        or (opened.st_dev, opened.st_ino) != (active["device"], active["inode"])
        or (current.st_dev, current.st_ino) != (active["device"], active["inode"])
        or (after_opened.st_dev, after_opened.st_ino) != (active["device"], active["inode"])
        or (after_current.st_dev, after_current.st_ino) != (active["device"], active["inode"])
        or opened.st_size != after_opened.st_size
        or current.st_size != after_current.st_size
        or current.st_nlink != 1
        or after_current.st_nlink != 1
        or content != active["encoded"]
        or hashlib.sha256(content).hexdigest() != active["payload_sha256"]
    ):
        raise CoordinatorError("active PLACO+ worker lease identity or authorization payload drifted")
    return active


def inherited_worker_lease_fds(existing: Any = ()) -> tuple[int, ...]:
    if ACTIVE_WORKER_LEASE_FD is None or ACTIVE_WORKER_LEASE_CONTEXT is None:
        raise CoordinatorError("PLACO+ R launch attempted without its durable worker lease")
    try:
        descriptors = tuple(int(value) for value in existing)
    except (TypeError, ValueError) as error:
        raise CoordinatorError("PLACO+ R launch has malformed inherited descriptors") from error
    return tuple(dict.fromkeys((*descriptors, ACTIVE_WORKER_LEASE_FD)))


def inherited_worker_environment(existing: Any = None) -> dict[str, str]:
    if ACTIVE_WORKER_LEASE_FD is None or ACTIVE_WORKER_LEASE_CONTEXT is None:
        raise CoordinatorError("PLACO+ R launch attempted without its durable worker lease")
    environment = dict(os.environ if existing is None else existing)
    # The frozen R worker invokes a tiny Python lease validator.  External
    # Python warning/debug settings must not change that validator's protocol
    # (for example, turning its otherwise-ignored ResourceWarning into a
    # second output line and rejecting a valid lease).
    environment.pop("PYTHONWARNINGS", None)
    environment.pop("PYTHONDEVMODE", None)
    environment[WORKER_LEASE_ENV] = str(ACTIVE_WORKER_LEASE_FD)
    environment[WORKER_LEASE_PAYLOAD_ENV] = str(ACTIVE_WORKER_LEASE_CONTEXT["payload_b64"])
    environment[WORKER_LEASE_DIGEST_ENV] = str(ACTIVE_WORKER_LEASE_CONTEXT["payload_sha256"])
    return environment


PUBLICATION_TEMP_WITH_HASH_RE = re.compile(
    r"^\.(?P<target>.+)\.(?P<sha256>[0-9a-f]{64})\.(?P<pid>[0-9]+)\.tmp$"
)
PUBLICATION_TEMP_LEGACY_RE = re.compile(
    r"^\.(?P<target>.+)\.(?P<pid>[0-9]+)\.tmp$"
)


def publish_temp_reconciliation_audit(records: list[dict[str, object]]) -> None:
    if not records:
        return
    payload = {
        "schema_version": "sleep-atlas-track-b-placo-publication-temp-reconciliation.1",
        "artifacts": records,
        "rule": "HASHED_COMPLETE_TEMP_MAY_BE_PROMOTED;EXACT_ALIAS_MAY_BE_UNLINKED;OTHERWISE_FAIL_CLOSED",
        "scientific_result_deleted": False,
        "coordinator_sha256": coordinator_sha256(),
    }
    content = json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n"
    digest = hashlib.sha256(content).hexdigest()
    path = ROOT / GATE.GATE_ROOT / "sequential_temp_reconciliation" / f"receipt.{digest}.json"
    safe_repo_path(path.parent, "publication-temp reconciliation audit parent")
    path.parent.mkdir(parents=True, exist_ok=True)
    safe_repo_path(path.parent, "publication-temp reconciliation audit parent")
    if artifact_present(path):
        identity = stable_identity(path)
        if identity["bytes"] != len(content) or identity["sha256"] != digest:
            raise CoordinatorError("publication-temp reconciliation audit drifted")
        return
    prefix = f".{path.name}."
    for stale_audit in sorted(path.parent.iterdir()):
        if not stale_audit.name.startswith(prefix) or not stale_audit.name.endswith(".tmp"):
            continue
        if re.fullmatch(re.escape(prefix) + r"[0-9a-f]{64}\.[0-9]+\.tmp", stale_audit.name) is None:
            raise CoordinatorError("malformed publication-temp audit temporary")
        observed = os.lstat(stale_audit)
        if not stat.S_ISREG(observed.st_mode):
            raise CoordinatorError("publication-temp audit temporary is not regular")
        stale_audit.unlink()
    temporary = path.with_name(f".{path.name}.{digest}.{os.getpid()}.tmp")
    with temporary.open("xb") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    if hashlib.sha256(temporary.read_bytes()).hexdigest() != digest:
        raise CoordinatorError("publication-temp audit temporary changed before commit")
    os.rename(temporary, path)
    fsync_directory(path.parent)


def reconcile_publication_temporaries_for_target(
    path: Path, content: bytes, *, derived: bool = False,
) -> None:
    if not path.parent.exists():
        return
    wanted_hash = hashlib.sha256(content).hexdigest()
    records: list[dict[str, object]] = []
    candidates: list[tuple[Path, dict[str, object], str]] = []
    prefix = f".{path.name}."
    for temporary in sorted(path.parent.iterdir()):
        if not temporary.name.startswith(prefix) or not temporary.name.endswith(".tmp"):
            continue
        hashed = PUBLICATION_TEMP_WITH_HASH_RE.fullmatch(temporary.name)
        legacy = PUBLICATION_TEMP_LEGACY_RE.fullmatch(temporary.name)
        match = hashed or legacy
        if match is None or match.group("target") != path.name:
            raise CoordinatorError(f"malformed publication temporary: {temporary}")
        identity = stable_identity(temporary, allow_empty=True)
        embedded_hash = hashed.group("sha256") if hashed is not None else None
        matches_current = identity["sha256"] == wanted_hash and identity["bytes"] == len(content)
        hash_complete = (
            (embedded_hash is not None and identity["sha256"] == embedded_hash)
            or (embedded_hash is None and matches_current)
        )
        if not hash_complete:
            action = "DELETE_TORN_OR_LEGACY_TEMP"
        elif derived:
            action = "DELETE_DERIVED_TEMP"
        elif not matches_current:
            action = "DELETE_STALE_COMPLETE_TEMP"
        else:
            action = "UNLINK_EXACT_ALIAS" if artifact_present(path) else "PROMOTE_COMPLETE_TEMP"
        if artifact_present(path):
            destination = stable_identity(path, allow_empty=True)
            if not derived and action == "UNLINK_EXACT_ALIAS" and (
                destination["bytes"] != len(content) or destination["sha256"] != wanted_hash
            ):
                action = "DELETE_TEMP_COEXISTING_WITH_DRIFTED_FINAL"
        elif embedded_hash is None and derived:
            # Legacy derived temps are safe to discard because the aggregate is
            # reconstructed from immutable source receipts.
            action = "DELETE_LEGACY_DERIVED_TEMP"
        candidates.append((temporary, identity, action))
        records.append({**identity, "target": str(path.relative_to(ROOT)), "action": action})
    publish_temp_reconciliation_audit(records)
    for temporary, identity, action in candidates:
        current = os.lstat(temporary)
        if (current.st_dev, current.st_ino, current.st_size) != (
            identity["device"], identity["inode"], identity["bytes"],
        ):
            raise CoordinatorError("publication temporary changed after reconciliation audit")
        if action == "PROMOTE_COMPLETE_TEMP":
            try:
                os.link(temporary, path, follow_symlinks=False)
            except FileExistsError as error:
                raise CoordinatorError("publication destination appeared during temp promotion") from error
            destination = os.lstat(path)
            if (destination.st_dev, destination.st_ino, destination.st_size) != (
                current.st_dev, current.st_ino, current.st_size,
            ):
                raise CoordinatorError("promoted publication temporary is not the exact inode")
            fsync_directory(path.parent)
        temporary.unlink()
        fsync_directory(temporary.parent)


def publish_bytes_no_replace(path: Path, content: bytes) -> None:
    safe_repo_path(path.parent, "publication parent")
    path.parent.mkdir(parents=True, exist_ok=True)
    safe_repo_path(path.parent, "publication parent")
    safe_repo_path(path, "publication destination")
    reconcile_publication_temporaries_for_target(path, content)
    if path.exists() or path.is_symlink():
        if stable_identity(path)["sha256"] != hashlib.sha256(content).hexdigest() or path.stat().st_size != len(content):
            raise CoordinatorError(f"immutable publication differs from expected bytes: {path}")
        return
    content_hash = hashlib.sha256(content).hexdigest()
    temporary = path.with_name(f".{path.name}.{content_hash}.{os.getpid()}.tmp")
    safe_repo_path(temporary, "publication temporary")
    published: tuple[int, int] | None = None
    try:
        with temporary.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        source = os.lstat(temporary)
        try:
            os.link(temporary, path, follow_symlinks=False)
        except FileExistsError as error:
            raise CoordinatorError(f"immutable destination appeared concurrently: {path}") from error
        published = (source.st_dev, source.st_ino)
        destination = os.lstat(path)
        if (
            not stat.S_ISREG(destination.st_mode)
            or (source.st_dev, source.st_ino) != (destination.st_dev, destination.st_ino)
        ):
            raise CoordinatorError(f"no-replace publication identity mismatch: {path}")
        fsync_directory(path.parent)
    except BaseException:
        if published is not None:
            try:
                current = os.lstat(path)
                if (current.st_dev, current.st_ino) == published:
                    path.unlink()
                    fsync_directory(path.parent)
            except FileNotFoundError:
                pass
        raise
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def read_json(path: Path) -> dict[str, Any]:
    safe_repo_path(path, "JSON artifact")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CoordinatorError(f"unreadable JSON artifact {path}: {error}") from error
    if not isinstance(value, dict):
        raise CoordinatorError(f"JSON artifact must contain one object: {path}")
    return value


def reference_archive_state_invariant(reference_state: dict[str, Any]) -> dict[str, Any]:
    guard_identity = stable_identity(ARCHIVE_STATE_GUARD_PATH)
    return {
        "accepted_states": [ARCHIVES_PRESENT, ARCHIVES_EVICTED],
        "guard_implementation": {
            key: guard_identity[key] for key in ("path", "bytes", "sha256")
        },
        "reference_provenance": reference_state["reference_provenance_identity"],
        "archive_family_sha256": reference_state["archive_family_sha256"],
        "extracted_payload_count": reference_state["extracted_payload_count"],
        "extracted_family_sha256": reference_state["extracted_family_sha256"],
        "verification": (
            "CURRENT_STATE_REVALIDATED_EVERY_CALL;MUTABLE_PRESENT_OR_EVICTED_TOKEN_"
            "INTENTIONALLY_EXCLUDED_FROM_THE_IMMUTABLE_TERMINAL_GATE"
        ),
    }


def verify_reference_archive_state(gate_lock: dict[str, Any]) -> dict[str, Any]:
    """Rehash the reference and bind it to the state-invariant terminal gate."""

    try:
        observed = ARCHIVE_STATE_GUARD.verify_reference_state(ROOT)
    except (ARCHIVE_STATE_GUARD.ArchiveEvictionError, OSError) as error:
        raise CoordinatorError(f"LAVA reference archive state failed verification: {error}") from error
    if not isinstance(observed, dict) or observed.get("state") not in ACCEPTED_ARCHIVE_STATES:
        raise CoordinatorError("LAVA reference archive guard returned a non-permitted state")
    if (
        observed.get("extracted_payload_count") != 44
        or not re.fullmatch(r"[0-9a-f]{64}", str(observed.get("archive_family_sha256", "")))
        or not re.fullmatch(r"[0-9a-f]{64}", str(observed.get("extracted_family_sha256", "")))
        or not isinstance(observed.get("reference_provenance_identity"), dict)
    ):
        raise CoordinatorError("LAVA reference archive guard returned incomplete scientific identity")
    try:
        frozen = gate_lock["lava_terminal"]["reference_archive_state_guard"]
    except (KeyError, TypeError) as error:
        raise CoordinatorError("terminal gate lacks the LAVA reference-state invariant") from error
    if frozen != reference_archive_state_invariant(observed):
        raise CoordinatorError("current LAVA reference family differs from the terminal gate")
    return observed


def execution_contract_payload(
    gate_lock: dict[str, Any], reference_state: dict[str, Any],
) -> dict[str, Any]:
    code = [
        Path(f"scripts/{number}_{name}")
        for number, name in (
            (123, "track_b_pleiotropy_contract.py"),
            (124, "build_track_b_pleiotropy_input_gate.py"),
            (125, "materialize_track_b_placo_pair.py"),
            (126, "run_track_b_placo_pair.R"),
            (127, "prepare_track_b_pleiotropy_ld.py"),
            (139, "build_track_b_placo_terminal_gate_v2.py"),
            (140, "materialize_track_b_placo_pair_v2.py"),
            (141, "run_track_b_placo_pair_v2.R"),
            (143, "run_track_b_placo_sequential_v2.py"),
            (146, "manage_lava_reference_archives.py"),
        )
    ]
    return {
        "schema_version": CONTRACT_SCHEMA,
        "analysis_id": gate_lock["analysis_id"],
        "pair_order": list(PAIR_ORDER),
        "pair_concurrency": 1,
        "workers": WORKERS,
        "worker_selection": "START_AND_REMAIN_AT_ONE;FULL_SCAN_MUST_MATCH_BENCHMARK",
        "benchmark_variants": BENCHMARK_VARIANTS,
        "benchmark_selection": "DETERMINISTIC_EVENLY_SPACED_INDICES_ACROSS_COMPLETE_ALIGNED_FAMILY",
        "benchmark_global_nuisance": "MUST_BE_FRESHLY_COMPUTED_ON_ALL_VALID_VARIANTS_AND_MEASURED",
        "ram_limit_bytes": RAM_LIMIT_BYTES,
        "ram_safety_reserve_bytes": int(gate_lock["resource_envelope"]["ram_safety_reserve_bytes"]),
        "ram_admission_rule": (
            "LIVE_OWNED_FRESH_SESSION_PGID_RSS_GUARD_INCLUDING_REPARENTED_DESCENDANTS;"
            "UNMEASURABLE_LIVE_RSS_TERMINATES_FAIL_CLOSED;"
            "MEASURED_PEAK_PLUS_FROZEN_RESERVE_MUST_NOT_EXCEED_MIN_8_GIB_AND_HOST_RAM;"
            "CANONICAL_LEDGER_REQUIRES_EXACT_SUCCESSFUL_FULL_SCAN_RAM_RECEIPT"
        ),
        "nuisance_scope": "ALL_VALID_ALIGNED_GENOME_WIDE_VARIANTS_PER_PAIR_BEFORE_SINGLE_VARIANT_TESTING",
        "test_decomposition": "ONLY_IDENTICAL_SINGLE_VARIANT_PLACO_PLUS_CALLS_AFTER_GLOBAL_NUISANCE_IS_FROZEN",
        "within_pair_bh": "ONCE_OVER_THE_COMPLETE_FIXED_ALIGNED_PAIR_FAMILY_WITH_FAILURE_P_EQUAL_ONE",
        "restart": (
            "PAIR_STABLE_FLOCK_INHERITED_BY_R_ACROSS_COORDINATOR_DEATH;"
            "COORDINATOR_AUTHORIZATION_BINDS_PAIR_RUN_FINGERPRINT_DEVICE_AND_INODE;"
            "WORKER_BOUND_NUISANCE_AND_RSS_EVIDENCE_BOUND_SHARD_CHECKPOINTS;"
            "UNMEASURED_SHARDS_AUDITED_AND_RECOMPUTED;DOCUMENTED_ZERO_OR_NONZERO_CRASH_TEMPS_AUDITED;"
            "EXACT_CANONICAL_HARDLINK_PREFIX_RECONCILED;NO_COMPLETED_PAIR_RECOMPUTATION"
        ),
        "storage": {
            "materialization_admission": "8_GIB_PLUS_600_BYTES_TIMES_SUM_OF_PAIR_SOURCE_ROWS",
            "post_materialization_runner_floor": "LOCKED_TASK_MINIMUM_FREE_BYTES_8_GIB",
            "reason": "DO_NOT_REAPPLY_MATERIALIZATION_TEMPORARY_SPACE_AFTER_EXACT_ALIGNED_FAMILY_EXISTS",
            "dynamic_recheck_before_each_stage": True,
            "sequential_forecast": "ACTUAL_RETAINED_CANONICAL_BYTES_PLUS_EXPLICIT_PENDING_PRIOR_LEDGER_BUDGET",
        },
        "cleanup": {
            "allowed_only_after": "DEEP_COMPLETE_LEDGER_VALIDATION_AND_IMMUTABLE_REPRODUCIBILITY_BUNDLE",
            "fixed_regenerable_paths": sorted(FIXED_REGENERABLE_ARTIFACTS),
            "shard_pattern": SHARD_RE.pattern,
            "known_work_directories": sorted(KNOWN_WORK_DIRECTORIES),
            "unknown_path_policy": "FAIL_CLOSED_NO_DELETE",
            "post_authorization_inventory_recheck": True,
            "reclaim_accounting": "UNIQUE_INODES_AND_EXTERNAL_HARDLINKS_AT_AUTHORIZATION",
            "canonical_results_deleted": False,
            "dense_sources_deleted": False,
        },
        "ld": {
            "extraction_during_full_p_scans": False,
            "locus_publication": "BLOCKED_UNTIL_ALL_FULL_P_SCANS_AND_SEPARATE_LD_COLLATOR",
        },
        "ram_benchmark_fields": RAM_FIELDS,
        "ram_receipt_rule": (
            "EXACT_JSON_BYTES_RUNTIME_STATUS_METHOD_AND_HASH_MUST_RECONSTRUCT_THE_TSV_ROW;"
            "ZERO_RSS_ALLOWED_ONLY_FOR_NONZERO_EXIT_BEFORE_ANY_PROCESS_RSS_WAS_OBSERVED"
        ),
        "terminal_gate_lock": stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK),
        "lava_reference_archive_state_guard": reference_archive_state_invariant(reference_state),
        "code": {str(path): stable_identity(ROOT / path) for path in code},
    }


def verify_execution_contract() -> dict[str, Any]:
    try:
        gate_lock = GATE.verify_gate()
    except (GATE.GateError, OSError) as error:
        raise CoordinatorError(str(error)) from error
    reference_state = verify_reference_archive_state(gate_lock)
    expected = execution_contract_payload(gate_lock, reference_state)
    observed = read_json(ROOT / EXECUTION_CONTRACT)
    if observed != expected:
        raise CoordinatorError("sequential PLACO+ execution contract is missing or drifted")
    return observed


def seal_execution_contract() -> dict[str, Any]:
    try:
        lock = BRIDGE.execution_lock("seal-sequential-contract", "A,B,CONTROL")
        with lock:
            gate_lock = GATE.verify_gate()
            reference_state = verify_reference_archive_state(gate_lock)
            policy = GATE.PLEIOTROPY.validate_policy(ROOT)
            evidence = GATE.placo_result_evidence(policy)
            path = ROOT / EXECUTION_CONTRACT
            if path.exists() or path.is_symlink():
                return verify_execution_contract()
            if evidence:
                raise CoordinatorError(
                    "refusing to freeze the sequential contract after PLACO result or checkpoint access: "
                    + ", ".join(str(item.relative_to(ROOT)) for item in evidence[:3])
                )
            payload = execution_contract_payload(gate_lock, reference_state)
            publish_bytes_no_replace(path, json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n")
            return verify_execution_contract()
    except (GATE.GateError, BRIDGE.BridgeError, SystemExit, OSError) as error:
        raise CoordinatorError(str(error)) from error


def pair_required_storage_bytes(context: dict[str, Any]) -> int:
    source_rows = sum(int(source["rows"]) for source in context["sources"])
    return int(ENGINE.MINIMUM_FREE_BYTES + 600 * source_rows)


def pending_canonical_planning_reserve_bytes(context: dict[str, Any]) -> int:
    """Conservative planning allowance; runtime admission still uses actual free bytes."""
    source_rows = sum(int(source["rows"]) for source in context["sources"])
    return 600 * source_rows


def artifact_present(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def forecast_with_required(free_bytes: int, required_bytes: int) -> dict[str, int]:
    return {
        "free_bytes": int(free_bytes),
        "required_bytes": int(required_bytes),
        "additional_bytes_required": max(0, int(required_bytes) - int(free_bytes)),
        "headroom_bytes": max(0, int(free_bytes) - int(required_bytes)),
    }


def pair_preflight_storage_forecast(
    context: dict[str, Any], free_bytes: int,
) -> tuple[dict[str, int], str]:
    """Forecast only storage not already allocated in a restartable pair state.

    This is a planning view, not a substitute for the authoritative runtime
    admission and deep validation.  Existing canonical/publication prefixes and
    cleanup receipts either resume without allocating another dense family or
    fail closed; they must not be charged as though materialization starts over.
    """
    cleanup_done = cleanup_complete_path(context)
    if artifact_present(cleanup_done):
        stable_identity(cleanup_done)
        return forecast_with_required(free_bytes, 0), "CLEANUP_COMPLETE"
    cleanup_plan = cleanup_bundle_path(context)
    if artifact_present(cleanup_plan):
        stable_identity(cleanup_plan)
        return forecast_with_required(free_bytes, 0), "CLEANUP_PLAN_READY"
    paths = ENGINE.materialized_paths(context)
    publication_state = canonical_state(paths)
    if publication_state != "ABSENT":
        for key in ("canonical_ledger", "canonical_provenance"):
            path = paths[key]
            if artifact_present(path):
                stable_identity(path)
        return (
            forecast_with_required(free_bytes, 0),
            f"CANONICAL_{publication_state}_REQUIRES_EXACT_RESUME_VALIDATION",
        )
    materialization = [paths["aligned"], paths["provenance"], paths["task"]]
    present = [artifact_present(path) for path in materialization]
    if all(present):
        for path in materialization:
            stable_identity(path)
        try:
            task = ENGINE.one_tsv_row(paths["task"])
            required = int(task["minimum_free_bytes"])
        except (SystemExit, KeyError, TypeError, ValueError) as error:
            raise CoordinatorError("materialized task lacks its frozen runner storage floor") from error
        if required != int(ENGINE.MINIMUM_FREE_BYTES):
            raise CoordinatorError("materialized task runner storage floor drifted")
        return (
            forecast_with_required(free_bytes, required),
            "EXISTING_MATERIALIZATION_REQUIRES_DEEP_RESUME_VALIDATION",
        )
    state = "PARTIAL_MATERIALIZATION_REQUIRES_AUDITED_REBUILD" if any(present) else "NEW_MATERIALIZATION"
    return forecast_with_required(free_bytes, pair_required_storage_bytes(context)), state


def pending_pair_canonical_reserve_bytes(context: dict[str, Any]) -> int:
    """Count only future canonical bytes that are not already allocated on disk."""
    if artifact_present(cleanup_complete_path(context)) or artifact_present(cleanup_bundle_path(context)):
        return 0
    paths = ENGINE.materialized_paths(context)
    if canonical_state(paths) != "ABSENT":
        return 0
    if artifact_present(paths["staged_ledger"]):
        stable_identity(paths["staged_ledger"])
        return 0  # Canonical publication is an exact hard link to this inode.
    return pending_canonical_planning_reserve_bytes(context)


def storage_forecast(context: dict[str, Any], free_bytes: int | None = None) -> dict[str, int]:
    free = shutil.disk_usage(ROOT).free if free_bytes is None else int(free_bytes)
    required = pair_required_storage_bytes(context)
    return forecast_with_required(free, required)


def require_storage(context: dict[str, Any]) -> dict[str, int]:
    forecast = storage_forecast(context)
    if forecast["additional_bytes_required"]:
        raise CoordinatorError(
            "BLOCKED_BY_COMPUTE insufficient storage for exact pair "
            f"{context['pair_id']}: required={forecast['required_bytes']} "
            f"free={forecast['free_bytes']} additional_required={forecast['additional_bytes_required']}"
        )
    return forecast


def require_runner_storage(context: dict[str, Any]) -> dict[str, int]:
    paths = ENGINE.materialized_paths(context)
    try:
        task = ENGINE.one_tsv_row(paths["task"])
        required = int(task["minimum_free_bytes"])
    except (SystemExit, KeyError, TypeError, ValueError) as error:
        raise CoordinatorError("materialized task lacks its frozen runner storage floor") from error
    free = shutil.disk_usage(ROOT).free
    if required != int(ENGINE.MINIMUM_FREE_BYTES):
        raise CoordinatorError("materialized task runner storage floor drifted")
    if free < required:
        raise CoordinatorError(
            "BLOCKED_BY_COMPUTE insufficient post-materialization runner storage: "
            f"pair={context['pair_id']} required={required} free={free} "
            f"additional_required={required - free}"
        )
    return {"free_bytes": free, "required_bytes": required, "additional_bytes_required": 0}


def admit_ram(
    benchmark: dict[str, str], gate_lock: dict[str, Any], *, physical_memory_bytes: int | None = None,
) -> dict[str, int]:
    try:
        workers = int(benchmark["workers"])
        peak = int(benchmark["peak_process_tree_rss_bytes"])
        reserve = int(gate_lock["resource_envelope"]["ram_safety_reserve_bytes"])
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise CoordinatorError("benchmark lacks numeric RAM admission evidence") from error
    physical = ENGINE.physical_memory_bytes() if physical_memory_bytes is None else int(physical_memory_bytes)
    limit = min(RAM_LIMIT_BYTES, physical)
    if workers != WORKERS:
        raise CoordinatorError("sequential PLACO+ production requires the one-worker benchmark")
    if peak <= 0 or reserve < 0 or peak + reserve > limit:
        raise CoordinatorError(
            "BLOCKED_BY_COMPUTE measured PLACO+ peak plus frozen reserve exceeds the 8-GiB envelope: "
            f"peak={peak} reserve={reserve} limit={limit} additional_required={max(0, peak + reserve - limit)}"
        )
    return {"peak_bytes": peak, "reserve_bytes": reserve, "limit_bytes": limit}


def record_directory(context: dict[str, Any]) -> Path:
    return ROOT / RECORD_ROOT / context["fingerprint"] / context["pair_id"]


def allowed_record_publication_target(base: Path, target: Path) -> bool:
    try:
        relative = target.relative_to(base)
    except ValueError:
        return False
    parts = relative.parts
    if len(parts) == 1:
        return parts[0] in {
            "reproducibility_and_cleanup.plan.json", "cleanup.complete.json",
        }
    if len(parts) != 2:
        return False
    directory, name = parts
    if directory == "ram":
        return name.endswith(".tsv") or name.endswith(".provenance.json")
    if directory == "attempts":
        return name.endswith(".log") or name.endswith(".json")
    if directory == "repairs":
        return name.endswith(".json")
    return False


def repair_record_publication_temporaries(context: dict[str, Any]) -> None:
    """Promote or remove only hash-complete record temps under the pair lease."""
    require_active_worker_lease(context)
    base = record_directory(context)
    if not base.exists():
        return
    safe_repo_path(base, "PLACO+ execution-record root")
    grouped: dict[Path, bytes] = {}
    for temporary in sorted(base.rglob("*")):
        if not temporary.name.startswith(".") or not temporary.name.endswith(".tmp"):
            continue
        if not temporary.is_file() or temporary.is_symlink():
            raise CoordinatorError(f"record publication temporary is not regular: {temporary}")
        hashed = PUBLICATION_TEMP_WITH_HASH_RE.fullmatch(temporary.name)
        legacy = PUBLICATION_TEMP_LEGACY_RE.fullmatch(temporary.name)
        match = hashed or legacy
        if match is None:
            raise CoordinatorError(f"malformed record publication temporary: {temporary}")
        target = temporary.with_name(match.group("target"))
        if not allowed_record_publication_target(base, target):
            raise CoordinatorError(f"record publication temporary targets an unknown artifact: {temporary}")
        if artifact_present(target):
            content = target.read_bytes()
        elif hashed is not None:
            content = temporary.read_bytes()
        else:
            content = b"\x00LEGACY_TEMP_WITHOUT_FINAL_IS_UNVERIFIABLE"
        prior = grouped.setdefault(target, content)
        if prior != content:
            raise CoordinatorError("record publication temporaries disagree for one destination")
    for target, content in grouped.items():
        reconcile_publication_temporaries_for_target(target, content)


def ram_receipt_path(context: dict[str, Any], phase: str) -> Path:
    return record_directory(context) / "ram" / f"{phase}.tsv"


def record_ram(
    context: dict[str, Any], phase: str, *, analysis: str, peak_bytes: int,
    runtime_seconds: float, exit_status: int, output_hash: str,
    rss_measurement_method: str,
) -> dict[str, object]:
    if not re.fullmatch(r"[0-9a-f]{64}", output_hash):
        raise CoordinatorError("RAM receipt output hash must be SHA-256")
    if (
        peak_bytes < 0
        or (exit_status == 0 and peak_bytes == 0)
        or not math.isfinite(runtime_seconds)
        or runtime_seconds <= 0
    ):
        raise CoordinatorError(
            "RAM receipt requires a nonnegative peak (positive on success) and positive runtime"
        )
    if rss_measurement_method not in RSS_MEASUREMENT_METHODS:
        raise CoordinatorError("RAM receipt uses an unknown measurement method")
    if (peak_bytes == 0) != (
        rss_measurement_method == "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE"
    ):
        raise CoordinatorError("zero-RSS failure evidence and its measurement method disagree")
    paths = ENGINE.materialized_paths(context)
    task = ENGINE.one_tsv_row(paths["task"])
    row: dict[str, object] = {
        "analysis": analysis,
        "pair": context["pair_id"],
        "locus": "GENOME_WIDE_VALID_VARIANT_FAMILY",
        "chromosome": "ALL",
        "n_snps": int(task["aligned_input_rows"]),
        "peak_ram_gb": f"{peak_bytes / 1024**3:.12g}",
        "runtime_sec": f"{runtime_seconds:.12g}",
        "exit_status": exit_status,
        "output_hash": output_hash,
    }
    receipt = ram_receipt_path(context, phase)
    content = tsv_bytes(RAM_FIELDS, [row])
    provenance = {
        "schema_version": "sleep-atlas-track-b-placo-ram-measurement.1",
        "analysis_id": context["policy"]["analysis_id"],
        "pair_id": context["pair_id"],
        "phase": phase,
        "run_fingerprint": context["fingerprint"],
        "workers": WORKERS,
        "peak_process_tree_rss_bytes": peak_bytes,
        "runtime_seconds": runtime_seconds,
        "exit_status": exit_status,
        "rss_measurement_method": rss_measurement_method,
        "sample_interval_seconds": 0.05,
        "output_hash": output_hash,
        "task": stable_identity(paths["task"]),
        "terminal_gate_lock_sha256": stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK)["sha256"],
        "coordinator_sha256": coordinator_sha256(),
        "row": row,
        "ram_tsv": {
            "path": str(receipt.relative_to(ROOT)),
            "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
            "fields": RAM_FIELDS,
        },
    }
    publish_bytes_no_replace(
        receipt.with_suffix(".provenance.json"),
        json.dumps(provenance, indent=2, sort_keys=True).encode() + b"\n",
    )
    publish_bytes_no_replace(receipt, content)
    rebuild_ram_aggregate()
    return row


def validate_ram_measurement_provenance(
    provenance: dict[str, Any], receipt: Path, *, context: dict[str, Any] | None = None,
    require_live_task: bool = True,
) -> tuple[dict[str, str], bytes]:
    """Validate both representations of one immutable RAM measurement.

    The exact byte count used for admission lives in JSON while the requested
    public table expresses GiB.  Treating those as unrelated fields would let a
    damaged resume record understate RAM while retaining an apparently valid
    TSV row, so both forms are reconstructed and compared here.
    """
    row = provenance.get("row")
    if not isinstance(row, dict) or set(row) != set(RAM_FIELDS):
        raise CoordinatorError("RAM measurement provenance lacks the exact requested row")
    try:
        peak = int(provenance["peak_process_tree_rss_bytes"])
        runtime = float(provenance["runtime_seconds"])
        exit_status = int(provenance["exit_status"])
        workers = int(provenance["workers"])
        sample_interval = float(provenance["sample_interval_seconds"])
        n_snps = int(row["n_snps"])
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise CoordinatorError("RAM measurement provenance has invalid numeric fields") from error
    method = provenance.get("rss_measurement_method")
    output_hash = provenance.get("output_hash")
    if (
        peak < 0
        or (exit_status == 0 and peak == 0)
        or not math.isfinite(runtime) or runtime <= 0
        or workers != WORKERS
        or not math.isclose(sample_interval, 0.05, rel_tol=0, abs_tol=0)
        or n_snps <= 0
        or method not in RSS_MEASUREMENT_METHODS
        or not isinstance(output_hash, str)
        or re.fullmatch(r"[0-9a-f]{64}", output_hash) is None
    ):
        raise CoordinatorError("RAM measurement provenance has invalid execution evidence")
    if (peak == 0) != (method == "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE"):
        raise CoordinatorError("zero-RSS measurement provenance has the wrong method")
    pair_id = provenance.get("pair_id")
    fingerprint = provenance.get("run_fingerprint")
    phase = provenance.get("phase")
    if (
        pair_id not in PAIR_ORDER
        or not isinstance(fingerprint, str)
        or re.fullmatch(r"[0-9a-f]{64}", fingerprint) is None
        or not isinstance(phase, str) or not phase
        or receipt.name != f"{phase}.tsv"
    ):
        raise CoordinatorError("RAM measurement provenance path identity drifted")
    wanted_row = {
        "analysis": row["analysis"],
        "pair": pair_id,
        "locus": "GENOME_WIDE_VALID_VARIANT_FAMILY",
        "chromosome": "ALL",
        "n_snps": str(n_snps),
        "peak_ram_gb": f"{peak / 1024**3:.12g}",
        "runtime_sec": f"{runtime:.12g}",
        "exit_status": str(exit_status),
        "output_hash": output_hash,
    }
    ordered_row = {field: str(row[field]) for field in RAM_FIELDS}
    if not str(wanted_row["analysis"]).strip() or ordered_row != wanted_row:
        raise CoordinatorError("RAM measurement TSV row disagrees with exact JSON execution evidence")
    content = tsv_bytes(RAM_FIELDS, [wanted_row])
    expected_task_path = str(
        ENGINE.PLACO_WORK_ROOT / fingerprint / pair_id / "task.tsv"
    )
    task_record = provenance.get("task")
    if (
        not isinstance(task_record, dict)
        or task_record.get("path") != expected_task_path
        or not isinstance(task_record.get("bytes"), int)
        or int(task_record["bytes"]) <= 0
        or re.fullmatch(r"[0-9a-f]{64}", str(task_record.get("sha256", ""))) is None
    ):
        raise CoordinatorError("RAM measurement provenance lacks its exact materialized task identity")
    expected = {
        "schema_version": "sleep-atlas-track-b-placo-ram-measurement.1",
        "analysis_id": "track-b-v1.0-pleiotropy",
        "pair_id": pair_id,
        "phase": phase,
        "run_fingerprint": fingerprint,
        "workers": WORKERS,
        "terminal_gate_lock_sha256": stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK)["sha256"],
        "coordinator_sha256": coordinator_sha256(),
    }
    if context is not None:
        expected.update({
            "analysis_id": context["policy"]["analysis_id"],
            "pair_id": context["pair_id"],
            "run_fingerprint": context["fingerprint"],
        })
    if any(provenance.get(key) != value for key, value in expected.items()):
        raise CoordinatorError("RAM measurement provenance identity drifted")
    tsv_record = provenance.get("ram_tsv", {})
    if (
        tsv_record.get("path") != str(receipt.relative_to(ROOT))
        or tsv_record.get("bytes") != len(content)
        or tsv_record.get("sha256") != hashlib.sha256(content).hexdigest()
        or tsv_record.get("fields") != RAM_FIELDS
    ):
        raise CoordinatorError("RAM measurement TSV binding drifted")
    if require_live_task:
        task = ROOT / expected_task_path
        if task_record != stable_identity(task):
            raise CoordinatorError("RAM measurement task identity drifted")
        task_row = ENGINE.one_tsv_row(task)
        if int(task_row.get("aligned_input_rows", -1)) != n_snps:
            raise CoordinatorError("RAM measurement variant count differs from its task")
    return wanted_row, content


def load_ram_measurement(
    context: dict[str, Any], phase: str, *, require_live_task: bool = True,
) -> dict[str, Any] | None:
    receipt = ram_receipt_path(context, phase)
    provenance_path = receipt.with_suffix(".provenance.json")
    receipt_present = receipt.exists() or receipt.is_symlink()
    provenance_present = provenance_path.exists() or provenance_path.is_symlink()
    if receipt_present and not provenance_present:
        raise CoordinatorError(f"RAM receipt lacks immutable provenance: {receipt}")
    if not provenance_present:
        return None
    provenance = read_json(provenance_path)
    ordered_row, content = validate_ram_measurement_provenance(
        provenance, receipt, context=context, require_live_task=require_live_task,
    )
    publish_bytes_no_replace(receipt, content)
    fields, rows = ENGINE.read_tsv(receipt)
    if fields != RAM_FIELDS or rows != [ordered_row]:
        raise CoordinatorError("RAM measurement TSV differs from immutable provenance")
    return provenance


def full_scan_rejection_phases(context: dict[str, Any]) -> list[str]:
    ram_root = record_directory(context) / "ram"
    if not ram_root.exists():
        return []
    safe_repo_path(ram_root, "PLACO+ RAM evidence root")
    phases: set[str] = set()
    for path in ram_root.iterdir():
        name = path.name
        if name.endswith(".provenance.json"):
            phase = name.removesuffix(".provenance.json")
        elif name.endswith(".tsv"):
            phase = name.removesuffix(".tsv")
        else:
            continue
        if phase.startswith(FULL_SCAN_REJECTION_PREFIX):
            phases.add(phase)
    return sorted(phases)


def purge_unmeasured_full_scan_checkpoints(context: dict[str, Any]) -> None:
    """Never bootstrap production from shards created without immutable RSS evidence."""
    successful = ram_receipt_path(context, "full_scan").with_suffix(".provenance.json")
    if artifact_present(successful):
        return
    paths = ENGINE.materialized_paths(context)
    checkpoints = sorted(paths["checkpoint_dir"].glob("shard_*.rds*"))
    if not checkpoints:
        return
    audit_and_unlink_regenerable(
        context, "unmeasured_full_scan_checkpoint_lineage", checkpoints,
        "NO_IMMUTABLE_PROCESS_TREE_RSS_EVIDENCE_FOR_COMPLETED_SHARDS;"
        "RECOMPUTE_IDENTICAL_SINGLE_VARIANT_SHARDS_WITH_FROZEN_GLOBAL_NUISANCE",
    )


def require_full_scan_ram_evidence(
    context: dict[str, Any], gate_lock: dict[str, Any], output_hash: str,
    *, require_live_task: bool = True,
) -> dict[str, Any]:
    """Require an admitted successful full-scan receipt bound to the exact ledger."""
    rejected = full_scan_rejection_phases(context)
    for phase in rejected:
        evidence = load_ram_measurement(
            context, phase, require_live_task=require_live_task,
        )
        if evidence is None:
            raise CoordinatorError("full-scan RAM rejection evidence is incomplete")
        if evidence.get("row", {}).get("analysis") == FULL_SCAN_ANALYSIS:
            raise CoordinatorError("full-scan rejection evidence falsely claims a successful analysis")
    if rejected:
        raise CoordinatorError(
            "BLOCKED_BY_COMPUTE immutable full-scan RAM-envelope rejection evidence exists: "
            + ",".join(rejected)
        )
    evidence = load_ram_measurement(
        context, "full_scan", require_live_task=require_live_task,
    )
    if evidence is None:
        raise CoordinatorError("canonical PLACO+ result lacks successful full-scan RAM evidence")
    if (
        evidence.get("phase") != "full_scan"
        or evidence.get("row", {}).get("analysis") != FULL_SCAN_ANALYSIS
        or int(evidence.get("exit_status", -1)) != 0
        or evidence.get("output_hash") != output_hash
    ):
        raise CoordinatorError("full-scan RAM evidence is not bound to the exact canonical ledger")
    admit_ram(
        {
            "workers": str(evidence["workers"]),
            "peak_process_tree_rss_bytes": str(evidence["peak_process_tree_rss_bytes"]),
        },
        gate_lock,
        physical_memory_bytes=RAM_LIMIT_BYTES,
    )
    return evidence


def atomic_replace_derived(path: Path, content: bytes) -> None:
    safe_repo_path(path.parent, "derived publication parent")
    path.parent.mkdir(parents=True, exist_ok=True)
    safe_repo_path(path.parent, "derived publication parent")
    reconcile_publication_temporaries_for_target(path, content, derived=True)
    if path.is_symlink():
        raise CoordinatorError(f"derived publication may not replace a symlink: {path}")
    if path.exists():
        observed = os.lstat(path)
        if not stat.S_ISREG(observed.st_mode) or observed.st_nlink != 1:
            raise CoordinatorError(f"derived publication is not a private regular file: {path}")
        if path.read_bytes() == content:
            return
    content_hash = hashlib.sha256(content).hexdigest()
    temporary = path.with_name(f".{path.name}.{content_hash}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def rebuild_ram_aggregate() -> None:
    rows: list[dict[str, str]] = []
    record_root = ROOT / RECORD_ROOT
    fingerprints_by_pair: dict[str, str] = {}
    receipt_paths: list[Path] = []
    provenance_paths: list[Path] = []
    if record_root.exists():
        safe_repo_path(record_root, "PLACO RAM record root")
        for fingerprint_root in sorted(record_root.iterdir()):
            observed_root = os.lstat(fingerprint_root)
            if (
                not stat.S_ISDIR(observed_root.st_mode)
                or re.fullmatch(r"[0-9a-f]{64}", fingerprint_root.name) is None
            ):
                raise CoordinatorError(f"invalid PLACO RAM fingerprint namespace: {fingerprint_root}")
            pair_directories = sorted(path for path in fingerprint_root.iterdir() if path.name in PAIR_ORDER)
            unknown = sorted(path.name for path in fingerprint_root.iterdir() if path.name not in PAIR_ORDER)
            if unknown or len(pair_directories) != 1:
                raise CoordinatorError(
                    f"PLACO RAM fingerprint namespace must contain exactly one known pair: {fingerprint_root}"
                )
            pair = pair_directories[0].name
            pair_stat = os.lstat(pair_directories[0])
            if not stat.S_ISDIR(pair_stat.st_mode) or pair in fingerprints_by_pair:
                raise CoordinatorError(f"duplicate or unsafe PLACO RAM pair namespace: {pair}")
            fingerprints_by_pair[pair] = fingerprint_root.name
            ram_root = pair_directories[0] / "ram"
            if not ram_root.exists():
                continue
            ram_stat = os.lstat(ram_root)
            if not stat.S_ISDIR(ram_stat.st_mode):
                raise CoordinatorError(f"PLACO RAM namespace is not a directory: {ram_root}")
            entries = sorted(ram_root.iterdir())
            unknown_ram = [
                path for path in entries
                if not (path.name.endswith(".tsv") or path.name.endswith(".provenance.json"))
            ]
            if unknown_ram:
                raise CoordinatorError(f"unknown PLACO RAM artifact: {unknown_ram[0]}")
            for provenance_path in [path for path in entries if path.name.endswith(".provenance.json")]:
                if not stat.S_ISREG(os.lstat(provenance_path).st_mode):
                    raise CoordinatorError(f"unsafe PLACO RAM provenance: {provenance_path}")
                phase = provenance_path.name.removesuffix(".provenance.json")
                receipt = provenance_path.with_name(f"{phase}.tsv")
                provenance = read_json(provenance_path)
                row, content = validate_ram_measurement_provenance(
                    provenance, receipt, require_live_task=False,
                )
                if provenance.get("pair_id") != pair or provenance.get("run_fingerprint") != fingerprint_root.name:
                    raise CoordinatorError("PLACO RAM record is stored under the wrong pair fingerprint")
                publish_bytes_no_replace(receipt, content)
                fields, observed = ENGINE.read_tsv(receipt)
                if fields != RAM_FIELDS or observed != [row]:
                    raise CoordinatorError(f"invalid immutable RAM receipt: {receipt}")
                rows.append(row)
                receipt_paths.append(receipt)
                provenance_paths.append(provenance_path)
            orphan_receipts = [
                path for path in ram_root.glob("*.tsv")
                if not path.with_suffix(".provenance.json").exists()
            ]
            if orphan_receipts:
                raise CoordinatorError(f"PLACO RAM receipt lacks provenance: {orphan_receipts[0]}")
    order = {pair: index for index, pair in enumerate(PAIR_ORDER)}
    combined = sorted(
        zip(rows, receipt_paths, provenance_paths, strict=True),
        key=lambda item: (order[item[0]["pair"]], item[1].name),
    )
    rows = [item[0] for item in combined]
    receipt_paths = [item[1] for item in combined]
    provenance_paths = [item[2] for item in combined]
    content = tsv_bytes(RAM_FIELDS, rows)
    aggregate = ROOT / RAM_AGGREGATE
    atomic_replace_derived(aggregate, content)
    family_fingerprint = hashlib.sha256(canonical_json_bytes({
        "pair_order": list(PAIR_ORDER),
        "pair_fingerprints": fingerprints_by_pair,
    })).hexdigest()
    provenance = {
        "schema_version": "sleep-atlas-track-b-placo-ram-benchmark-derived.1",
        "analysis_id": "track-b-v1.0-pleiotropy",
        "row_count": len(rows),
        "pair_order": list(PAIR_ORDER),
        "execution_fingerprints_by_pair": fingerprints_by_pair,
        "execution_family_fingerprint": family_fingerprint,
        "sequential_contract": stable_identity(ROOT / EXECUTION_CONTRACT),
        "source_receipts": [stable_identity(path) for path in receipt_paths],
        "source_measurement_provenance": [stable_identity(path) for path in provenance_paths],
        "aggregate": {
            "path": str(RAM_AGGREGATE),
            "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
            "schema": RAM_FIELDS,
        },
        "note": "Additive PLACO measurements; federation into the source-fingerprint-bound top-level table is separate.",
    }
    atomic_replace_derived(
        ROOT / RAM_AGGREGATE_PROVENANCE,
        json.dumps(provenance, indent=2, sort_keys=True).encode() + b"\n",
    )


def benchmark_bundle_hash(paths: dict[str, Path], benchmark: dict[str, str]) -> str:
    benchmark_identity = stable_identity(paths["benchmark"])
    nuisance_identity = stable_identity(paths["checkpoint_dir"] / "nuisance.rds")
    payload = {
        "benchmark": {
            key: benchmark_identity[key] for key in ("path", "bytes", "sha256")
        },
        "nuisance": {
            key: nuisance_identity[key] for key in ("path", "bytes", "sha256")
        },
        "input_sha256": benchmark["input_sha256"],
        "workers": benchmark["workers"],
    }
    return hashlib.sha256(canonical_json_bytes(payload)).hexdigest()


def validate_admit_fresh_benchmark(
    context: dict[str, Any], gate_lock: dict[str, Any],
    *, require_ram_receipt: bool,
) -> dict[str, str]:
    try:
        row = BRIDGE.validate_benchmark(context, gate_lock, required_workers=WORKERS)
    except (BRIDGE.BridgeError, SystemExit) as error:
        raise CoordinatorError(str(error)) from error
    if row.get("global_nuisance_checkpoint_reused") != "FALSE":
        raise CoordinatorError(
            "matched PLACO+ benchmark must freshly measure global nuisance estimation"
        )
    admit_ram(row, gate_lock)
    if require_ram_receipt:
        evidence = load_ram_measurement(context, "matched_worker_benchmark")
        paths = ENGINE.materialized_paths(context)
        output_hash = benchmark_bundle_hash(paths, row)
        if (
            evidence is None
            or evidence.get("row", {}).get("analysis") != MATCHED_BENCHMARK_ANALYSIS
            or int(evidence.get("exit_status", -1)) != 0
            or int(evidence.get("peak_process_tree_rss_bytes", -1))
            != int(row["peak_process_tree_rss_bytes"])
            or not math.isclose(
                float(evidence.get("runtime_seconds", -1)),
                float(row["wrapper_wall_seconds"]), rel_tol=0, abs_tol=0,
            )
            or evidence.get("output_hash") != output_hash
        ):
            raise CoordinatorError("matched benchmark lacks its exact admitted RAM receipt")
    return row


def audit_and_unlink_regenerable(
    context: dict[str, Any], label: str, paths: list[Path], reason: str,
    *, allow_empty: bool = False,
) -> None:
    base = ENGINE.materialized_paths(context)["base"]
    identities: list[dict[str, object]] = []
    for path in paths:
        if not (path.exists() or path.is_symlink()):
            continue
        try:
            relative = path.relative_to(base).as_posix()
        except ValueError as error:
            raise CoordinatorError("repair candidate escapes the pair work root") from error
        if not allowed_cleanup_relative(relative):
            raise CoordinatorError(f"repair candidate is not explicitly regenerable: {relative}")
        identities.append(stable_identity(path, allow_empty=allow_empty))
    if not identities:
        return
    receipt = {
        "schema_version": "sleep-atlas-track-b-placo-regenerable-repair.1",
        "pair_id": context["pair_id"],
        "run_fingerprint": context["fingerprint"],
        "label": label,
        "reason": reason,
        "scientific_result": False,
        "removed": identities,
        "source_or_canonical_deletion": False,
        "coordinator_sha256": coordinator_sha256(),
    }
    digest = hashlib.sha256(canonical_json_bytes(receipt)).hexdigest()
    audit = record_directory(context) / "repairs" / f"{label}.{digest}.json"
    publish_bytes_no_replace(audit, json.dumps(receipt, indent=2, sort_keys=True).encode() + b"\n")
    for record in identities:
        path = ROOT / str(record["path"])
        safe_repo_path(path, "regenerable repair candidate")
        current = os.lstat(path)
        if (
            not stat.S_ISREG(current.st_mode)
            or (current.st_dev, current.st_ino, current.st_size)
            != (record["device"], record["inode"], record["bytes"])
        ):
            raise CoordinatorError("regenerable repair candidate changed before unlink")
        path.unlink()
        fsync_directory(path.parent)


def repair_crash_temporaries(context: dict[str, Any]) -> None:
    """Audit and remove only documented, non-scientific crash temporaries.

    The pair-stable worker lease is held by the production caller and inherited
    by its R child, so a matching PID-named file cannot belong to a live prior
    PLACO+ execution.  Empty files are valid crash evidence and are
    hashed/audited rather than treated as scientific artifacts.
    """
    base = ENGINE.materialized_paths(context)["base"]
    if not base.exists():
        return
    safe_repo_path(base, "PLACO pair work root")
    temporaries: list[Path] = []
    for path in sorted(base.rglob("*")):
        observed = os.lstat(path)
        if stat.S_ISDIR(observed.st_mode):
            continue
        relative = path.relative_to(base).as_posix()
        if KNOWN_TEMP_RE.fullmatch(relative) is None:
            continue
        if not stat.S_ISREG(observed.st_mode):
            raise CoordinatorError(f"PLACO crash temporary is not a regular file: {path}")
        temporaries.append(path)
    audit_and_unlink_regenerable(
        context, "documented_crash_temporaries", temporaries,
        "STALE_PID_SCOPED_TEMPORARIES_AFTER_PROCESS_DEATH;NEVER_SCIENTIFIC_RESULTS",
        allow_empty=True,
    )


def repair_partial_materialization(context: dict[str, Any]) -> None:
    paths = ENGINE.materialized_paths(context)
    family = [paths["aligned"], paths["provenance"], paths["task"]]
    present = [path for path in family if path.exists() or path.is_symlink()]
    if not present or len(present) == len(family):
        return
    downstream = [
        paths["benchmark_raw"], paths["benchmark"], paths["staged_ledger"],
        paths["run_summary"], paths["staged_provenance"],
        paths["canonical_ledger"], paths["canonical_provenance"],
    ]
    if any(path.exists() or path.is_symlink() for path in downstream):
        raise CoordinatorError("partial materialization coexists with downstream or canonical evidence")
    if paths["checkpoint_dir"].exists() and any(paths["checkpoint_dir"].iterdir()):
        raise CoordinatorError("partial materialization coexists with PLACO+ checkpoints")
    audit_and_unlink_regenerable(
        context, "partial_materialization", present,
        "PROCESS_DEATH_DURING_NO_REPLACE_MATERIALIZATION_FAMILY;EXACT_DENSE_ALIGNMENT_IS_REGENERABLE",
    )


def repair_partial_execution_families(context: dict[str, Any]) -> None:
    paths = ENGINE.materialized_paths(context)
    if canonical_state(paths) != "ABSENT":
        return
    checkpoint = paths["checkpoint_dir"]
    nuisance = [checkpoint / "nuisance.rds", checkpoint / "nuisance.sha256.tsv"]
    nuisance_present = [path for path in nuisance if path.exists() or path.is_symlink()]
    if len(nuisance_present) == 1:
        downstream = [paths["benchmark_raw"], paths["benchmark"], paths["staged_ledger"], paths["run_summary"]]
        downstream.extend(checkpoint.glob("shard_*.rds*"))
        if any(path.exists() or path.is_symlink() for path in downstream):
            raise CoordinatorError("partial nuisance family coexists with downstream PLACO+ evidence")
        audit_and_unlink_regenerable(
            context, "partial_nuisance", nuisance_present,
            "PROCESS_DEATH_DURING_GLOBAL_NUISANCE_CHECKPOINT_PUBLICATION;RECOMPUTE_ON_ALL_VALID_VARIANTS",
        )
    if checkpoint.exists():
        shard_numbers: set[str] = set()
        for path in checkpoint.glob("shard_*.rds*"):
            match = re.fullmatch(r"shard_([0-9]{6})\.rds(?:\.sha256\.tsv)?", path.name)
            if match is None:
                raise CoordinatorError(f"unknown shard checkpoint artifact: {path.name}")
            shard_numbers.add(match.group(1))
        for number in sorted(shard_numbers):
            family = [
                checkpoint / f"shard_{number}.rds",
                checkpoint / f"shard_{number}.rds.sha256.tsv",
            ]
            present = [path for path in family if path.exists() or path.is_symlink()]
            if len(present) == 1:
                if paths["staged_ledger"].exists() or paths["run_summary"].exists():
                    raise CoordinatorError("partial shard family coexists with staged PLACO+ evidence")
                audit_and_unlink_regenerable(
                    context, f"partial_shard_{number}", present,
                    "PROCESS_DEATH_DURING_SINGLE_VARIANT_SHARD_PUBLICATION;SHARD_IS_REGENERABLE_WITH_FROZEN_NUISANCE",
                )
    stage = [paths["staged_ledger"], paths["run_summary"]]
    stage_present = [path for path in stage if path.exists() or path.is_symlink()]
    if len(stage_present) == 1:
        if paths["staged_provenance"].exists() or paths["staged_provenance"].is_symlink():
            raise CoordinatorError("partial staged result coexists with publication provenance")
        audit_and_unlink_regenerable(
            context, "partial_stage", stage_present,
            "PROCESS_DEATH_DURING_STAGED_LEDGER_FAMILY_PUBLICATION;REBUILD_FROM_COMPLETE_WORKER_BOUND_SHARDS",
        )


def ensure_materialized(context: dict[str, Any]) -> dict[str, Any]:
    paths = ENGINE.materialized_paths(context)
    family = [paths["aligned"], paths["provenance"], paths["task"]]
    present = [path.exists() or path.is_symlink() for path in family]
    try:
        if all(present):
            return ENGINE.verify_materialized(context)
        if any(present):
            repair_partial_materialization(context)
        require_storage(context)
        ENGINE.resource_preflight(context)
        return ENGINE.materialize_pair(context)
    except SystemExit as error:
        raise CoordinatorError(str(error)) from error


def repair_partial_benchmark(context: dict[str, Any]) -> None:
    paths = ENGINE.materialized_paths(context)
    raw, final = paths["benchmark_raw"], paths["benchmark"]
    successful_measurement = ram_receipt_path(
        context, "matched_worker_benchmark",
    ).with_suffix(".provenance.json")
    if artifact_present(final) and artifact_present(successful_measurement):
        return
    if artifact_present(successful_measurement):
        raise CoordinatorError(
            "matched-benchmark RAM provenance exists without its immutable benchmark"
        )
    forbidden = [paths["staged_ledger"], paths["run_summary"], paths["staged_provenance"]]
    forbidden.extend(paths["checkpoint_dir"].glob("shard_*.rds*"))
    if any(path.exists() or path.is_symlink() for path in forbidden):
        raise CoordinatorError("benchmark wrapper is missing after full-scan evidence appeared")
    candidates = [
        raw,
        final,
        paths["checkpoint_dir"] / "nuisance.rds",
        paths["checkpoint_dir"] / "nuisance.sha256.tsv",
    ]
    present = [path for path in candidates if artifact_present(path)]
    audit_and_unlink_regenerable(
        context, "unadmitted_benchmark_lineage", present,
        "MATCHED_BENCHMARK_WRAPPER_NOT_IMMUTABLY_ATTESTED;"
        "RECOMPUTE_GLOBAL_NUISANCE_ON_ALL_VALID_VARIANTS_UNDER_THE_LIVE_RAM_GUARD",
    )


def failure_event_log_content(process_log: bytes, metadata: dict[str, object]) -> bytes:
    trailer = "".join(f"{key}={metadata[key]}\n" for key in sorted(metadata)).encode()
    return process_log + (b"\n" if process_log and not process_log.endswith(b"\n") else b"") + trailer


def failure_group_termination_status(
    termination_reason: str | None, *, process_started: bool,
) -> str:
    if not process_started:
        return "NOT_STARTED"
    if termination_reason is None:
        return "NOT_REQUIRED_PROCESS_ALREADY_EXITED"
    if termination_reason.endswith("_PROCESS_GROUP_TERMINATION_UNCONFIRMED"):
        return "UNCONFIRMED"
    return "CONFIRMED"


def publish_failure_event(
    context: dict[str, Any], *, failure_class: str, phase_prefix: str,
    analysis: str, attempt_log_label: str, process_log: bytes,
    peak_bytes: int, runtime_seconds: float, public_exit_status: int,
    owned_child_exit_status: int | None, rss_measurement_method: str,
    output_hash: str | None, termination_reason: str | None,
    process_group_termination_status: str, terminal_compute_rejection: bool,
    details: dict[str, object] | None = None,
) -> dict[str, Any]:
    """Commit one self-contained failure event before any derived log/RAM row."""
    require_active_worker_lease(context)
    if (
        not isinstance(failure_class, str)
        or re.fullmatch(r"[A-Z0-9_]+", failure_class) is None
        or not isinstance(phase_prefix, str)
        or re.fullmatch(r"[a-z0-9_]+_", phase_prefix) is None
        or not isinstance(analysis, str) or not analysis.strip()
        or not isinstance(attempt_log_label, str)
        or re.fullmatch(r"[A-Za-z0-9_]+", attempt_log_label) is None
        or not isinstance(process_log, bytes)
        or not isinstance(public_exit_status, int)
        or isinstance(public_exit_status, bool)
        or public_exit_status == 0
        or not isinstance(peak_bytes, int) or isinstance(peak_bytes, bool)
        or peak_bytes < 0
        or not isinstance(runtime_seconds, (int, float))
        or isinstance(runtime_seconds, bool)
        or not math.isfinite(runtime_seconds)
        or runtime_seconds <= 0
        or not isinstance(terminal_compute_rejection, bool)
        or not isinstance(details or {}, dict)
        or rss_measurement_method not in RSS_MEASUREMENT_METHODS
        or process_group_termination_status not in PROCESS_GROUP_TERMINATION_STATUSES
        or (
            output_hash is not None
            and (
                not isinstance(output_hash, str)
                or re.fullmatch(r"[0-9a-f]{64}", output_hash) is None
            )
        )
        or (terminal_compute_rejection and public_exit_status != 78)
        or (
            terminal_compute_rejection
            and not str(termination_reason or "").startswith("RAM_LIMIT_EXCEEDED")
        )
    ):
        raise CoordinatorError("failure event lacks exact non-success execution evidence")
    if (peak_bytes == 0) != (
        rss_measurement_method == "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE"
    ):
        raise CoordinatorError("failure event zero-RSS evidence is inconsistent")
    paths = ENGINE.materialized_paths(context)
    task = ENGINE.one_tsv_row(paths["task"])
    metadata: dict[str, object] = {
        "failure_class": failure_class,
        "live_ram_guard_triggered": str(bool(
            termination_reason and termination_reason.startswith("RAM_LIMIT_EXCEEDED")
        )).upper(),
        "owned_child_exit_status": (
            owned_child_exit_status if owned_child_exit_status is not None else "UNKNOWN"
        ),
        "peak_process_tree_rss_bytes": peak_bytes,
        "process_group_termination_status": process_group_termination_status,
        "public_exit_status": public_exit_status,
        "rss_measurement_method": rss_measurement_method,
        "runtime_seconds": f"{runtime_seconds:.12g}",
        "terminal_compute_rejection": str(terminal_compute_rejection).upper(),
        "termination_reason": termination_reason or "NONE",
    }
    log = failure_event_log_content(process_log, metadata)
    log_hash = hashlib.sha256(log).hexdigest()
    bound_output_hash = log_hash if output_hash is None else output_hash
    core: dict[str, Any] = {
        "analysis": analysis,
        "analysis_id": context["policy"]["analysis_id"],
        "attempt_log": {
            "base64": base64.b64encode(log).decode("ascii"),
            "bytes": len(log),
            "label": attempt_log_label,
            "sha256": log_hash,
        },
        "coordinator_sha256": coordinator_sha256(),
        "details": details or {},
        "failure_class": failure_class,
        "input_sha256": task["aligned_input_sha256"],
        "n_snps": int(task["aligned_input_rows"]),
        "output_hash": bound_output_hash,
        "owned_child_exit_status": owned_child_exit_status,
        "pair_id": context["pair_id"],
        "peak_process_tree_rss_bytes": peak_bytes,
        "phase_prefix": phase_prefix,
        "process_group_termination_status": process_group_termination_status,
        "public_exit_status": public_exit_status,
        "rss_measurement_method": rss_measurement_method,
        "run_fingerprint": context["fingerprint"],
        "runtime_seconds": runtime_seconds,
        "sequential_contract": stable_identity(ROOT / EXECUTION_CONTRACT),
        "task": stable_identity(paths["task"]),
        "terminal_compute_rejection": terminal_compute_rejection,
        "terminal_gate_lock": stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK),
        "termination_reason": termination_reason,
    }
    event_id = hashlib.sha256(canonical_json_bytes(core)).hexdigest()
    event = {
        "schema_version": FAILURE_EVENT_SCHEMA,
        "event_id": event_id,
        "event": core,
    }
    path = record_directory(context) / "attempts" / f"failure_event.{event_id}.json"
    publish_bytes_no_replace(
        path, json.dumps(event, indent=2, sort_keys=True).encode() + b"\n",
    )
    reconcile_failure_event(context, path)
    return event


def validate_failure_event(
    context: dict[str, Any], path: Path,
) -> tuple[dict[str, Any], str, Path, bytes]:
    match = FAILURE_EVENT_RE.fullmatch(path.name)
    if match is None or path.parent != record_directory(context) / "attempts":
        raise CoordinatorError(f"malformed PLACO+ failure-event path: {path}")
    payload = read_json(path)
    if set(payload) != {"schema_version", "event_id", "event"}:
        raise CoordinatorError("failure event has an unexpected schema")
    core = payload.get("event")
    if not isinstance(core, dict) or set(core) != FAILURE_EVENT_CORE_FIELDS:
        raise CoordinatorError("failure event lacks its self-contained evidence")
    event_id = hashlib.sha256(canonical_json_bytes(core)).hexdigest()
    if (
        payload.get("schema_version") != FAILURE_EVENT_SCHEMA
        or payload.get("event_id") != event_id
        or match.group(1) != event_id
        or core.get("analysis_id") != context["policy"]["analysis_id"]
        or core.get("pair_id") != context["pair_id"]
        or core.get("run_fingerprint") != context["fingerprint"]
        or core.get("coordinator_sha256") != coordinator_sha256()
        or core.get("public_exit_status") == 0
        or not isinstance(core.get("terminal_compute_rejection"), bool)
    ):
        raise CoordinatorError("failure event identity or failure status drifted")
    if core.get("sequential_contract") != stable_identity(ROOT / EXECUTION_CONTRACT):
        raise CoordinatorError("failure event sequential-contract identity drifted")
    if core.get("terminal_gate_lock") != stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK):
        raise CoordinatorError("failure event terminal-gate identity drifted")
    try:
        peak = int(core["peak_process_tree_rss_bytes"])
        runtime = float(core["runtime_seconds"])
        status = int(core["public_exit_status"])
        n_snps = int(core["n_snps"])
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise CoordinatorError("failure event has malformed numeric evidence") from error
    method = core.get("rss_measurement_method")
    group_status = core.get("process_group_termination_status")
    if (
        peak < 0 or status == 0 or n_snps <= 0
        or not math.isfinite(runtime) or runtime <= 0
        or method not in RSS_MEASUREMENT_METHODS
        or group_status not in PROCESS_GROUP_TERMINATION_STATUSES
        or (peak == 0) != (method == "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE")
        or not isinstance(core.get("analysis"), str) or not core["analysis"].strip()
        or re.fullmatch(r"[0-9a-f]{64}", str(core.get("output_hash", ""))) is None
        or not isinstance(core.get("phase_prefix"), str) or not core["phase_prefix"]
        or re.fullmatch(r"[a-z0-9_]+_", core["phase_prefix"]) is None
        or not isinstance(core.get("failure_class"), str)
        or re.fullmatch(r"[A-Z0-9_]+", core["failure_class"]) is None
        or not isinstance(core.get("details"), dict)
        or (
            core["terminal_compute_rejection"]
            and not str(core["phase_prefix"]).endswith("ram_rejected_")
        )
        or (core["terminal_compute_rejection"] and status != 78)
        or (
            core["terminal_compute_rejection"]
            and not str(core.get("termination_reason") or "").startswith(
                "RAM_LIMIT_EXCEEDED"
            )
        )
        or (
            not core["terminal_compute_rejection"]
            and str(core["phase_prefix"]).endswith("ram_rejected_")
        )
    ):
        raise CoordinatorError("failure event has invalid execution evidence")
    log_record = core.get("attempt_log")
    if not isinstance(log_record, dict) or set(log_record) != {"base64", "bytes", "label", "sha256"}:
        raise CoordinatorError("failure event lacks its exact process log")
    try:
        log = base64.b64decode(log_record["base64"], validate=True)
    except (TypeError, ValueError) as error:
        raise CoordinatorError("failure event process log is not valid base64") from error
    label = str(log_record.get("label", ""))
    if (
        not label or re.fullmatch(r"[A-Za-z0-9_]+", label) is None
        or len(log) != log_record.get("bytes")
        or hashlib.sha256(log).hexdigest() != log_record.get("sha256")
    ):
        raise CoordinatorError("failure event process log binding drifted")
    phase = str(core["phase_prefix"]) + event_id[:16]
    log_path = record_directory(context) / "attempts" / f"{label}.{log_record['sha256']}.log"
    task_record = core.get("task")
    expected_task_path = str(
        ENGINE.PLACO_WORK_ROOT / context["fingerprint"] / context["pair_id"] / "task.tsv"
    )
    try:
        task_bytes = int(task_record.get("bytes", 0)) if isinstance(task_record, dict) else 0
    except (TypeError, ValueError, OverflowError) as error:
        raise CoordinatorError("failure event task/input identity drifted") from error
    if (
        not isinstance(task_record, dict)
        or set(task_record) != {"path", "device", "inode", "bytes", "sha256"}
        or task_record.get("path") != expected_task_path
        or task_bytes <= 0
        or re.fullmatch(r"[0-9a-f]{64}", str(task_record.get("sha256", ""))) is None
        or re.fullmatch(r"[0-9a-f]{64}", str(core.get("input_sha256", ""))) is None
    ):
        raise CoordinatorError("failure event task/input identity drifted")
    return core, phase, log_path, log


def reconcile_failure_event(context: dict[str, Any], path: Path) -> None:
    core, phase, log_path, log = validate_failure_event(context, path)
    publish_bytes_no_replace(log_path, log)
    provenance_path = ram_receipt_path(context, phase).with_suffix(".provenance.json")
    task_path = ENGINE.materialized_paths(context)["task"]
    provenance_preexisted = artifact_present(provenance_path)
    if not provenance_preexisted:
        if not artifact_present(task_path) or not identity_matches(task_path, core["task"]):
            raise CoordinatorError("unreconciled failure event lost its exact live task")
        task = ENGINE.one_tsv_row(task_path)
        if (
            int(task.get("aligned_input_rows", -1)) != int(core["n_snps"])
            or task.get("aligned_input_sha256") != core["input_sha256"]
        ):
            raise CoordinatorError("failure event input family differs from its exact task")
        record_ram(
            context, phase, analysis=str(core["analysis"]),
            peak_bytes=int(core["peak_process_tree_rss_bytes"]),
            runtime_seconds=float(core["runtime_seconds"]),
            exit_status=int(core["public_exit_status"]),
            output_hash=str(core["output_hash"]),
            rss_measurement_method=str(core["rss_measurement_method"]),
        )
    evidence = load_ram_measurement(
        context, phase, require_live_task=artifact_present(task_path),
    )
    if (
        evidence is None
        or evidence.get("row", {}).get("analysis") != core["analysis"]
        or int(evidence.get("peak_process_tree_rss_bytes", -1))
        != int(core["peak_process_tree_rss_bytes"])
        or not math.isclose(
            float(evidence.get("runtime_seconds", -1)),
            float(core["runtime_seconds"]), rel_tol=0, abs_tol=0,
        )
        or int(evidence.get("exit_status", 0)) != int(core["public_exit_status"])
        or evidence.get("rss_measurement_method") != core["rss_measurement_method"]
        or evidence.get("output_hash") != core["output_hash"]
    ):
        raise CoordinatorError("derived RAM evidence differs from its durable failure event")
    # record_ram() rebuilds this derived table on the first reconciliation.  A
    # coordinator SIGKILL after provenance/TSV publication but before that
    # rebuild leaves a pre-existing provenance file on resume, so explicitly
    # repair the aggregate in that exact crash state as well.
    if provenance_preexisted:
        rebuild_ram_aggregate()


def reconcile_failure_events(context: dict[str, Any]) -> None:
    require_active_worker_lease(context)
    attempts = record_directory(context) / "attempts"
    if not attempts.exists():
        return
    safe_repo_path(attempts, "PLACO+ failure-event namespace")
    for path in sorted(attempts.iterdir()):
        if not path.name.startswith("failure_event."):
            continue
        if FAILURE_EVENT_RE.fullmatch(path.name) is None:
            raise CoordinatorError(f"malformed PLACO+ failure event: {path}")
        reconcile_failure_event(context, path)


def benchmark_failure_measurement(
    context: dict[str, Any], error: BaseException, *, peak_bytes: int,
    runtime_seconds: float, termination_reason: str | None,
    guard_limit_bytes: int, process_started: bool, exit_status: int | None,
    process_log: bytes = b"",
) -> None:
    if peak_bytes <= 0:
        if process_started:
            observed = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            peak_bytes = observed if sys.platform == "darwin" else observed * 1024
        method = (
            "CHILD_RUSAGE_SESSION_MAXRSS_FALLBACK_PS_UNAVAILABLE"
            if peak_bytes > 0
            else "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE"
        )
    else:
        method = "AGGREGATE_PROCESS_TREE_PS_SAMPLED"
    wrapper_log = (
        f"{type(error).__name__}: {error}\n"
        f"termination_reason={termination_reason or 'NONE'}\n"
        f"owned_child_exit_status={exit_status if exit_status is not None else 'UNKNOWN'}\n"
        f"live_ram_guard_limit_bytes={guard_limit_bytes}\n"
        f"live_ram_guard_triggered={str(bool(termination_reason and termination_reason.startswith('RAM_LIMIT_EXCEEDED'))).upper()}\n"
    ).encode()
    log = process_log + (b"\n" if process_log and not process_log.endswith(b"\n") else b"") + wrapper_log
    measured_runtime = max(runtime_seconds, sys.float_info.epsilon)
    ram_rejected = bool(
        termination_reason and termination_reason.startswith("RAM_LIMIT_EXCEEDED")
    )
    public_exit_status = 78 if ram_rejected else (
        int(exit_status)
        if exit_status is not None and int(exit_status) != 0
        else 70
    )
    phase_prefix = (
        "matched_worker_benchmark_ram_rejected_"
        if ram_rejected else "matched_worker_benchmark_failed_"
    )
    analysis = (
        "PLACO_PLUS_MATCHED_WORKER_BENCHMARK_RAM_ENVELOPE_REJECTED"
        if ram_rejected else "PLACO_PLUS_MATCHED_WORKER_BENCHMARK_FAILED_ATTEMPT"
    )
    publish_failure_event(
        context,
        failure_class="MATCHED_BENCHMARK_PROCESS_FAILURE",
        phase_prefix=phase_prefix,
        analysis=analysis,
        peak_bytes=peak_bytes,
        runtime_seconds=measured_runtime,
        public_exit_status=public_exit_status,
        owned_child_exit_status=exit_status,
        rss_measurement_method=method,
        output_hash=None,
        termination_reason=termination_reason,
        process_group_termination_status=failure_group_termination_status(
            termination_reason, process_started=process_started,
        ),
        terminal_compute_rejection=ram_rejected,
        attempt_log_label="matched_worker_benchmark_failed",
        process_log=log,
        details={"live_ram_guard_limit_bytes": guard_limit_bytes},
    )


def run_benchmark_with_live_guard(
    context: dict[str, Any], gate_lock: dict[str, Any],
) -> dict[str, str]:
    """Run the frozen benchmark while observing and enforcing its one-worker cap."""
    maximum = ram_worker_limit_bytes(gate_lock)
    original_sampler = ENGINE.process_tree_rss_bytes
    original_subprocess = ENGINE.subprocess
    original_popen = original_subprocess.Popen
    owned_process: subprocess.Popen[Any] | None = None
    benchmark_output = tempfile.TemporaryFile(mode="w+b")
    peak = 0
    termination_reason: str | None = None

    def output_bytes() -> bytes:
        benchmark_output.flush()
        benchmark_output.seek(0)
        return benchmark_output.read()

    class FileBackedBenchmarkProcess:
        """Popen facade whose completion never waits for descendant-held pipe EOF."""

        def __init__(self, process: subprocess.Popen[Any], *, encoding: str, errors: str):
            self._process = process
            self._encoding = encoding
            self._errors = errors

        @property
        def pid(self) -> int:
            return int(self._process.pid)

        @property
        def returncode(self) -> int | None:
            return self._process.returncode

        def poll(self) -> int | None:
            return self._process.poll()

        def wait(self, timeout: float | None = None) -> int:
            return int(self._process.wait(timeout=timeout))

        def kill(self) -> None:
            self._process.kill()

        def communicate(
            self, input: Any = None, timeout: float | None = None,
        ) -> tuple[str, str]:
            nonlocal termination_reason
            if input is not None:
                raise CoordinatorError("matched benchmark may not write to worker stdin")
            self._process.wait(timeout=timeout)
            group_existed, confirmed = terminate_owned_process_group_state(self._process)
            if group_existed:
                reason = termination_reason or "BENCHMARK_LEADER_EXITED_WITH_OWNED_DESCENDANTS"
                termination_reason = termination_result(
                    reason, confirmed,
                )
                raise CoordinatorError(
                    "matched benchmark leader exited while its owned process group remained: "
                    + termination_reason
                )
            return output_bytes().decode(self._encoding, errors=self._errors), ""

    def isolated_popen(*args: Any, **kwargs: Any) -> FileBackedBenchmarkProcess:
        nonlocal owned_process
        if owned_process is not None:
            raise CoordinatorError("matched benchmark attempted more than one owned R launch")
        if kwargs.get("start_new_session") is False:
            raise CoordinatorError("matched benchmark may not disable its owned process group")
        encoding = str(kwargs.get("encoding") or "utf-8")
        errors = str(kwargs.get("errors") or "replace")
        kwargs["start_new_session"] = True
        kwargs["pass_fds"] = inherited_worker_lease_fds(kwargs.get("pass_fds", ()))
        kwargs["env"] = inherited_worker_environment(kwargs.get("env"))
        # The inherited V1 engine calls communicate() after the leader exits.
        # Pipes make that wait for EOF from every descendant.  A seekable file
        # preserves all bytes without coupling coordinator progress to pipe EOF.
        kwargs["stdout"] = benchmark_output
        kwargs["stderr"] = original_subprocess.STDOUT
        owned_process = original_popen(*args, **kwargs)
        return FileBackedBenchmarkProcess(owned_process, encoding=encoding, errors=errors)

    class BenchmarkSubprocessProxy:
        """Wrap only the R launch; RSS sampler `ps` calls use untouched subprocess.run."""

        PIPE = original_subprocess.PIPE
        Popen = staticmethod(isolated_popen)
        run = staticmethod(original_subprocess.run)

    def terminate_benchmark(reason: str) -> None:
        nonlocal termination_reason
        termination_reason = reason
        if owned_process is None:
            termination_reason = termination_result(reason, False)
            return
        termination_reason = termination_result(
            reason, terminate_owned_process_group(owned_process),
        )

    def guarded_sampler(pid: int) -> int:
        nonlocal peak, termination_reason
        try:
            if owned_process is None or pid != owned_process.pid:
                raise CoordinatorError("matched benchmark sampled a process outside its owned group")
            observed = int(process_group_rss_bytes(owned_process.pid))
        except Exception as error:
            terminate_benchmark(f"RSS_SAMPLER_ERROR_{type(error).__name__}")
            raise
        peak = max(peak, observed)
        process_is_live = owned_process is not None and owned_process.poll() is None
        if observed <= 0 and process_is_live and termination_reason is None:
            terminate_benchmark("RSS_SAMPLER_RETURNED_NONPOSITIVE")
            raise CoordinatorError("matched benchmark RSS became unmeasurable while its process was live")
        if observed > maximum and termination_reason is None:
            if owned_process is None:
                raise CoordinatorError("matched benchmark process ownership was lost")
            terminate_benchmark("RAM_LIMIT_EXCEEDED")
        return observed

    started = time.monotonic()
    ENGINE.process_tree_rss_bytes = guarded_sampler
    ENGINE.subprocess = BenchmarkSubprocessProxy
    try:
        result = ENGINE.benchmark_pair(context, BENCHMARK_VARIANTS, WORKERS)
        if owned_process is not None:
            group_existed, confirmed = terminate_owned_process_group_state(owned_process)
            if group_existed:
                termination_reason = termination_result(
                    "BENCHMARK_RETURNED_WITH_OWNED_PROCESS_GROUP", confirmed,
                )
        if termination_reason is not None:
            raise CoordinatorError(
                f"matched benchmark returned after fail-closed termination: {termination_reason}"
            )
        return result
    except BaseException as error:
        if owned_process is not None:
            group_existed, confirmed = terminate_owned_process_group_state(owned_process)
            if group_existed:
                reason = termination_reason or "BENCHMARK_WRAPPER_ABORTED_WITH_OWNED_PROCESS_GROUP"
                termination_reason = termination_result(reason, confirmed)
        benchmark_failure_measurement(
            context, error, peak_bytes=peak,
            runtime_seconds=time.monotonic() - started,
            termination_reason=termination_reason,
            guard_limit_bytes=maximum,
            process_started=owned_process is not None,
            exit_status=(
                None if owned_process is None or owned_process.poll() is None
                else int(owned_process.returncode)
            ),
            process_log=output_bytes(),
        )
        raise
    finally:
        if owned_process is not None:
            terminate_owned_process_group(owned_process)
        ENGINE.subprocess = original_subprocess
        ENGINE.process_tree_rss_bytes = original_sampler
        benchmark_output.close()


def publish_benchmark_rejection_evidence(
    context: dict[str, Any], error: BaseException,
) -> Path:
    """Preserve failure evidence even when a malformed benchmark has no usable RSS row."""
    paths = ENGINE.materialized_paths(context)
    artifact_identities: dict[str, dict[str, object]] = {}
    artifact_failures: dict[str, str] = {}
    candidates = {
        "benchmark_raw": paths["benchmark_raw"],
        "benchmark": paths["benchmark"],
        "nuisance": paths["checkpoint_dir"] / "nuisance.rds",
        "nuisance_lock": paths["checkpoint_dir"] / "nuisance.sha256.tsv",
        "task": paths["task"],
    }
    for label, path in candidates.items():
        if not (path.exists() or path.is_symlink()):
            continue
        try:
            artifact_identities[label] = stable_identity(path)
        except (CoordinatorError, OSError) as artifact_error:
            artifact_failures[label] = f"{type(artifact_error).__name__}: {artifact_error}"
    payload = {
        "schema_version": "sleep-atlas-track-b-placo-benchmark-rejection.1",
        "analysis_id": context["policy"]["analysis_id"],
        "pair_id": context["pair_id"],
        "run_fingerprint": context["fingerprint"],
        "error_type": type(error).__name__,
        "error": str(error),
        "artifacts": artifact_identities,
        "artifact_identity_failures": artifact_failures,
        "terminal_gate_lock_sha256": stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK)["sha256"],
        "sequential_contract_sha256": stable_identity(ROOT / EXECUTION_CONTRACT)["sha256"],
        "coordinator_sha256": coordinator_sha256(),
        "scientific_result_accepted": False,
    }
    digest = hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
    destination = record_directory(context) / "attempts" / f"benchmark_rejected.{digest}.json"
    publish_bytes_no_replace(
        destination, json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n",
    )
    return destination


def record_benchmark_rejection(
    context: dict[str, Any], row: dict[str, str], error: BaseException,
) -> None:
    paths = ENGINE.materialized_paths(context)
    artifact_identities: dict[str, dict[str, object]] = {}
    for label, path in {
        "benchmark_raw": paths["benchmark_raw"],
        "benchmark": paths["benchmark"],
        "nuisance": paths["checkpoint_dir"] / "nuisance.rds",
        "nuisance_lock": paths["checkpoint_dir"] / "nuisance.sha256.tsv",
    }.items():
        if artifact_present(path):
            artifact_identities[label] = stable_identity(path)
    bundle_hash = benchmark_bundle_hash(paths, row)
    ram_rejected = "BLOCKED_BY_COMPUTE" in str(error)
    publish_failure_event(
        context,
        failure_class=(
            "MATCHED_BENCHMARK_RAM_ENVELOPE_REJECTION"
            if ram_rejected else "MATCHED_BENCHMARK_VALIDATION_REJECTION"
        ),
        phase_prefix=(
            "matched_worker_benchmark_ram_rejected_"
            if ram_rejected else "matched_worker_benchmark_rejected_"
        ),
        analysis=(
            "PLACO_PLUS_MATCHED_WORKER_BENCHMARK_RAM_ENVELOPE_REJECTED"
            if ram_rejected
            else "PLACO_PLUS_MATCHED_WORKER_BENCHMARK_VALIDATION_REJECTED"
        ),
        attempt_log_label="matched_worker_benchmark_rejected",
        process_log=f"{type(error).__name__}: {error}\n".encode(),
        peak_bytes=int(row["peak_process_tree_rss_bytes"]),
        runtime_seconds=float(row["wrapper_wall_seconds"]),
        public_exit_status=78 if ram_rejected else 65,
        owned_child_exit_status=0,
        rss_measurement_method=row["rss_measurement_method"],
        output_hash=bundle_hash,
        termination_reason="RAM_LIMIT_EXCEEDED_AFTER_BENCHMARK_VALIDATION" if ram_rejected else None,
        process_group_termination_status="NOT_REQUIRED_PROCESS_ALREADY_EXITED",
        terminal_compute_rejection=ram_rejected,
        details={
            "artifacts": artifact_identities,
            "benchmark_bundle_sha256": bundle_hash,
            "error": str(error),
            "error_type": type(error).__name__,
        },
    )


def ensure_benchmark(context: dict[str, Any], gate_lock: dict[str, Any]) -> dict[str, str]:
    paths = ENGINE.materialized_paths(context)
    ram_root = record_directory(context) / "ram"
    prior_ram_rejections = sorted(ram_root.glob("matched_worker_benchmark_ram_rejected_*.provenance.json"))
    for provenance_path in prior_ram_rejections:
        phase = provenance_path.name.removesuffix(".provenance.json")
        if load_ram_measurement(context, phase) is None:
            raise CoordinatorError("matched benchmark RAM rejection evidence is incomplete")
    if prior_ram_rejections:
        raise CoordinatorError(
            "BLOCKED_BY_COMPUTE immutable matched-benchmark RAM-envelope rejection evidence exists: "
            + str(prior_ram_rejections[0].relative_to(ROOT))
        )
    repair_partial_benchmark(context)
    benchmark_exists = paths["benchmark"].exists() or paths["benchmark"].is_symlink()
    if not benchmark_exists:
        if any(path.exists() for path in paths["checkpoint_dir"].glob("shard_*.rds*")):
            raise CoordinatorError("full-scan shards exist without an admitted immutable benchmark")
        try:
            require_runner_storage(context)
            run_benchmark_with_live_guard(context, gate_lock)
        except SystemExit as error:
            raise CoordinatorError(str(error)) from error
        except (OSError, BRIDGE.BridgeError) as error:
            raise CoordinatorError(str(error)) from error
    try:
        row = validate_admit_fresh_benchmark(
            context, gate_lock, require_ram_receipt=False,
        )
    except CoordinatorError as error:
        try:
            fields, rows = ENGINE.read_tsv(paths["benchmark"])
            if fields != ENGINE.BENCHMARK_RAW_FIELDS + [
                "peak_process_tree_rss_bytes", "rss_measurement_method", "wrapper_wall_seconds",
                "rss_sample_interval_seconds", "host_physical_memory_bytes", "host_free_bytes_at_benchmark",
                "peak_rss_fraction_of_physical_memory",
            ] or len(rows) != 1:
                raise CoordinatorError("rejected benchmark lacks its complete measured row")
            record_benchmark_rejection(context, rows[0], error)
        except (CoordinatorError, SystemExit, OSError, KeyError, TypeError, ValueError) as record_error:
            raise CoordinatorError(
                f"{error}; durable measured failure event could not be preserved: {record_error}"
            ) from error
        raise CoordinatorError(str(error)) from error
    paths = ENGINE.materialized_paths(context)
    output_hash = benchmark_bundle_hash(paths, row)
    record_ram(
        context, "matched_worker_benchmark",
        analysis=MATCHED_BENCHMARK_ANALYSIS,
        peak_bytes=int(row["peak_process_tree_rss_bytes"]),
        runtime_seconds=float(row["wrapper_wall_seconds"]),
        exit_status=0, output_hash=output_hash,
        rss_measurement_method=row["rss_measurement_method"],
    )
    validate_admit_fresh_benchmark(context, gate_lock, require_ram_receipt=True)
    return row


def ram_worker_limit_bytes(gate_lock: dict[str, Any]) -> int:
    reserve = int(gate_lock["resource_envelope"]["ram_safety_reserve_bytes"])
    limit = min(RAM_LIMIT_BYTES, ENGINE.physical_memory_bytes()) - reserve
    if reserve < 0 or limit <= 0:
        raise CoordinatorError("invalid PLACO+ live RAM guard envelope")
    return limit


def process_group_rss_bytes(process_group_id: int) -> int:
    """Measure every process in the owned fresh-session process group.

    PPID trees can lose a still-live descendant after reparenting. Every R
    launch owns a fresh session whose PGID is the leader PID, so PGID
    aggregation preserves the termination and measurement boundary exactly.
    """

    if process_group_id <= 0:
        raise CoordinatorError("invalid owned process-group identifier")
    result = subprocess.run(
        ["ps", "-axo", "pid=,pgid=,rss="], check=False,
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise CoordinatorError(
            f"ps failed while measuring owned process group: {result.stderr.strip()}"
        )
    total_kib = 0
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            _, pgid, rss_kib = (int(value) for value in fields)
        except ValueError:
            continue
        if pgid == process_group_id and rss_kib >= 0:
            total_kib += rss_kib
    return total_kib * 1024


def terminate_owned_process_group_state(process: subprocess.Popen[Any]) -> tuple[bool, bool]:
    """Terminate a fresh session, including descendants after leader exit.

    Signal-zero process-group probes are not portable under managed macOS
    execution.  Repeated real signals give us a stronger result: ESRCH proves
    the group disappeared; EPERM or a group still accepting SIGKILL after the
    bounded grace period is an unconfirmed termination and must fail closed.
    """
    process_group_id = process.pid
    try:
        os.killpg(process_group_id, signal.SIGTERM)
    except ProcessLookupError:
        process.poll()
        return False, True
    except PermissionError:
        return True, False
    deadline = time.monotonic() + 0.5
    while time.monotonic() < deadline:
        process.poll()
        time.sleep(0.05)
    try:
        os.killpg(process_group_id, signal.SIGKILL)
    except ProcessLookupError:
        process.poll()
        return True, True
    except PermissionError:
        return True, False
    deadline = time.monotonic() + 1.0
    while time.monotonic() < deadline:
        process.poll()
        try:
            os.killpg(process_group_id, signal.SIGKILL)
        except ProcessLookupError:
            process.poll()
            return True, True
        except PermissionError:
            return True, False
        time.sleep(0.05)
    process.poll()
    return True, False


def terminate_owned_process_group(process: subprocess.Popen[Any]) -> bool:
    return terminate_owned_process_group_state(process)[1]


def termination_result(reason: str, confirmed: bool) -> str:
    return reason if confirmed else f"{reason}_PROCESS_GROUP_TERMINATION_UNCONFIRMED"


def monitored_command(
    command: list[str], *, maximum_process_tree_rss_bytes: int | None = None,
) -> dict[str, Any]:
    started = time.monotonic()
    peak = 0
    termination_reason: str | None = None
    process: subprocess.Popen[bytes] | None = None
    with tempfile.TemporaryFile(mode="w+b") as output:
        try:
            process = subprocess.Popen(
                command, cwd=ROOT, stdout=output, stderr=subprocess.STDOUT,
                start_new_session=True,
                pass_fds=inherited_worker_lease_fds(),
                env=inherited_worker_environment(),
            )
        except OSError as error:
            termination_reason = f"PROCESS_START_ERROR_{type(error).__name__}"
            returncode = 70
            output.write(f"{type(error).__name__}: {error}\n".encode())
        if process is not None:
            while process.poll() is None:
                try:
                    observed = int(process_group_rss_bytes(process.pid))
                except Exception as error:
                    termination_reason = f"RSS_SAMPLER_ERROR_{type(error).__name__}"
                    output.write(f"{type(error).__name__}: {error}\n".encode())
                    termination_reason = termination_result(
                        termination_reason, terminate_owned_process_group(process),
                    )
                    break
                if observed <= 0 and process.poll() is None:
                    termination_reason = "RSS_SAMPLER_RETURNED_NONPOSITIVE"
                    termination_reason = termination_result(
                        termination_reason, terminate_owned_process_group(process),
                    )
                    break
                peak = max(peak, observed)
                if (
                    maximum_process_tree_rss_bytes is not None
                    and peak > maximum_process_tree_rss_bytes
                ):
                    termination_reason = "RAM_LIMIT_EXCEEDED"
                    termination_reason = termination_result(
                        termination_reason, terminate_owned_process_group(process),
                    )
                    break
                time.sleep(0.05)
            if termination_reason is None:
                group_existed, confirmed = terminate_owned_process_group_state(process)
                if group_existed:
                    termination_reason = termination_result(
                        "OWNED_PROCESS_GROUP_OUTLIVED_LEADER", confirmed,
                    )
            try:
                returncode = int(process.wait(timeout=1.0))
            except subprocess.TimeoutExpired:
                reason = termination_reason or "OWNED_PROCESS_GROUP_OUTLIVED_MONITOR"
                termination_reason = termination_result(
                    reason, terminate_owned_process_group(process),
                )
                try:
                    returncode = int(process.wait(timeout=1.0))
                except subprocess.TimeoutExpired:
                    returncode = 70
            if termination_reason is not None and returncode == 0:
                returncode = 70
        wall = time.monotonic() - started
        output.seek(0)
        log = output.read()
    if termination_reason is not None:
        log += (
            "\nTRACK_B_PLACO_LIVE_RSS_MONITOR_TERMINATED "
            f"reason={termination_reason} peak_bytes={peak} "
            f"limit_bytes={maximum_process_tree_rss_bytes}\n"
        ).encode()
    if peak <= 0:
        if process is None:
            peak = 0
            method = "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE"
        else:
            observed = int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
            peak = observed if sys.platform == "darwin" else observed * 1024
            method = (
                "CHILD_RUSAGE_SESSION_MAXRSS_FALLBACK_PS_UNAVAILABLE"
                if peak > 0
                else "NO_PROCESS_RSS_OBSERVED_BEFORE_FAILURE"
            )
    else:
        method = "AGGREGATE_PROCESS_TREE_PS_SAMPLED"
    return {
        "exit_status": returncode,
        "peak_bytes": peak,
        "runtime_seconds": wall,
        "rss_measurement_method": method,
        "ram_guard_triggered": bool(
            termination_reason and termination_reason.startswith("RAM_LIMIT_EXCEEDED")
        ),
        "termination_reason": termination_reason,
        "ram_guard_limit_bytes": maximum_process_tree_rss_bytes,
        "log": log,
    }


def publish_attempt_log(context: dict[str, Any], phase: str, measurement: dict[str, Any]) -> Path:
    content = bytes(measurement["log"])
    if not content:
        content = b"NO_PROCESS_OUTPUT\n"
    digest = hashlib.sha256(content).hexdigest()
    path = record_directory(context) / "attempts" / f"{phase}.{digest}.log"
    publish_bytes_no_replace(path, content)
    return path


def publish_staged_pair(
    context: dict[str, Any], gate_lock: dict[str, Any], workers: int,
) -> dict[str, Any]:
    paths = ENGINE.materialized_paths(context)
    published: list[tuple[Path, int, int]] = []
    original_exclusive_family = ENGINE.exclusive_family

    def publication_boundary(sources: list[Path], destinations: list[Path]) -> None:
        BRIDGE.exclusive_family_v2(
            sources,
            destinations,
            before_publication=lambda: BRIDGE.revalidate_publication_state(
                context, gate_lock, paths, workers,
            ),
            after_publication=lambda: BRIDGE.revalidate_publication_state(
                context, gate_lock, paths, workers,
            ),
            committed=published,
        )

    try:
        BRIDGE.revalidate_publication_state(context, gate_lock, paths, workers)
        ENGINE.exclusive_family = publication_boundary
        provenance = ENGINE.publish_run(ROOT, paths["task"], paths["run_summary"])
        if provenance.get("reference_sha256") != (
            "NOT_APPLICABLE_TO_PLACO_FULL_P_SCAN_LD_REQUIRED_FOR_LATER_CLUMPING"
        ):
            raise CoordinatorError("full-P result incorrectly claims LD-backed locus publication")
        if provenance.get("execution_resource_envelope", {}).get("observed_workers") != workers:
            raise CoordinatorError("published PLACO+ result worker count differs from benchmark")
        BRIDGE.revalidate_publication_state(context, gate_lock, paths, workers)
        return provenance
    except (SystemExit, GATE.GateError, BRIDGE.BridgeError) as error:
        try:
            BRIDGE.rollback_v2_publication(published)
        except BRIDGE.BridgeError as rollback_error:
            raise CoordinatorError(str(rollback_error)) from error
        raise CoordinatorError(str(error)) from error
    except BaseException:
        BRIDGE.rollback_v2_publication(published)
        raise
    finally:
        ENGINE.exclusive_family = original_exclusive_family


def run_full_scan_monitored(
    context: dict[str, Any], gate_lock: dict[str, Any], benchmark: dict[str, str],
) -> dict[str, Any]:
    paths = ENGINE.materialized_paths(context)
    prior_rejections = full_scan_rejection_phases(context)
    for phase in prior_rejections:
        if load_ram_measurement(context, phase) is None:
            raise CoordinatorError("full-scan RAM rejection evidence is incomplete")
    if prior_rejections:
        raise CoordinatorError(
            "BLOCKED_BY_COMPUTE immutable full-scan RAM-envelope rejection evidence exists: "
            + ",".join(prior_rejections)
        )
    stage_present = [
        paths["staged_ledger"].exists() or paths["staged_ledger"].is_symlink(),
        paths["run_summary"].exists() or paths["run_summary"].is_symlink(),
    ]
    if all(stage_present):
        existing = load_ram_measurement(context, "full_scan")
        if existing is not None:
            try:
                BRIDGE.validate_stage_against_benchmark(paths, benchmark, WORKERS)
            except BRIDGE.BridgeError as error:
                raise CoordinatorError(str(error)) from error
            staged_hash = stable_identity(paths["staged_ledger"])["sha256"]
            require_full_scan_ram_evidence(context, gate_lock, str(staged_hash))
            rebuild_ram_aggregate()
            return publish_staged_pair(context, gate_lock, WORKERS)
        rejected_receipts = sorted(
            (record_directory(context) / "ram").glob("full_scan_ram_rejected_*.tsv")
        )
        if rejected_receipts:
            raise CoordinatorError(
                "BLOCKED_BY_COMPUTE complete staged PLACO+ output has immutable prior "
                f"RAM-envelope rejection evidence: {rejected_receipts[0].relative_to(ROOT)}"
            )
        if (
            paths["staged_provenance"].exists()
            or paths["staged_provenance"].is_symlink()
            or canonical_state(paths) != "ABSENT"
        ):
            raise CoordinatorError("unmeasured staged PLACO+ family coexists with publication evidence")
        audit_and_unlink_regenerable(
            context, "complete_stage_without_monitor_receipt",
            [paths["staged_ledger"], paths["run_summary"]],
            "COORDINATOR_DIED_AFTER_EXACT_STAGE_LINKS_BUT_BEFORE_RAM_RECEIPT;"
            "REBUILD_FROM_FROZEN_GLOBAL_NUISANCE_AND_COMPLETE_SHARDS_UNDER_LIVE_MONITOR",
        )
    purge_unmeasured_full_scan_checkpoints(context)
    try:
        require_runner_storage(context)
        ENGINE.verify_materialized(context)
    except SystemExit as error:
        raise CoordinatorError(str(error)) from error
    command = [
        str(ROOT / ".r-env/bin/Rscript"), str(ROOT / ENGINE.RUNNER),
        ENGINE.relative(ROOT, paths["task"]), "--root", str(ROOT), "--execute",
        "--workers", str(WORKERS), "--stage-only",
    ]
    measurement = monitored_command(
        command,
        maximum_process_tree_rss_bytes=ram_worker_limit_bytes(gate_lock),
    )
    if measurement["exit_status"] != 0:
        ram_rejected = bool(measurement.get("ram_guard_triggered"))
        phase_prefix = FULL_SCAN_REJECTION_PREFIX if ram_rejected else "full_scan_failed_"
        analysis = (
            "PLACO_PLUS_FULL_P_SCAN_RAM_ENVELOPE_REJECTED"
            if ram_rejected else "PLACO_PLUS_FULL_P_SCAN_FAILED_ATTEMPT"
        )
        event = publish_failure_event(
            context,
            failure_class="FULL_SCAN_PROCESS_FAILURE",
            phase_prefix=phase_prefix,
            analysis=analysis,
            attempt_log_label="full_scan_failed",
            process_log=bytes(measurement["log"]),
            peak_bytes=int(measurement["peak_bytes"]),
            runtime_seconds=float(measurement["runtime_seconds"]),
            public_exit_status=(78 if ram_rejected else int(measurement["exit_status"])),
            owned_child_exit_status=int(measurement["exit_status"]),
            rss_measurement_method=str(measurement["rss_measurement_method"]),
            output_hash=None,
            termination_reason=measurement.get("termination_reason"),
            process_group_termination_status=failure_group_termination_status(
                measurement.get("termination_reason"),
                process_started=not str(measurement.get("termination_reason", "")).startswith(
                    "PROCESS_START_ERROR_"
                ),
            ),
            terminal_compute_rejection=ram_rejected,
            details={
                "ram_guard_limit_bytes": measurement.get("ram_guard_limit_bytes"),
            },
        )
        raise CoordinatorError(
            f"PLACO+ full-P runner failed for {context['pair_id']}; "
            f"durable_failure_event={event['event_id']}"
        )
    try:
        BRIDGE.validate_stage_against_benchmark(paths, benchmark, WORKERS)
    except BRIDGE.BridgeError as error:
        publish_failure_event(
            context,
            failure_class="FULL_SCAN_INVALID_STAGE",
            phase_prefix="full_scan_invalid_stage_",
            analysis="PLACO_PLUS_FULL_P_SCAN_INVALID_STAGE_FAILED_ATTEMPT",
            attempt_log_label="full_scan_invalid_stage",
            process_log=(
                bytes(measurement["log"])
                + f"\n{type(error).__name__}: {error}\n".encode()
            ),
            peak_bytes=int(measurement["peak_bytes"]),
            runtime_seconds=float(measurement["runtime_seconds"]),
            public_exit_status=65,
            owned_child_exit_status=0,
            rss_measurement_method=str(measurement["rss_measurement_method"]),
            output_hash=None,
            termination_reason=None,
            process_group_termination_status="NOT_REQUIRED_PROCESS_ALREADY_EXITED",
            terminal_compute_rejection=False,
            details={"validation_error": str(error)},
        )
        raise CoordinatorError(str(error)) from error
    staged_hash = stable_identity(paths["staged_ledger"])["sha256"]
    observed_benchmark = {
        **benchmark,
        "workers": str(WORKERS),
        "peak_process_tree_rss_bytes": str(measurement["peak_bytes"]),
    }
    try:
        admission = admit_ram(observed_benchmark, gate_lock)
    except CoordinatorError as error:
        publish_failure_event(
            context,
            failure_class="FULL_SCAN_RAM_ENVELOPE_REJECTION",
            phase_prefix=FULL_SCAN_REJECTION_PREFIX,
            analysis="PLACO_PLUS_FULL_P_SCAN_RAM_ENVELOPE_REJECTED",
            attempt_log_label="full_scan_ram_rejected",
            process_log=(
                bytes(measurement["log"])
                + f"\n{type(error).__name__}: {error}\n".encode()
            ),
            peak_bytes=int(measurement["peak_bytes"]),
            runtime_seconds=float(measurement["runtime_seconds"]),
            public_exit_status=78,
            owned_child_exit_status=0,
            rss_measurement_method=str(measurement["rss_measurement_method"]),
            output_hash=str(staged_hash),
            termination_reason="RAM_LIMIT_EXCEEDED_AFTER_FULL_SCAN",
            process_group_termination_status="NOT_REQUIRED_PROCESS_ALREADY_EXITED",
            terminal_compute_rejection=True,
            details={"staged_output_sha256": str(staged_hash)},
        )
        raise
    publish_attempt_log(context, "full_scan", measurement)
    record_ram(
        context, "full_scan",
        analysis=FULL_SCAN_ANALYSIS,
        peak_bytes=admission["peak_bytes"],
        runtime_seconds=float(measurement["runtime_seconds"]),
        exit_status=0, output_hash=str(staged_hash),
        rss_measurement_method=str(measurement["rss_measurement_method"]),
    )
    return publish_staged_pair(context, gate_lock, WORKERS)


def canonical_state(paths: dict[str, Path]) -> str:
    present = [
        paths["canonical_ledger"].exists() or paths["canonical_ledger"].is_symlink(),
        paths["canonical_provenance"].exists() or paths["canonical_provenance"].is_symlink(),
    ]
    if all(present):
        return "COMPLETE"
    if any(present):
        return "PARTIAL"
    return "ABSENT"


def validate_completed_pair(
    context: dict[str, Any], gate_lock: dict[str, Any], benchmark: dict[str, str] | None = None,
) -> dict[str, Any]:
    paths = ENGINE.materialized_paths(context)
    if canonical_state(paths) != "COMPLETE":
        raise CoordinatorError("canonical PLACO+ pair family is not complete")
    try:
        ENGINE.verify_materialized(context)
        validated_benchmark = validate_admit_fresh_benchmark(
            context, gate_lock, require_ram_receipt=True,
        )
        if benchmark is not None and benchmark != validated_benchmark:
            raise CoordinatorError("caller benchmark differs from the freshly validated benchmark")
        benchmark = validated_benchmark
        BRIDGE.validate_stage_against_benchmark(paths, benchmark, WORKERS)
        task = ENGINE.one_tsv_row(paths["task"])
        summary = ENGINE.one_tsv_row(paths["run_summary"])
        counts = ENGINE.validate_staged_ledger(ROOT, task, summary)
    except (SystemExit, BRIDGE.BridgeError) as error:
        raise CoordinatorError(str(error)) from error
    ledger = stable_identity(paths["canonical_ledger"])
    staged = stable_identity(paths["staged_ledger"])
    canonical_provenance = stable_identity(paths["canonical_provenance"])
    staged_provenance = stable_identity(paths["staged_provenance"])
    if (ledger["device"], ledger["inode"]) != (staged["device"], staged["inode"]):
        raise CoordinatorError("canonical ledger is not the exact deeply validated staged inode")
    if (canonical_provenance["device"], canonical_provenance["inode"]) != (
        staged_provenance["device"], staged_provenance["inode"],
    ):
        raise CoordinatorError("canonical provenance is not the exact staged provenance inode")
    full_scan_ram = require_full_scan_ram_evidence(
        context, gate_lock, str(ledger["sha256"]),
    )
    provenance = read_json(paths["canonical_provenance"])
    expected = {
        "schema_version": "sleep-atlas-track-b-placo-result.1",
        "analysis_id": task["analysis_id"],
        "pair_id": context["pair_id"],
        "run_fingerprint": context["fingerprint"],
        "input_sha256": task["aligned_input_sha256"],
        "input_rows": int(task["aligned_input_rows"]),
        "output": task["canonical_ledger"],
        "output_sha256": ledger["sha256"],
        "output_bytes": ledger["bytes"],
        "output_rows": counts["rows"],
        "exact_schema": ",".join(ENGINE.LEDGER_FIELDS),
        "within_pair_bh_family_n": counts["rows"],
        "terminal_result_state": counts["terminal_status"],
        "qc_status": "PASS",
        "reference_sha256": "NOT_APPLICABLE_TO_PLACO_FULL_P_SCAN_LD_REQUIRED_FOR_LATER_CLUMPING",
    }
    if any(provenance.get(key) != value for key, value in expected.items()):
        raise CoordinatorError("canonical PLACO+ provenance differs from the validated complete ledger")
    if provenance.get("complete_family_counts") != counts:
        raise CoordinatorError("canonical PLACO+ complete-family counts drifted")
    nuisance = provenance.get("nuisance_estimation", {})
    if (
        nuisance.get("scope") != "ALL_VALID_ALIGNED_GENOME_WIDE_VARIANTS_BEFORE_ANY_SINGLE_VARIANT_SHARD"
        or nuisance.get("input_rows") != counts["rows"]
        or provenance.get("execution_resource_envelope", {}).get("observed_workers") != WORKERS
    ):
        raise CoordinatorError("canonical PLACO+ nuisance scope or worker provenance drifted")
    try:
        if GATE.verify_gate() != gate_lock:
            raise CoordinatorError("terminal gate changed before storage cleanup")
    except GATE.GateError as error:
        raise CoordinatorError(str(error)) from error
    return {"provenance": provenance, "counts": counts, "canonical_ledger": ledger,
            "canonical_provenance": canonical_provenance, "full_scan_ram": full_scan_ram}


def same_file_identity(left: dict[str, Any], right: dict[str, Any]) -> bool:
    return all(left.get(key) == right.get(key) for key in ("device", "inode", "bytes", "sha256"))


def publication_repair_receipt(
    context: dict[str, Any], action: str, artifacts: list[dict[str, object]],
) -> Path:
    payload = {
        "schema_version": "sleep-atlas-track-b-placo-publication-repair.1",
        "analysis_id": context["policy"]["analysis_id"],
        "pair_id": context["pair_id"],
        "run_fingerprint": context["fingerprint"],
        "action": action,
        "artifacts": artifacts,
        "scientific_data_deleted": False,
        "rule": "ONLY_EXACT_HARDLINKS_TO_DEEPLY_VALIDATABLE_STAGED_ARTIFACTS_ARE_RECONCILED",
        "terminal_gate_lock_sha256": stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK)["sha256"],
        "sequential_contract_sha256": stable_identity(ROOT / EXECUTION_CONTRACT)["sha256"],
        "coordinator_sha256": coordinator_sha256(),
    }
    digest = hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
    path = record_directory(context) / "repairs" / f"canonical_publication.{digest}.json"
    publish_bytes_no_replace(path, json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n")
    return path


def reconcile_partial_canonical_publication(
    context: dict[str, Any], gate_lock: dict[str, Any], contract: dict[str, Any],
) -> str:
    """Recover only the two possible SIGKILL prefixes of V1's three-link commit."""
    if read_json(ROOT / EXECUTION_CONTRACT) != contract:
        raise CoordinatorError("sequential contract changed before canonical reconciliation")
    paths = ENGINE.materialized_paths(context)
    canonical_ledger = paths["canonical_ledger"]
    canonical_provenance = paths["canonical_provenance"]
    staged_provenance = paths["staged_provenance"]
    present = tuple(
        path.exists() or path.is_symlink()
        for path in (canonical_ledger, canonical_provenance, staged_provenance)
    )
    if present == (False, False, False):
        return "ABSENT"
    if present == (True, True, True):
        return "COMPLETE"
    # The frozen publisher links ledger, canonical provenance, then staged
    # provenance.  Every other pattern is corruption rather than a crash prefix.
    if present not in {(True, False, False), (True, True, False)}:
        raise CoordinatorError(
            f"canonical PLACO+ publication is not an exact crash prefix for {context['pair_id']}"
        )
    try:
        benchmark = validate_admit_fresh_benchmark(
            context, gate_lock, require_ram_receipt=True,
        )
        BRIDGE.revalidate_publication_state(context, gate_lock, paths, WORKERS)
    except (BRIDGE.BridgeError, GATE.GateError, SystemExit, CoordinatorError) as error:
        raise CoordinatorError(f"partial canonical publication failed prerequisite validation: {error}") from error
    stage_ledger_identity = stable_identity(paths["staged_ledger"])
    canonical_ledger_identity = stable_identity(canonical_ledger)
    if not same_file_identity(stage_ledger_identity, canonical_ledger_identity):
        raise CoordinatorError("partial canonical ledger is not the exact staged hardlink")
    require_full_scan_ram_evidence(
        context, gate_lock, str(stage_ledger_identity["sha256"]),
    )
    if present == (True, False, False):
        publication_repair_receipt(
            context, "UNLINK_EXACT_FIRST_HARDLINK_AND_REENTER_DEEP_PUBLICATION",
            [stage_ledger_identity, canonical_ledger_identity],
        )
        safe_repo_path(canonical_ledger, "partial canonical ledger")
        current = os.lstat(canonical_ledger)
        if (current.st_dev, current.st_ino, current.st_size) != (
            canonical_ledger_identity["device"], canonical_ledger_identity["inode"],
            canonical_ledger_identity["bytes"],
        ):
            raise CoordinatorError("partial canonical ledger changed before safe reconciliation")
        canonical_ledger.unlink()
        fsync_directory(canonical_ledger.parent)
        if not identity_matches(paths["staged_ledger"], stage_ledger_identity):
            raise CoordinatorError("staged ledger changed while reconciling its exact canonical hardlink")
        return "ABSENT"

    canonical_provenance_identity = stable_identity(canonical_provenance)
    repair = publication_repair_receipt(
        context, "LINK_MISSING_THIRD_HARDLINK_FROM_EXACT_CANONICAL_PROVENANCE",
        [stage_ledger_identity, canonical_ledger_identity, canonical_provenance_identity],
    )
    safe_repo_path(canonical_provenance, "canonical provenance repair source")
    safe_repo_path(staged_provenance.parent, "staged provenance repair parent")
    try:
        os.link(canonical_provenance, staged_provenance, follow_symlinks=False)
    except FileExistsError as error:
        raise CoordinatorError("staged provenance appeared during canonical reconciliation") from error
    linked = os.lstat(staged_provenance)
    if (
        not stat.S_ISREG(linked.st_mode)
        or (linked.st_dev, linked.st_ino, linked.st_size) != (
            canonical_provenance_identity["device"], canonical_provenance_identity["inode"],
            canonical_provenance_identity["bytes"],
        )
    ):
        try:
            staged_provenance.unlink()
        except FileNotFoundError:
            pass
        raise CoordinatorError("reconciled staged provenance is not the exact canonical hardlink")
    fsync_directory(staged_provenance.parent)
    try:
        validate_completed_pair(context, gate_lock, benchmark)
        try:
            if GATE.verify_gate() != gate_lock:
                raise CoordinatorError("terminal gate changed during canonical reconciliation")
        except GATE.GateError as error:
            raise CoordinatorError(str(error)) from error
    except BaseException:
        current = os.lstat(staged_provenance)
        if (current.st_dev, current.st_ino) == (linked.st_dev, linked.st_ino):
            staged_provenance.unlink()
            fsync_directory(staged_provenance.parent)
        raise
    if not repair.exists():
        raise CoordinatorError("canonical reconciliation audit receipt disappeared")
    return "COMPLETE"


def allowed_cleanup_relative(relative: str) -> bool:
    return (
        relative in FIXED_REGENERABLE_ARTIFACTS
        or SHARD_RE.fullmatch(relative) is not None
        or KNOWN_TEMP_RE.fullmatch(relative) is not None
    )


def collect_cleanup_candidates(base: Path) -> list[dict[str, object]]:
    safe_repo_path(base, "PLACO pair work root")
    candidates: list[dict[str, object]] = []
    if not base.exists():
        return candidates
    for path in sorted(base.rglob("*")):
        observed = os.lstat(path)
        if stat.S_ISDIR(observed.st_mode):
            relative_directory = path.relative_to(base).as_posix()
            if relative_directory not in KNOWN_WORK_DIRECTORIES:
                raise CoordinatorError(
                    f"unknown PLACO work directory is not cleanup-authorized: {relative_directory}"
                )
            continue
        if not stat.S_ISREG(observed.st_mode):
            raise CoordinatorError(f"PLACO work tree contains a non-regular cleanup candidate: {path}")
        relative = path.relative_to(base).as_posix()
        if not allowed_cleanup_relative(relative):
            raise CoordinatorError(f"unknown PLACO work artifact is not cleanup-authorized: {relative}")
        candidates.append(stable_identity(path, allow_empty=KNOWN_TEMP_RE.fullmatch(relative) is not None))
    return candidates


def cleanup_reclaim_accounting(
    candidates: list[dict[str, object]],
) -> dict[str, object]:
    """Describe logical deletion and physical reclaim without counting hard links twice."""
    groups: dict[tuple[int, int], dict[str, object]] = {}
    logical_bytes = 0
    for record in candidates:
        path = ROOT / str(record["path"])
        current = os.lstat(path)
        if (
            not stat.S_ISREG(current.st_mode)
            or (current.st_dev, current.st_ino, current.st_size)
            != (record["device"], record["inode"], record["bytes"])
        ):
            raise CoordinatorError("cleanup reclaim candidate changed during authorization")
        logical_bytes += int(record["bytes"])
        key = int(record["device"]), int(record["inode"])
        group = groups.setdefault(key, {
            "device": key[0], "inode": key[1], "bytes": int(record["bytes"]),
            "authorized_work_paths": [], "link_count_at_authorization": int(current.st_nlink),
        })
        if (
            group["bytes"] != int(record["bytes"])
            or group["link_count_at_authorization"] != int(current.st_nlink)
        ):
            raise CoordinatorError("hard-linked cleanup candidates have inconsistent identities")
        group["authorized_work_paths"].append(str(record["path"]))
    unique_bytes = sum(int(group["bytes"]) for group in groups.values())
    reclaimable = 0
    externally_retained = 0
    rows: list[dict[str, object]] = []
    for group in sorted(groups.values(), key=lambda item: (int(item["device"]), int(item["inode"]))):
        paths = list(group["authorized_work_paths"])
        fully_reclaimed = int(group["link_count_at_authorization"]) == len(paths)
        group["reclaimable_after_all_authorized_unlinks"] = fully_reclaimed
        if fully_reclaimed:
            reclaimable += int(group["bytes"])
        else:
            externally_retained += int(group["bytes"])
        rows.append(group)
    return {
        "authorized_logical_bytes": logical_bytes,
        "unique_work_inode_bytes": unique_bytes,
        "logical_hardlink_alias_bytes": logical_bytes - unique_bytes,
        "physically_reclaimable_bytes_at_authorization": reclaimable,
        "externally_retained_hardlink_bytes_at_authorization": externally_retained,
        "inode_groups": rows,
    }


def validate_cleanup_reclaim_accounting(
    candidates: list[dict[str, Any]], accounting: Any,
) -> None:
    if not isinstance(accounting, dict) or not isinstance(accounting.get("inode_groups"), list):
        raise CoordinatorError("cleanup bundle lacks hardlink-aware reclaim accounting")
    candidate_groups: dict[tuple[int, int], dict[str, object]] = {}
    logical_bytes = 0
    for record in candidates:
        try:
            key = int(record["device"]), int(record["inode"])
            size = int(record["bytes"])
            path = str(record["path"])
        except (KeyError, TypeError, ValueError) as error:
            raise CoordinatorError("cleanup reclaim accounting has a malformed candidate") from error
        logical_bytes += size
        group = candidate_groups.setdefault(key, {"bytes": size, "paths": []})
        if int(group["bytes"]) != size:
            raise CoordinatorError("cleanup candidates disagree about one hard-linked inode")
        group["paths"].append(path)
    observed_groups: dict[tuple[int, int], dict[str, Any]] = {}
    for row in accounting["inode_groups"]:
        if not isinstance(row, dict):
            raise CoordinatorError("cleanup reclaim accounting contains a malformed inode group")
        try:
            key = int(row["device"]), int(row["inode"])
            size = int(row["bytes"])
            link_count = int(row["link_count_at_authorization"])
            paths = list(row["authorized_work_paths"])
        except (KeyError, TypeError, ValueError) as error:
            raise CoordinatorError("cleanup reclaim accounting contains invalid fields") from error
        if key in observed_groups or any(not isinstance(path, str) for path in paths):
            raise CoordinatorError("cleanup reclaim accounting repeats an inode or has invalid paths")
        wanted = candidate_groups.get(key)
        if (
            wanted is None
            or size != int(wanted["bytes"])
            or sorted(paths) != sorted(wanted["paths"])
            or link_count < len(paths)
            or row.get("reclaimable_after_all_authorized_unlinks") != (link_count == len(paths))
        ):
            raise CoordinatorError("cleanup reclaim accounting differs from authorized candidates")
        observed_groups[key] = row
    if set(observed_groups) != set(candidate_groups):
        raise CoordinatorError("cleanup reclaim accounting omits an authorized inode")
    unique_bytes = sum(int(group["bytes"]) for group in candidate_groups.values())
    reclaimable = sum(
        int(row["bytes"]) for row in observed_groups.values()
        if row["reclaimable_after_all_authorized_unlinks"] is True
    )
    externally_retained = unique_bytes - reclaimable
    expected = {
        "authorized_logical_bytes": logical_bytes,
        "unique_work_inode_bytes": unique_bytes,
        "logical_hardlink_alias_bytes": logical_bytes - unique_bytes,
        "physically_reclaimable_bytes_at_authorization": reclaimable,
        "externally_retained_hardlink_bytes_at_authorization": externally_retained,
    }
    if any(accounting.get(key) != value for key, value in expected.items()):
        raise CoordinatorError("cleanup reclaim byte totals drifted")


def cleanup_bundle_path(context: dict[str, Any]) -> Path:
    return record_directory(context) / "reproducibility_and_cleanup.plan.json"


def cleanup_complete_path(context: dict[str, Any]) -> Path:
    return record_directory(context) / "cleanup.complete.json"


def retained_work_storage_forecast(
    context: dict[str, Any], gate_lock: dict[str, Any],
) -> dict[str, object]:
    base = ENGINE.materialized_paths(context)["base"]
    logical_retained = 0
    inode_groups: dict[tuple[int, int], dict[str, int]] = {}
    unknown_artifacts: list[str] = []
    unsafe_artifacts: list[str] = []
    if base.exists():
        for path in base.rglob("*"):
            observed = os.lstat(path)
            relative = path.relative_to(base).as_posix()
            if stat.S_ISDIR(observed.st_mode):
                if relative not in KNOWN_WORK_DIRECTORIES:
                    unknown_artifacts.append(relative + "/")
                continue
            if stat.S_ISREG(observed.st_mode):
                logical_retained += int(observed.st_size)
                key = int(observed.st_dev), int(observed.st_ino)
                group = inode_groups.setdefault(key, {
                    "bytes": int(observed.st_size), "work_links": 0,
                    "total_links": int(observed.st_nlink),
                })
                group["work_links"] += 1
                if not allowed_cleanup_relative(relative):
                    unknown_artifacts.append(relative)
            else:
                unsafe_artifacts.append(relative)
    unique_physical = sum(group["bytes"] for group in inode_groups.values())
    physically_reclaimable = sum(
        group["bytes"] for group in inode_groups.values()
        if group["work_links"] == group["total_links"]
    )
    externally_retained = unique_physical - physically_reclaimable
    index = PAIR_ORDER.index(context["pair_id"])
    next_pair = PAIR_ORDER[index + 1] if index + 1 < len(PAIR_ORDER) else None
    free = shutil.disk_usage(ROOT).free
    required = 0
    if next_pair is not None:
        trait1, trait2, _ = ENGINE.PAIR_IDENTITIES[next_pair]
        dense = gate_lock.get("frozen_v1", {}).get("dense_sources", {})
        try:
            source_rows = int(dense[trait1]["rows"]) + int(dense[trait2]["rows"])
        except (KeyError, TypeError, ValueError) as error:
            raise CoordinatorError("terminal gate lacks next-pair storage identities") from error
        required = int(ENGINE.MINIMUM_FREE_BYTES + 600 * source_rows)
    return {
        "retained_work_logical_bytes": logical_retained,
        "retained_work_unique_inode_bytes": unique_physical,
        "logical_hardlink_alias_bytes": logical_retained - unique_physical,
        "physically_reclaimable_work_bytes": physically_reclaimable,
        "hardlinked_bytes_not_reclaimable": externally_retained,
        "unknown_work_artifact_count": len(unknown_artifacts),
        "unknown_work_artifacts": unknown_artifacts,
        "unsafe_work_artifact_count": len(unsafe_artifacts),
        "unsafe_work_artifacts": unsafe_artifacts,
        "free_bytes": free,
        "next_pair": next_pair or "NONE",
        "next_pair_required_bytes": required,
        "additional_storage_required_if_cleanup_remains_blocked": max(0, required - free),
    }


def cleanup_failure(
    error: BaseException, context: dict[str, Any], gate_lock: dict[str, Any],
) -> CoordinatorError:
    forecast = retained_work_storage_forecast(context, gate_lock)
    detail = " ".join(f"{key}={value}" for key, value in forecast.items())
    return CoordinatorError(f"cleanup failed closed; no unlisted data removed: {error}; {detail}")


def build_cleanup_bundle(
    context: dict[str, Any], validated: dict[str, Any], contract: dict[str, Any],
) -> dict[str, Any]:
    paths = ENGINE.materialized_paths(context)
    full_scan_ram = validated.get("full_scan_ram")
    if (
        not isinstance(full_scan_ram, dict)
        or full_scan_ram.get("phase") != "full_scan"
        or full_scan_ram.get("row", {}).get("analysis") != FULL_SCAN_ANALYSIS
        or full_scan_ram.get("exit_status") != 0
        or full_scan_ram.get("output_hash") != validated.get("canonical_ledger", {}).get("sha256")
    ):
        raise CoordinatorError("validated pair lacks exact successful full-scan RAM evidence")
    embedded: dict[str, dict[str, object]] = {}
    total = 0
    for relative in sorted(PRESERVED_WORK_ARTIFACTS | OPTIONAL_PRESERVED_WORK_ARTIFACTS):
        path = paths["base"] / relative
        if not path.exists():
            if relative in OPTIONAL_PRESERVED_WORK_ARTIFACTS:
                continue
            raise CoordinatorError(f"required reproducibility artifact is missing before cleanup: {relative}")
        identity = stable_identity(path)
        content = path.read_bytes()
        if hashlib.sha256(content).hexdigest() != identity["sha256"]:
            raise CoordinatorError(f"reproducibility artifact changed while embedding: {relative}")
        total += len(content)
        if total > MAX_EMBEDDED_PROVENANCE_BYTES:
            raise CoordinatorError("small-artifact reproducibility bundle unexpectedly exceeds 64 MiB")
        embedded[relative] = {
            "bytes": identity["bytes"], "sha256": identity["sha256"],
            "base64": base64.b64encode(content).decode("ascii"),
        }
    rebuild_ram_aggregate()
    ram_receipts = sorted((record_directory(context) / "ram").glob("*.tsv"))
    if not ram_receipts:
        raise CoordinatorError("RAM evidence is missing before storage cleanup")
    ram_provenance = [path.with_suffix(".provenance.json") for path in ram_receipts]
    if any(not path.exists() or path.is_symlink() for path in ram_provenance):
        raise CoordinatorError("RAM evidence lacks immutable measurement provenance")
    for receipt, provenance_path in zip(ram_receipts, ram_provenance, strict=True):
        provenance = read_json(provenance_path)
        row, content = validate_ram_measurement_provenance(
            provenance, receipt, context=context, require_live_task=True,
        )
        if receipt.read_bytes() != content:
            raise CoordinatorError("RAM receipt bytes differ from exact measurement provenance")
        fields, rows = ENGINE.read_tsv(receipt)
        if fields != RAM_FIELDS or rows != [row]:
            raise CoordinatorError("RAM receipt row differs from exact measurement provenance")
    cleanup_candidates = collect_cleanup_candidates(paths["base"])
    reclaim_accounting = cleanup_reclaim_accounting(cleanup_candidates)
    return {
        "schema_version": REPRODUCIBILITY_SCHEMA,
        "analysis_id": context["policy"]["analysis_id"],
        "pair_id": context["pair_id"],
        "pair_order": list(PAIR_ORDER),
        "run_fingerprint": context["fingerprint"],
        "sequential_contract_sha256": stable_identity(ROOT / EXECUTION_CONTRACT)["sha256"],
        "terminal_gate_lock_sha256": stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK)["sha256"],
        "coordinator_sha256": coordinator_sha256(),
        "deep_validation": {
            "complete_aligned_family": True,
            "all_numerical_failures_retained_with_p_one": True,
            "bh_validated_once_over_complete_pair_family": True,
            "global_nuisance_all_valid_variants": True,
            "fresh_global_nuisance_in_matched_benchmark": True,
            "successful_full_scan_ram_evidence": True,
            "canonical_locus_claim_allowed": False,
            "counts": validated["counts"],
        },
        "canonical": {
            "ledger": validated["canonical_ledger"],
            "provenance": validated["canonical_provenance"],
        },
        "embedded_reproducibility_artifacts": embedded,
        "ram_receipts": [stable_identity(path) for path in ram_receipts],
        "ram_measurement_provenance": [stable_identity(path) for path in ram_provenance],
        "full_scan_ram_evidence": {
            "phase": "full_scan",
            "analysis": FULL_SCAN_ANALYSIS,
            "peak_process_tree_rss_bytes": int(full_scan_ram["peak_process_tree_rss_bytes"]),
            "rss_measurement_method": full_scan_ram["rss_measurement_method"],
            "output_hash": full_scan_ram["output_hash"],
            "receipt": stable_identity(ram_receipt_path(context, "full_scan")),
            "provenance": stable_identity(
                ram_receipt_path(context, "full_scan").with_suffix(".provenance.json")
            ),
        },
        "cleanup_candidates": cleanup_candidates,
        "cleanup_reclaim_accounting": reclaim_accounting,
        "cleanup_rule": "ONLY_LISTED_EXACT_INODES;ABSENT_IS_RESUMED_SUCCESS;IDENTITY_DRIFT_FAILS_CLOSED",
        "regeneration": "FROZEN_DENSE_SOURCES_PLUS_GATE_CODE_AND_EMBEDDED_MATERIALIZATION_TASK_NUISANCE_BENCHMARK_SUMMARY",
        "source_or_canonical_deletion": False,
    }


def validate_cleanup_bundle(context: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    bundle_path = cleanup_bundle_path(context)
    bundle = read_json(bundle_path)
    if (
        bundle.get("schema_version") != REPRODUCIBILITY_SCHEMA
        or bundle.get("analysis_id") != context["policy"]["analysis_id"]
        or bundle.get("pair_id") != context["pair_id"]
        or bundle.get("pair_order") != list(PAIR_ORDER)
        or bundle.get("run_fingerprint") != context["fingerprint"]
        or bundle.get("sequential_contract_sha256") != stable_identity(ROOT / EXECUTION_CONTRACT)["sha256"]
        or bundle.get("terminal_gate_lock_sha256") != stable_identity(ROOT / GATE.TERMINAL_GATE_LOCK)["sha256"]
        or bundle.get("coordinator_sha256") != coordinator_sha256()
        or bundle.get("source_or_canonical_deletion") is not False
        or contract.get("schema_version") != CONTRACT_SCHEMA
    ):
        raise CoordinatorError("immutable PLACO+ cleanup bundle identity drifted")
    canonical = bundle.get("canonical", {})
    expected_canonical_paths = {
        "ledger": ENGINE.materialized_paths(context)["canonical_ledger"],
        "provenance": ENGINE.materialized_paths(context)["canonical_provenance"],
    }
    for key in ("ledger", "provenance"):
        record = canonical.get(key)
        if not isinstance(record, dict):
            raise CoordinatorError("cleanup bundle lacks canonical evidence")
        path = ROOT / str(record.get("path", ""))
        if path != expected_canonical_paths[key]:
            raise CoordinatorError("cleanup bundle points to the wrong canonical PLACO+ artifact")
        if not identity_matches(path, record):
            raise CoordinatorError("canonical PLACO+ evidence changed after cleanup authorization")
    deep_validation = bundle.get("deep_validation")
    if (
        not isinstance(deep_validation, dict)
        or deep_validation.get("complete_aligned_family") is not True
        or deep_validation.get("all_numerical_failures_retained_with_p_one") is not True
        or deep_validation.get("bh_validated_once_over_complete_pair_family") is not True
        or deep_validation.get("global_nuisance_all_valid_variants") is not True
        or deep_validation.get("fresh_global_nuisance_in_matched_benchmark") is not True
        or deep_validation.get("successful_full_scan_ram_evidence") is not True
        or deep_validation.get("canonical_locus_claim_allowed") is not False
        or not isinstance(deep_validation.get("counts"), dict)
    ):
        raise CoordinatorError("cleanup bundle lost its exact scientific-validation attestation")
    canonical_provenance = read_json(expected_canonical_paths["provenance"])
    if deep_validation["counts"] != canonical_provenance.get("complete_family_counts"):
        raise CoordinatorError("cleanup validation counts differ from canonical PLACO+ provenance")
    embedded = bundle.get("embedded_reproducibility_artifacts")
    if not isinstance(embedded, dict) or not PRESERVED_WORK_ARTIFACTS.issubset(embedded):
        raise CoordinatorError("cleanup bundle lacks required reproducibility artifacts")
    embedded_content: dict[str, bytes] = {}
    for relative, record in embedded.items():
        if relative not in PRESERVED_WORK_ARTIFACTS | OPTIONAL_PRESERVED_WORK_ARTIFACTS:
            raise CoordinatorError("cleanup bundle contains an unexpected embedded artifact")
        try:
            content = base64.b64decode(record["base64"], validate=True)
        except (KeyError, TypeError, ValueError) as error:
            raise CoordinatorError("cleanup bundle contains invalid embedded bytes") from error
        if len(content) != record.get("bytes") or hashlib.sha256(content).hexdigest() != record.get("sha256"):
            raise CoordinatorError("cleanup bundle embedded artifact hash drifted")
        embedded_content[relative] = content
    try:
        benchmark_reader = csv.DictReader(
            io.StringIO(embedded_content["benchmark.tsv"].decode("utf-8")), delimiter="\t",
        )
        benchmark_fields = list(benchmark_reader.fieldnames or [])
        benchmark_rows = list(benchmark_reader)
    except (KeyError, UnicodeError, csv.Error) as error:
        raise CoordinatorError("cleanup bundle embedded benchmark is unreadable") from error
    if (
        benchmark_fields != ENGINE.BENCHMARK_RAW_FIELDS + [
            "peak_process_tree_rss_bytes", "rss_measurement_method", "wrapper_wall_seconds",
            "rss_sample_interval_seconds", "host_physical_memory_bytes", "host_free_bytes_at_benchmark",
            "peak_rss_fraction_of_physical_memory",
        ]
        or len(benchmark_rows) != 1
        or benchmark_rows[0].get("global_nuisance_checkpoint_reused") != "FALSE"
        or benchmark_rows[0].get("workers") != str(WORKERS)
    ):
        raise CoordinatorError("cleanup bundle benchmark did not freshly measure global nuisance")
    embedded_benchmark = benchmark_rows[0]
    ram_records_by_field: dict[str, dict[str, dict[str, Any]]] = {}
    for field in ("ram_receipts", "ram_measurement_provenance"):
        records = bundle.get(field)
        if not isinstance(records, list) or not records:
            raise CoordinatorError("cleanup bundle lacks complete RAM measurement evidence")
        indexed: dict[str, dict[str, Any]] = {}
        for record in records:
            if not isinstance(record, dict) or not identity_matches(ROOT / str(record.get("path", "")), record):
                raise CoordinatorError("RAM measurement evidence changed after cleanup authorization")
            record_path = str(record["path"])
            if record_path in indexed:
                raise CoordinatorError("cleanup bundle repeats RAM measurement evidence")
            indexed[record_path] = record
        ram_records_by_field[field] = indexed
    provenance_records = bundle["ram_measurement_provenance"]
    provenance_by_phase: dict[str, dict[str, Any]] = {}
    for provenance_record in provenance_records:
        provenance_path = ROOT / str(provenance_record["path"])
        suffix = ".provenance.json"
        if not provenance_path.name.endswith(suffix):
            raise CoordinatorError("cleanup bundle has a malformed RAM provenance path")
        receipt = provenance_path.with_name(
            provenance_path.name.removesuffix(suffix) + ".tsv"
        )
        receipt_relative = str(receipt.relative_to(ROOT))
        if receipt_relative not in ram_records_by_field["ram_receipts"]:
            raise CoordinatorError("cleanup bundle RAM provenance lacks its bound TSV receipt")
        provenance = read_json(provenance_path)
        row, content = validate_ram_measurement_provenance(
            provenance, receipt, context=context, require_live_task=False,
        )
        if receipt.read_bytes() != content:
            raise CoordinatorError("cleanup bundle RAM receipt bytes drifted")
        fields, rows = ENGINE.read_tsv(receipt)
        if fields != RAM_FIELDS or rows != [row]:
            raise CoordinatorError("cleanup bundle RAM receipt row drifted")
        phase = str(provenance["phase"])
        if phase in provenance_by_phase:
            raise CoordinatorError("cleanup bundle repeats one RAM measurement phase")
        provenance_by_phase[phase] = provenance
    expected_provenance_paths = {
        str((ROOT / relative).with_suffix(".provenance.json").relative_to(ROOT))
        for relative in ram_records_by_field["ram_receipts"]
    }
    if set(ram_records_by_field["ram_measurement_provenance"]) != expected_provenance_paths:
        raise CoordinatorError("cleanup bundle RAM receipt/provenance families are not exact")
    rejected = sorted(
        phase for phase in provenance_by_phase if phase.startswith(FULL_SCAN_REJECTION_PREFIX)
    )
    if rejected:
        raise CoordinatorError("cleanup bundle contains prior full-scan RAM rejection evidence")
    full_scan = provenance_by_phase.get("full_scan")
    matched_benchmark = provenance_by_phase.get("matched_worker_benchmark")
    attestation = bundle.get("full_scan_ram_evidence")
    if not isinstance(full_scan, dict) or not isinstance(attestation, dict):
        raise CoordinatorError("cleanup bundle lacks exact successful full-scan RAM evidence")
    full_receipt = ram_receipt_path(context, "full_scan")
    expected_attestation = {
        "phase": "full_scan",
        "analysis": FULL_SCAN_ANALYSIS,
        "peak_process_tree_rss_bytes": int(full_scan["peak_process_tree_rss_bytes"]),
        "rss_measurement_method": full_scan["rss_measurement_method"],
        "output_hash": full_scan["output_hash"],
        "receipt": stable_identity(full_receipt),
        "provenance": stable_identity(full_receipt.with_suffix(".provenance.json")),
    }
    if (
        attestation != expected_attestation
        or full_scan.get("row", {}).get("analysis") != FULL_SCAN_ANALYSIS
        or int(full_scan.get("exit_status", -1)) != 0
        or full_scan.get("output_hash") != canonical["ledger"]["sha256"]
    ):
        raise CoordinatorError("cleanup full-scan RAM evidence differs from the canonical ledger")
    try:
        reserve = int(contract["ram_safety_reserve_bytes"])
        limit = int(contract["ram_limit_bytes"])
    except (KeyError, TypeError, ValueError) as error:
        raise CoordinatorError("sequential contract lacks the cleanup RAM envelope") from error
    if limit != RAM_LIMIT_BYTES:
        raise CoordinatorError("sequential cleanup contract RAM limit drifted")
    admit_ram(
        {
            "workers": str(full_scan["workers"]),
            "peak_process_tree_rss_bytes": str(full_scan["peak_process_tree_rss_bytes"]),
        },
        {"resource_envelope": {"ram_safety_reserve_bytes": reserve}},
        physical_memory_bytes=limit,
    )
    base = ENGINE.materialized_paths(context)["base"]
    benchmark_record = embedded["benchmark.tsv"]
    nuisance_record = embedded["checkpoints/nuisance.rds"]
    expected_benchmark_hash = hashlib.sha256(canonical_json_bytes({
        "benchmark": {
            "path": str((base / "benchmark.tsv").relative_to(ROOT)),
            "bytes": benchmark_record["bytes"],
            "sha256": benchmark_record["sha256"],
        },
        "nuisance": {
            "path": str((base / "checkpoints/nuisance.rds").relative_to(ROOT)),
            "bytes": nuisance_record["bytes"],
            "sha256": nuisance_record["sha256"],
        },
        "input_sha256": embedded_benchmark["input_sha256"],
        "workers": embedded_benchmark["workers"],
    })).hexdigest()
    if (
        not isinstance(matched_benchmark, dict)
        or matched_benchmark.get("row", {}).get("analysis") != MATCHED_BENCHMARK_ANALYSIS
        or int(matched_benchmark.get("exit_status", -1)) != 0
        or int(matched_benchmark.get("peak_process_tree_rss_bytes", -1))
        != int(embedded_benchmark["peak_process_tree_rss_bytes"])
        or not math.isclose(
            float(matched_benchmark.get("runtime_seconds", -1)),
            float(embedded_benchmark["wrapper_wall_seconds"]), rel_tol=0, abs_tol=0,
        )
        or matched_benchmark.get("output_hash") != expected_benchmark_hash
    ):
        raise CoordinatorError("cleanup bundle lacks exact matched benchmark RAM evidence")
    admit_ram(
        {
            "workers": str(matched_benchmark["workers"]),
            "peak_process_tree_rss_bytes": str(
                matched_benchmark["peak_process_tree_rss_bytes"]
            ),
        },
        {"resource_envelope": {"ram_safety_reserve_bytes": reserve}},
        physical_memory_bytes=limit,
    )
    candidates = bundle.get("cleanup_candidates")
    if not isinstance(candidates, list) or not candidates:
        raise CoordinatorError("cleanup bundle lacks exact work candidates")
    base = ENGINE.materialized_paths(context)["base"]
    indexed_candidates: dict[str, dict[str, Any]] = {}
    for record in candidates:
        if not isinstance(record, dict):
            raise CoordinatorError("cleanup bundle has a malformed candidate")
        path = ROOT / str(record.get("path", ""))
        try:
            relative = path.relative_to(base).as_posix()
        except ValueError as error:
            raise CoordinatorError("cleanup bundle candidate escapes the pair work root") from error
        if not allowed_cleanup_relative(relative):
            raise CoordinatorError("cleanup bundle candidate is not explicitly regenerable")
        if relative in indexed_candidates:
            raise CoordinatorError("cleanup bundle repeats a work candidate")
        indexed_candidates[relative] = record
        if path.exists() or path.is_symlink():
            if not identity_matches(
                path, record, include_inode=True,
                allow_empty=KNOWN_TEMP_RE.fullmatch(relative) is not None,
            ):
                raise CoordinatorError("cleanup candidate changed after immutable authorization")
    validate_cleanup_reclaim_accounting(candidates, bundle.get("cleanup_reclaim_accounting"))
    if base.exists():
        for path in sorted(base.rglob("*")):
            observed = os.lstat(path)
            if stat.S_ISDIR(observed.st_mode):
                relative_directory = path.relative_to(base).as_posix()
                if relative_directory not in KNOWN_WORK_DIRECTORIES:
                    raise CoordinatorError(
                        f"unlisted work directory appeared after cleanup authorization: {relative_directory}"
                    )
                continue
            relative = path.relative_to(base).as_posix()
            if relative not in indexed_candidates:
                raise CoordinatorError(
                    f"unlisted work artifact appeared after cleanup authorization: {relative}"
                )
            if not stat.S_ISREG(observed.st_mode):
                raise CoordinatorError(f"cleanup work artifact is not regular: {relative}")
    return bundle


def complete_cleanup(context: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    bundle = validate_cleanup_bundle(context, contract)
    for record in bundle["cleanup_candidates"]:
        path = ROOT / record["path"]
        if not (path.exists() or path.is_symlink()):
            continue
        safe_repo_path(path, "cleanup candidate")
        current = os.lstat(path)
        if (
            not stat.S_ISREG(current.st_mode)
            or (current.st_dev, current.st_ino, current.st_size)
            != (record["device"], record["inode"], record["bytes"])
        ):
            raise CoordinatorError("cleanup candidate identity changed immediately before unlink")
        path.unlink()
        fsync_directory(path.parent)
    base = ENGINE.materialized_paths(context)["base"]
    if base.exists():
        remaining_files = [
            path for path in base.rglob("*")
            if not stat.S_ISDIR(os.lstat(path).st_mode)
        ]
        if remaining_files:
            raise CoordinatorError(
                f"unlisted PLACO+ work artifact remains after cleanup: {remaining_files[0]}"
            )
    if base.exists():
        remaining_directories = sorted(
            path.relative_to(base).as_posix()
            for path in base.rglob("*")
            if stat.S_ISDIR(os.lstat(path).st_mode)
        )
        unknown_directories = [
            relative for relative in remaining_directories
            if relative not in KNOWN_WORK_DIRECTORIES
        ]
        if unknown_directories:
            raise CoordinatorError(
                f"unlisted PLACO+ work directory remains after cleanup: {unknown_directories[0]}"
            )
        for relative in sorted(KNOWN_WORK_DIRECTORIES, reverse=True):
            directory = base / relative
            try:
                directory.rmdir()
            except FileNotFoundError:
                pass
            except OSError as error:
                raise CoordinatorError(
                    f"authorized PLACO+ work directory could not be retired: {relative}: {error}"
                ) from error
        try:
            base.rmdir()
        except OSError as error:
            raise CoordinatorError(
                f"PLACO+ work root changed before retirement: {base}: {error}"
            ) from error
    if base.exists() or base.is_symlink():
        raise CoordinatorError("PLACO+ work root still exists after authorized cleanup")
    for record in bundle["cleanup_candidates"]:
        path = ROOT / record["path"]
        if path.exists() or path.is_symlink():
            raise CoordinatorError("authorized PLACO+ work cleanup is incomplete")
    completion = {
        "schema_version": CLEANUP_COMPLETE_SCHEMA,
        "analysis_id": context["policy"]["analysis_id"],
        "pair_id": context["pair_id"],
        "run_fingerprint": context["fingerprint"],
        "cleanup_bundle_sha256": stable_identity(cleanup_bundle_path(context))["sha256"],
        "removed_artifact_count": len(bundle["cleanup_candidates"]),
        "removed_logical_bytes": sum(int(item["bytes"]) for item in bundle["cleanup_candidates"]),
        "reclaim_accounting_at_authorization": bundle["cleanup_reclaim_accounting"],
        "canonical_ledger_retained": bundle["canonical"]["ledger"],
        "canonical_provenance_retained": bundle["canonical"]["provenance"],
        "source_or_canonical_deletion": False,
        "ld_extracted": False,
        "locus_claim_status": "BLOCKED_PENDING_FULL_P_FAMILY_AND_SEPARATE_LD_COLLATOR",
        "coordinator_sha256": coordinator_sha256(),
    }
    path = cleanup_complete_path(context)
    content = json.dumps(completion, indent=2, sort_keys=True).encode() + b"\n"
    publish_bytes_no_replace(path, content)
    if read_json(path) != completion:
        raise CoordinatorError("cleanup completion receipt drifted")
    return completion


def validate_clean_completed_pair(context: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    bundle = validate_cleanup_bundle(context, contract)
    completion = read_json(cleanup_complete_path(context))
    expected = {
        "schema_version": CLEANUP_COMPLETE_SCHEMA,
        "analysis_id": context["policy"]["analysis_id"],
        "pair_id": context["pair_id"],
        "run_fingerprint": context["fingerprint"],
        "cleanup_bundle_sha256": stable_identity(cleanup_bundle_path(context))["sha256"],
        "removed_artifact_count": len(bundle["cleanup_candidates"]),
        "removed_logical_bytes": sum(int(item["bytes"]) for item in bundle["cleanup_candidates"]),
        "reclaim_accounting_at_authorization": bundle["cleanup_reclaim_accounting"],
        "canonical_ledger_retained": bundle["canonical"]["ledger"],
        "canonical_provenance_retained": bundle["canonical"]["provenance"],
        "source_or_canonical_deletion": False,
        "ld_extracted": False,
        "locus_claim_status": "BLOCKED_PENDING_FULL_P_FAMILY_AND_SEPARATE_LD_COLLATOR",
        "coordinator_sha256": coordinator_sha256(),
    }
    if completion != expected:
        raise CoordinatorError("cleanup completion receipt is missing or drifted")
    for record in bundle["cleanup_candidates"]:
        path = ROOT / record["path"]
        if path.exists() or path.is_symlink():
            raise CoordinatorError("cleanup-complete receipt coexists with a retired work artifact")
    rebuild_ram_aggregate()
    return completion


def _process_pair_under_worker_lease(
    context: dict[str, Any], gate_lock: dict[str, Any], contract: dict[str, Any],
) -> str:
    require_active_worker_lease(context)
    repair_record_publication_temporaries(context)
    reconcile_failure_events(context)
    cleanup_plan = cleanup_bundle_path(context)
    cleanup_done = cleanup_complete_path(context)
    if cleanup_done.exists() or cleanup_done.is_symlink():
        validate_clean_completed_pair(context, contract)
        return "RESUMED_COMPLETE"
    if cleanup_plan.exists() or cleanup_plan.is_symlink():
        try:
            complete_cleanup(context, contract)
        except CoordinatorError as error:
            raise cleanup_failure(error, context, gate_lock) from error
        return "RESUMED_CLEANUP_COMPLETE"
    repair_crash_temporaries(context)
    paths = ENGINE.materialized_paths(context)
    state = reconcile_partial_canonical_publication(context, gate_lock, contract)
    if state == "COMPLETE":
        benchmark = validate_admit_fresh_benchmark(
            context, gate_lock, require_ram_receipt=True,
        )
        validated = validate_completed_pair(context, gate_lock, benchmark)
    else:
        ensure_materialized(context)
        repair_partial_execution_families(context)
        benchmark = ensure_benchmark(context, gate_lock)
        run_full_scan_monitored(context, gate_lock, benchmark)
        validated = validate_completed_pair(context, gate_lock, benchmark)
    try:
        bundle = build_cleanup_bundle(context, validated, contract)
        publish_bytes_no_replace(
            cleanup_plan, json.dumps(bundle, indent=2, sort_keys=True).encode() + b"\n",
        )
        complete_cleanup(context, contract)
    except CoordinatorError as error:
        raise cleanup_failure(error, context, gate_lock) from error
    return "COMPLETE_AND_CLEANED"


def process_pair(context: dict[str, Any], gate_lock: dict[str, Any], contract: dict[str, Any]) -> str:
    with worker_lease(context):
        return _process_pair_under_worker_lease(context, gate_lock, contract)


def _preflight_rows_snapshot(
    contract: dict[str, Any], configured: list[tuple[dict[str, Any], dict[str, Any]]],
) -> list[dict[str, object]]:
    free = shutil.disk_usage(ROOT).free
    canonical_bytes_by_pair: dict[str, int] = {}
    seen_canonical_inodes: set[tuple[int, int]] = set()
    for context, _ in configured:
        total = 0
        paths = ENGINE.materialized_paths(context)
        for path in (paths["canonical_ledger"], paths["canonical_provenance"]):
            if not (path.exists() or path.is_symlink()):
                continue
            record = stable_identity(path)
            inode = int(record["device"]), int(record["inode"])
            if inode not in seen_canonical_inodes:
                total += int(record["bytes"])
                seen_canonical_inodes.add(inode)
        canonical_bytes_by_pair[context["pair_id"]] = total
    rows: list[dict[str, object]] = []
    for index, (context, gate_lock) in enumerate(configured):
        pair = context["pair_id"]
        paths = ENGINE.materialized_paths(context)
        state = canonical_state(paths)
        cleanup_done = artifact_present(cleanup_complete_path(context))
        forecast, storage_resume_state = pair_preflight_storage_forecast(context, free)
        prior_reserves = [
            (prior_context["pair_id"], pending_pair_canonical_reserve_bytes(prior_context))
            for prior_context, _ in configured[:index]
        ]
        pending_prior_pairs = [
            prior_pair for prior_pair, reserve in prior_reserves if reserve > 0
        ]
        actual_prior_retained = sum(
            canonical_bytes_by_pair[prior_context["pair_id"]]
            for prior_context, _ in configured[:index]
        )
        pending_prior_reserve = sum(reserve for _, reserve in prior_reserves)
        projected_free = max(0, free - pending_prior_reserve)
        if cleanup_done:
            sequential_status = "RESUMED_COMPLETE_NO_PAIR_WORK_STORAGE_REQUIRED"
        else:
            if forecast["additional_bytes_required"]:
                sequential_status = "BLOCKED_CURRENT_STORAGE"
            elif pending_prior_pairs and projected_free < int(forecast["required_bytes"]):
                sequential_status = "CONDITIONAL_PROJECTED_STORAGE_SHORTFALL_RECHECK_ACTUAL_OUTPUTS"
            elif pending_prior_pairs:
                sequential_status = "CONDITIONAL_ON_PENDING_PRIOR_CANONICAL_LEDGER_BYTES"
            else:
                sequential_status = "PASS_CURRENT_STORAGE_WITH_ACTUAL_RETAINED_CANONICALS"
        rows.append({
            "pair": pair,
            "order": index + 1,
            "source_rows": sum(int(source["rows"]) for source in context["sources"]),
            **forecast,
            "actual_retained_canonical_bytes_all_pairs": sum(canonical_bytes_by_pair.values()),
            "actual_prior_retained_canonical_bytes": actual_prior_retained,
            "pending_prior_pairs": ",".join(pending_prior_pairs) if pending_prior_pairs else "NONE",
            "pending_prior_canonical_planning_reserve_bytes": pending_prior_reserve,
            "projected_free_after_pending_prior_planning_reserve": projected_free,
            "projected_additional_storage_required": max(
                0, int(forecast["required_bytes"]) - projected_free,
            ),
            "planning_reserve_rule": (
                "600_BYTES_PER_PENDING_PRIOR_PAIR_SOURCE_ROW_ONLY_WHEN_CANONICAL_OR_STAGED_LEDGER_BYTES_"
                "ARE_NOT_ALREADY_ALLOCATED;DYNAMIC_ACTUAL_FREE_BYTES_RECHECK_IS_AUTHORITATIVE"
            ),
            "maximum_cumulative_pending_prior_canonical_bytes_before_admission": (
                max(0, free - int(forecast["required_bytes"])) if pending_prior_pairs else 0
            ),
            "sequential_storage_status": sequential_status,
            "ram_limit_bytes": contract["ram_limit_bytes"],
            "ram_reserve_bytes": contract["ram_safety_reserve_bytes"],
            "workers": WORKERS,
            "state": state,
            "storage_resume_state": storage_resume_state,
            "cleanup_complete": cleanup_done,
        })
        if gate_lock["reference_and_publication"]["ld_locus_publication_allowed"] is not False:
            raise CoordinatorError("terminal gate unexpectedly allows LD/locus publication")
    return rows


def preflight_rows() -> list[dict[str, object]]:
    """Take one race-free storage snapshot, including orphan-worker exclusion."""
    with BRIDGE.execution_lock("sequential-full-p-family", "A,B,CONTROL"):
        contract = verify_execution_contract()
        configured: list[tuple[dict[str, Any], dict[str, Any]]] = []
        for pair in PAIR_ORDER:
            try:
                configured.append(BRIDGE.configure_context(pair))
            except BRIDGE.BridgeError as error:
                raise CoordinatorError(str(error)) from error
        for context, _ in configured:
            with worker_lease(context):
                pass
        return _preflight_rows_snapshot(contract, configured)


def run_sequential() -> list[dict[str, str]]:
    outcomes: list[dict[str, str]] = []
    with BRIDGE.execution_lock("sequential-full-p-family", "A,B,CONTROL"):
        contract = verify_execution_contract()
        configured: list[tuple[dict[str, Any], dict[str, Any]]] = []
        for pair in PAIR_ORDER:
            try:
                configured.append(BRIDGE.configure_context(pair))
            except BRIDGE.BridgeError as error:
                raise CoordinatorError(str(error)) from error
        for context, _ in configured:
            with worker_lease(context):
                repair_record_publication_temporaries(context)
        # Only after every pair's no-replace temps are repaired may a failure
        # event derive RAM evidence, because record_ram rebuilds the global table.
        for context, _ in configured:
            with worker_lease(context):
                reconcile_failure_events(context)
        for context, gate_lock in configured:
            outcome = process_pair(context, gate_lock, contract)
            outcomes.append({"pair": context["pair_id"], "outcome": outcome})
    return outcomes


def verify_family() -> list[dict[str, str]]:
    outcomes: list[dict[str, str]] = []
    with BRIDGE.execution_lock("sequential-full-p-family", "A,B,CONTROL"):
        contract = verify_execution_contract()
        configured: list[tuple[dict[str, Any], dict[str, Any]]] = []
        for pair in PAIR_ORDER:
            try:
                configured.append(BRIDGE.configure_context(pair))
            except BRIDGE.BridgeError as error:
                raise CoordinatorError(str(error)) from error
        for context, _ in configured:
            with worker_lease(context):
                repair_record_publication_temporaries(context)
        for context, _ in configured:
            with worker_lease(context):
                reconcile_failure_events(context)
        for context, _ in configured:
            with worker_lease(context):
                validate_clean_completed_pair(context, contract)
            outcomes.append({"pair": context["pair_id"], "outcome": "VERIFIED_COMPLETE_CLEAN"})
    return outcomes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--seal-contract", action="store_true")
    action.add_argument("--preflight", action="store_true")
    action.add_argument("--run", action="store_true")
    action.add_argument("--verify", action="store_true")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.seal_contract:
        payload = seal_execution_contract()
        print(
            "TRACK_B_PLACO_SEQUENTIAL_CONTRACT_SEALED "
            f"workers={payload['workers']} ram_limit={payload['ram_limit_bytes']} order=A,B,CONTROL"
        )
    elif args.preflight:
        print(json.dumps({"status": "PREFLIGHT", "pairs": preflight_rows()}, indent=2, sort_keys=True))
    elif args.verify:
        print(json.dumps({"status": "VERIFIED", "pairs": verify_family()}, indent=2, sort_keys=True))
    else:
        if not args.execute:
            raise CoordinatorError("explicit --execute is required for sequential PLACO+ production")
        print(json.dumps({"status": "COMPLETE", "pairs": run_sequential()}, indent=2, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except CoordinatorError as error:
        raise SystemExit(f"ERROR: {error}") from error
