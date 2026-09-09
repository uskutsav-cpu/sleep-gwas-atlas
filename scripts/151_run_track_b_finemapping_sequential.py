#!/usr/bin/env python3
"""Sequential, restart-safe Track B SuSiE-RSS and trait-coloc coordinator.

This additive coordinator consumes only the fully materialized manifest sealed
by script 148.  One pair x official LAVA block is launched in one fresh OS
session at a time.  It never changes the frozen locus, SNP family, signed LD,
model parameters, or priors in response to a result or resource failure.

The public verification API is deliberately side-effect free:

* ``verify_materialization_measurement``
* ``verify_sealed_run``
* ``validate_canonical_family``

No work occurs without ``--execute``; ``--verify`` only reads and validates.
Production ``--execute`` must run outside a sandbox that denies ``ps`` process
inventory, because successful parent-plus-PGID RAM measurement fails closed
without that OS permission.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import ctypes
import ctypes.util
import datetime as dt
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
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
RUNTIME_REL = Path("scripts/148_track_b_finemapping_runtime.py")
R_ENGINE_REL = Path("scripts/150_run_track_b_susie_coloc.R")
RSCRIPT_REL = Path(".r-env/bin/Rscript")
RUN_ROOT_REL = Path("results/track_b/finemapping/runs")
MATERIALIZATION_METRICS_REL = Path("results/track_b/finemapping/materialization_metrics")
MATERIALIZATION_ATTEMPTS_REL = Path("results/track_b/finemapping/materialization_attempts")
STAGE_ROOT_REL = Path("work/track_b_finemapping/runs")
QUARANTINE_ROOT_REL = Path("work/track_b_finemapping/quarantine")
LEASE_ROOT_REL = Path("work/track_b_finemapping/leases")
PUBLICATION_STAGE_REL = Path("work/track_b_finemapping/publications")
PUBLICATION_ROOT_REL = Path("results/track_b/finemapping/publications")
RESULT_PROVENANCE_REL = Path("results/track_b/finemapping/results.provenance.json")
RUN_INDEX_REL = Path("results/track_b/finemapping/run_index.tsv")
RESOURCE_METRICS_REL = Path("results/track_b/finemapping/resource_metrics.tsv")
COMPONENT_NAMESPACE_REL = Path("results/track_b/finemapping/RAM_BY_LOCUS.namespace.json")

RUN_SCHEMA = "sleep-atlas-track-b-finemapping-sealed-locus-run.2"
MEASUREMENT_SCHEMA = "sleep-atlas-track-b-finemapping-materialization-measurement.2"
PUBLICATION_SCHEMA = "sleep-atlas-track-b-finemapping-canonical-publication.1"
RAM_AMENDMENT_SCHEMA = "sleep-atlas-track-b-finemapping-pgid-ram-amendment.1"
ANALYSIS_ID = "track-b-v1.0-finemapping-trait-coloc"
PGID_BACKEND = "PARENT_PLUS_PGID_AGGREGATE_PS_RSS_BYTES"
RUSAGE_FAILURE_BACKEND = "RUSAGE_MAXRSS_PLATFORM_NORMALIZED_TO_BYTES"
PHYSICAL_ENVELOPE_BYTES = 8 * 1024**3
RESERVE_BYTES = 1024**3
POLL_SECONDS = 0.025
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")

RUN_INDEX_FIELDS = [
    "analysis_id", "run_fingerprint", "locus_entry_id", "pair_id", "CHR",
    "ld_block_id", "variant_count", "terminal_status", "run_directory",
    "run_seal_sha256", "output_manifest_sha256", "error",
]
COMPONENT_RAM_FIELDS = [
    "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
    "runtime_sec", "exit_status", "output_hash",
]


class CoordinatorError(RuntimeError):
    """Fail-closed orchestration, measurement, or publication error."""


def _load_runtime(root: Path = ROOT):
    path = root / RUNTIME_REL
    specification = importlib.util.spec_from_file_location("track_b_finemap_runtime_148", path)
    if specification is None or specification.loader is None:
        raise CoordinatorError(f"cannot load runtime: {path}")
    module = importlib.util.module_from_spec(specification)
    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        specification.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


def _utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def _canonical_json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical_json(value)).hexdigest()


def _normalize_ru_maxrss(value: int) -> int:
    if value < 0:
        raise CoordinatorError("negative ru_maxrss")
    return value if sys.platform == "darwin" else value * 1024


def ram_measurement_amendment(runtime: Any, root: Path = ROOT) -> dict[str, Any]:
    """Bind the user's additive PGID requirement without changing policy bytes."""

    _, policy_identity = runtime.load_policy(root)
    return {
        "schema_version": RAM_AMENDMENT_SCHEMA,
        "analysis_id": ANALYSIS_ID,
        "supersedes_only": "ram_aware_execution_contract.allowed_rss_measurement_backends",
        "successful_run_backend": PGID_BACKEND,
        "success_requirement": (
            "POSITIVE_LIVE_DEDUPLICATED_SUM_FROM_ONE_PS_SNAPSHOT_OF_COORDINATOR_"
            "PARENT_PLUS_EVERY_PROCESS_IN_OWNED_FRESH_SESSION_PGID"
        ),
        "failure_only_fallback": RUSAGE_FAILURE_BACKEND,
        "scientific_model_changed": False,
        "input_family_changed": False,
        "classification_thresholds_changed": False,
        "policy_identity": policy_identity,
        "sampler_script": runtime.stable_identity(root, SCRIPT),
    }


def process_group_rss_bytes(
    process_group_id: int, *, coordinator_process_id: int | None = None,
) -> int:
    """Sum coordinator RSS plus every member of the owned fresh-session PGID.

    Membership is selected from one ``ps`` snapshot and PIDs are deduplicated,
    so a coordinator that happens to share the PGID can never be counted twice.
    PGID rather than PPID membership retains descendants that have reparented.
    """

    if process_group_id <= 0:
        raise CoordinatorError("invalid owned process-group ID")
    completed = subprocess.run(
        ["ps", "-axo", "pid=,pgid=,rss="], capture_output=True, text=True, check=False,
    )
    if completed.returncode != 0:
        raise CoordinatorError(f"PGID RSS sampler failed: {completed.stderr.strip()}")
    total_kib = 0
    included: set[int] = set()
    coordinator = os.getpid() if coordinator_process_id is None else int(coordinator_process_id)
    if coordinator <= 0:
        raise CoordinatorError("invalid coordinator process ID")
    for line in completed.stdout.splitlines():
        values = line.split()
        if len(values) != 3:
            continue
        try:
            pid, pgid, rss_kib = (int(value) for value in values)
        except ValueError:
            continue
        if pid == coordinator or pgid == process_group_id:
            if rss_kib < 0:
                raise CoordinatorError("PGID sampler returned negative RSS")
            if pid not in included:
                total_kib += rss_kib
                included.add(pid)
    # A valid sample must include both the persistent coordinator and at least
    # one member of the separately owned worker PGID. Verify both explicitly.
    coordinator_seen = coordinator in included
    worker_seen = any(
        len(line.split()) == 3
        and line.split()[0].lstrip("-").isdigit()
        and line.split()[1].lstrip("-").isdigit()
        and int(line.split()[1]) == process_group_id
        for line in completed.stdout.splitlines()
    )
    return total_kib * 1024 if coordinator_seen and worker_seen else 0


def _terminate_process_group(process: subprocess.Popen[Any]) -> tuple[bool, bool]:
    """Terminate all descendants and prove the fresh-session PGID disappeared."""

    existed = False
    for requested_signal, grace in ((signal.SIGTERM, 0.5), (signal.SIGKILL, 1.0)):
        try:
            os.killpg(process.pid, requested_signal)
            existed = True
        except ProcessLookupError:
            process.poll()
            return existed, True
        except PermissionError:
            return True, False
        deadline = time.monotonic() + grace
        while time.monotonic() < deadline:
            process.poll()
            time.sleep(0.025)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        process.poll()
        return True, True
    except PermissionError:
        return True, False
    return True, False


def monitored_command(
    command: Sequence[str], *, root: Path = ROOT, inherited_fds: Sequence[int] = (),
    maximum_owned_pgid_rss_bytes: int, minimum_host_reserve_bytes: int = RESERVE_BYTES,
    extra_environment: Mapping[str, str] | None = None,
    memory_available: Any | None = None,
) -> dict[str, Any]:
    """Launch one fresh session and measure/guard its aggregate PGID RSS.

    A successful receipt is impossible without at least one positive live PGID
    observation. RUSAGE is retained only as conservative evidence on a failed
    launch and can never authorize scientific publication.
    """

    if not command or maximum_owned_pgid_rss_bytes <= 0 or minimum_host_reserve_bytes < 0:
        raise CoordinatorError("invalid monitored-command resource contract")
    started_utc = _utc_now()
    started = time.monotonic()
    peak = 0
    samples = 0
    failure_reason: str | None = None
    process: subprocess.Popen[bytes] | None = None
    environment = os.environ.copy()
    if extra_environment:
        environment.update({str(key): str(value) for key, value in extra_environment.items()})
    before_rusage = _normalize_ru_maxrss(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    available_memory = memory_available or _load_runtime(root)._memory_available_bytes
    with tempfile.TemporaryFile(mode="w+b") as log_handle:
        try:
            process = subprocess.Popen(
                list(command), cwd=root, stdout=log_handle, stderr=subprocess.STDOUT,
                start_new_session=True, pass_fds=tuple(inherited_fds), env=environment,
            )
        except OSError as error:
            failure_reason = f"PROCESS_START_ERROR:{type(error).__name__}:{error}"
            return_code = 70
        if process is not None:
            while process.poll() is None:
                try:
                    observed = process_group_rss_bytes(process.pid)
                    available = int(available_memory())
                except BaseException as error:
                    failure_reason = f"PGID_SAMPLER_ERROR:{type(error).__name__}:{error}"
                    _, confirmed = _terminate_process_group(process)
                    if not confirmed:
                        failure_reason += ":PROCESS_GROUP_TERMINATION_UNCONFIRMED"
                    break
                if observed <= 0:
                    # A just-execing child can be briefly absent from one ps
                    # snapshot. Keep sampling while it is demonstrably live.
                    time.sleep(POLL_SECONDS)
                    continue
                samples += 1
                peak = max(peak, observed)
                if observed > maximum_owned_pgid_rss_bytes or available < minimum_host_reserve_bytes:
                    failure_reason = "PHYSICAL_RAM_CAP_OR_HOST_RESERVE_EXCEEDED"
                    _, confirmed = _terminate_process_group(process)
                    if not confirmed:
                        failure_reason += ":PROCESS_GROUP_TERMINATION_UNCONFIRMED"
                    break
                time.sleep(POLL_SECONDS)
            if failure_reason is None:
                group_existed, confirmed = _terminate_process_group(process)
                if group_existed:
                    failure_reason = "OWNED_PROCESS_GROUP_OUTLIVED_LEADER"
                    if not confirmed:
                        failure_reason += ":PROCESS_GROUP_TERMINATION_UNCONFIRMED"
            try:
                return_code = int(process.wait(timeout=2.0))
            except subprocess.TimeoutExpired:
                failure_reason = failure_reason or "PROCESS_LEADER_DID_NOT_REAP"
                _, confirmed = _terminate_process_group(process)
                if not confirmed:
                    failure_reason += ":PROCESS_GROUP_TERMINATION_UNCONFIRMED"
                return_code = 70
            if failure_reason and return_code == 0:
                return_code = 70
        elapsed = time.monotonic() - started
        log_handle.seek(0)
        log = log_handle.read()
    after_rusage = _normalize_ru_maxrss(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    success_measurement = return_code == 0 and failure_reason is None and samples > 0 and peak > 0
    if return_code == 0 and not success_measurement:
        return_code = 70
        failure_reason = failure_reason or "NO_POSITIVE_LIVE_PGID_RSS_SAMPLE"
    backend = PGID_BACKEND if samples > 0 else RUSAGE_FAILURE_BACKEND
    if samples == 0:
        peak = max(0, after_rusage, before_rusage)
    return {
        "command": list(command), "started_utc": started_utc, "finished_utc": _utc_now(),
        "elapsed_monotonic_seconds": elapsed, "peak_rss_bytes": int(peak),
        "rss_measurement_backend": backend, "positive_pgid_samples": samples,
        "exit_code": return_code, "guard_failure": failure_reason,
        "success_measurement": success_measurement, "log": log,
    }


@contextlib.contextmanager
def worker_lease(
    runtime: Any, root: Path, pair_id: str, run_fingerprint: str,
) -> Iterator[int]:
    """Hold a pair-stable flock inherited by the complete owned PGID."""

    if pair_id not in runtime.PAIR_ORDER or not HEX64.fullmatch(run_fingerprint):
        raise CoordinatorError("worker lease requires a frozen pair and fingerprint")
    directory = runtime.ensure_directory(root, LEASE_ROOT_REL)
    path = directory / f"{pair_id}.lock"
    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        opened = os.fstat(descriptor)
        named = os.lstat(path)
        if (
            not stat.S_ISREG(opened.st_mode) or not stat.S_ISREG(named.st_mode)
            or opened.st_nlink != 1
            or (opened.st_dev, opened.st_ino) != (named.st_dev, named.st_ino)
        ):
            raise CoordinatorError("worker lease is not one private stable regular file")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno in {errno.EACCES, errno.EAGAIN}:
                raise CoordinatorError(f"pair {pair_id} has an active/orphaned worker") from error
            raise
        payload = {
            "schema_version": "sleep-atlas-track-b-finemapping-worker-lease.1",
            "pair_id": pair_id, "run_fingerprint": run_fingerprint,
            "authorization_nonce": os.urandom(32).hex(),
            "coordinator_sha256": runtime.stable_identity(root, SCRIPT)["sha256"],
        }
        content = _canonical_json(payload) + b"\n"
        os.ftruncate(descriptor, 0)
        os.lseek(descriptor, 0, os.SEEK_SET)
        offset = 0
        while offset < len(content):
            written = os.write(descriptor, content[offset:])
            if written <= 0:
                raise CoordinatorError("short worker-lease write")
            offset += written
        os.fsync(descriptor)
        runtime.fsync_directory(path.parent)
        # pass_fds makes the descriptor inheritable only for the requested
        # exec. Script 148 explicitly forwards it to its R LD child as well.
        yield descriptor
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


@contextlib.contextmanager
def coordinator_lease(runtime: Any, root: Path) -> Iterator[None]:
    directory = runtime.ensure_directory(root, LEASE_ROOT_REL)
    path = directory / "coordinator.lock"
    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        opened, named = os.fstat(descriptor), os.lstat(path)
        if (
            not stat.S_ISREG(opened.st_mode) or opened.st_nlink != 1
            or (opened.st_dev, opened.st_ino) != (named.st_dev, named.st_ino)
        ):
            raise CoordinatorError("coordinator lease is unsafe")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno in {errno.EACCES, errno.EAGAIN}:
                raise CoordinatorError("another Track B fine-mapping coordinator is active") from error
            raise
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def publish_directory_no_replace(runtime: Any, root: Path, stage: Path, destination: Path) -> None:
    """Kernel-atomic same-filesystem directory rename with no replacement."""

    runtime.safe_path(root, stage, "publication stage", must_exist=True)
    runtime.safe_path(root, destination.parent, "publication parent", must_exist=True)
    runtime.fsync_directory(stage)
    try:
        libc = ctypes.CDLL(ctypes.util.find_library("c") or None, use_errno=True)
    except OSError as error:
        raise CoordinatorError(f"cannot load libc for no-replace publication: {error}") from error
    source_raw, destination_raw = os.fsencode(stage), os.fsencode(destination)
    if sys.platform == "darwin":
        try:
            rename = libc.renamex_np
        except AttributeError as error:
            raise CoordinatorError("renamex_np unavailable for atomic no-replace publication") from error
        rename.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        rename.restype = ctypes.c_int
        result = rename(source_raw, destination_raw, 0x00000004)
    elif sys.platform.startswith("linux") and hasattr(libc, "renameat2"):
        rename = libc.renameat2
        rename.argtypes = [
            ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint,
        ]
        rename.restype = ctypes.c_int
        result = rename(-100, source_raw, -100, destination_raw, 0x00000001)
    else:
        raise CoordinatorError("kernel atomic no-replace directory publication is unavailable")
    if result:
        number = ctypes.get_errno()
        if number in {errno.EEXIST, errno.ENOTEMPTY}:
            raise CoordinatorError("no-replace directory destination already exists")
        raise CoordinatorError(f"atomic directory publication failed: {os.strerror(number)}")
    runtime.fsync_directory(destination.parent)


def _stream_publish_no_replace(
    runtime: Any, root: Path, source: Path, destination_rel: Path,
    expected: Mapping[str, Any],
) -> dict[str, object]:
    """Copy one package member exclusively without retaining it in RAM."""

    source_identity = runtime.stable_identity(root, source)
    if source_identity["bytes"] != expected["bytes"] or source_identity["sha256"] != expected["sha256"]:
        raise CoordinatorError("publication package member differs from its seal")
    destination = runtime.safe_path(root, destination_rel, "canonical publication")
    if destination.exists() or destination.is_symlink():
        observed = runtime.stable_identity(root, destination_rel)
        if observed["bytes"] != expected["bytes"] or observed["sha256"] != expected["sha256"]:
            raise CoordinatorError(f"preexisting canonical artifact conflicts: {destination_rel}")
        return observed
    parent = runtime.ensure_directory(root, destination_rel.parent)
    temporary = parent / f".{destination.name}.{os.getpid()}.{os.urandom(8).hex()}.tmp"
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    descriptor = os.open(temporary, flags, 0o600)
    digest = hashlib.sha256()
    count = 0
    try:
        with source.open("rb") as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                digest.update(block)
                count += len(block)
                offset = 0
                while offset < len(block):
                    offset += os.write(descriptor, block[offset:])
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    if count != expected["bytes"] or digest.hexdigest() != expected["sha256"]:
        temporary.unlink(missing_ok=True)
        raise CoordinatorError("publication source changed during streaming copy")
    try:
        os.link(temporary, destination, follow_symlinks=False)
        runtime.fsync_directory(parent)
    except FileExistsError:
        temporary.unlink(missing_ok=True)
        observed = runtime.stable_identity(root, destination_rel)
        if observed["bytes"] != expected["bytes"] or observed["sha256"] != expected["sha256"]:
            raise CoordinatorError("canonical destination raced with conflicting content")
        return observed
    temporary.unlink()
    runtime.fsync_directory(parent)
    return runtime.stable_identity(root, destination_rel)


def _new_stage(runtime: Any, root: Path, parent_rel: Path, prefix: str) -> Path:
    if not SAFE_ID.fullmatch(prefix):
        raise CoordinatorError("unsafe stage prefix")
    parent = runtime.ensure_directory(root, parent_rel)
    stage = parent / f".{prefix}.{os.getpid()}.{os.urandom(12).hex()}.building"
    os.mkdir(stage, 0o700)
    runtime.fsync_directory(parent)
    return stage


def _write_file(runtime: Any, path: Path, content: bytes) -> None:
    runtime._write_new_file(path, content)  # exact O_EXCL/fsync primitive from 148
    runtime.fsync_directory(path.parent)


def _publish_quarantine_metadata(
    runtime: Any, root: Path, path: Path, content: bytes,
) -> dict[str, object]:
    """Install a complete quarantine metadata file by no-replace hard link."""

    return runtime.publish_bytes_no_replace(root, path.relative_to(root), content)


def materialization_fingerprint(
    runtime: Any, family: Mapping[str, Any], pre_row: Mapping[str, str], root: Path = ROOT,
) -> str:
    scripts = {
        str(relative): runtime.stable_identity(root, relative)
        for relative in (RUNTIME_REL, Path("scripts/149_extract_track_b_signed_ld.R"), SCRIPT)
    }
    return _digest({
        "schema_version": MEASUREMENT_SCHEMA,
        "locus_entry_id": pre_row["locus_entry_id"],
        "pre_row": dict(pre_row), "pre_family_lock": family["lock_identity"],
        "policy": family["policy_identity"], "scripts": scripts,
        "ram_amendment": ram_measurement_amendment(runtime, root),
        "execution_amendment": runtime.EXECUTION_AMENDMENT,
        "execution_amendment_sha256": runtime.EXECUTION_AMENDMENT_SHA256,
        "physical_envelope_bytes": PHYSICAL_ENVELOPE_BYTES,
        "minimum_host_reserve_bytes": RESERVE_BYTES,
    })


def _quarantine_directory(
    runtime: Any, root: Path, source: Path, *, category: str, locus_entry_id: str,
    reason: str,
) -> dict[str, object]:
    """Atomically retain and seal every byte of an interrupted directory."""

    runtime.safe_path(root, source, "interrupted directory", must_exist=True)
    observed = os.lstat(source)
    if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
        raise CoordinatorError("interrupted state is not a real directory")
    category_id = re.sub(r"[^A-Za-z0-9_.-]", "_", category)
    destination_parent = runtime.ensure_directory(root, QUARANTINE_ROOT_REL / category_id)
    token = hashlib.sha256(
        f"{source}:{observed.st_dev}:{observed.st_ino}:{observed.st_ctime_ns}:{reason}".encode()
    ).hexdigest()[:20]
    destination = destination_parent / f"{locus_entry_id}.{token}"
    if destination.exists() or destination.is_symlink():
        raise CoordinatorError("quarantine destination collision")
    reserved = {
        "quarantine.evidence.json", "quarantine.receipt.json", "quarantine.seal.json",
    }
    if any((source / name).exists() or (source / name).is_symlink() for name in reserved):
        raise CoordinatorError("interrupted source collides with quarantine metadata names")
    evidence_records = _quarantine_evidence_records(runtime, root, source)
    publish_directory_no_replace(runtime, root, source, destination)
    if _quarantine_evidence_records(runtime, root, destination) != evidence_records:
        raise CoordinatorError("quarantined evidence changed across atomic move")
    evidence = {
        "schema_version": "sleep-atlas-track-b-finemapping-quarantine-evidence.1",
        "analysis_id": ANALYSIS_ID, "records": evidence_records,
    }
    evidence_content = json.dumps(evidence, indent=2, sort_keys=True).encode() + b"\n"
    evidence_identity = _publish_quarantine_metadata(
        runtime, root, destination / "quarantine.evidence.json", evidence_content,
    )
    receipt = {
        "schema_version": "sleep-atlas-track-b-finemapping-interruption-quarantine.1",
        "analysis_id": ANALYSIS_ID, "locus_entry_id": locus_entry_id,
        "category": category, "reason": reason, "moved_utc": _utc_now(),
        "source_path": str(source.relative_to(root)),
        "quarantine_path": str(destination.relative_to(root)),
        "scientific_output_deleted": False,
        "recovered_after_interruption": False,
        "evidence_manifest": evidence_identity,
    }
    receipt_identity = _publish_quarantine_metadata(
        runtime, root, destination / "quarantine.receipt.json",
        json.dumps(receipt, indent=2, sort_keys=True).encode() + b"\n",
    )
    seal = {
        "schema_version": "sleep-atlas-track-b-finemapping-interruption-quarantine.1",
        "analysis_id": ANALYSIS_ID, "locus_entry_id": locus_entry_id,
        "category": category, "evidence_manifest": evidence_identity,
        "receipt": receipt_identity,
    }
    return _publish_quarantine_metadata(
        runtime, root, destination / "quarantine.seal.json",
        json.dumps(seal, indent=2, sort_keys=True).encode() + b"\n",
    )


def _quarantine_evidence_records(
    runtime: Any, root: Path, directory: Path,
) -> dict[str, dict[str, object]]:
    """Hash every retained regular file, excluding quarantine metadata itself."""

    reserved = {
        "quarantine.evidence.json", "quarantine.receipt.json", "quarantine.seal.json",
    }
    records: dict[str, dict[str, object]] = {}
    for path in sorted(directory.rglob("*"), key=lambda item: str(item.relative_to(directory))):
        relative = path.relative_to(directory)
        observed = os.lstat(path)
        if len(relative.parts) == 1 and relative.name in reserved:
            continue
        if stat.S_ISLNK(observed.st_mode) or not (
            stat.S_ISREG(observed.st_mode) or stat.S_ISDIR(observed.st_mode)
        ):
            raise CoordinatorError("quarantine evidence contains an unsafe filesystem entry")
        if stat.S_ISDIR(observed.st_mode):
            continue
        identity = runtime.stable_identity(root, path, allow_empty=True)
        records[str(relative)] = {
            "bytes": identity["bytes"], "sha256": identity["sha256"],
        }
    return records


def verify_quarantine_directory(
    runtime: Any, root: Path, directory: Path, *, category: str, locus_entry_id: str,
) -> dict[str, object]:
    """Side-effect-free validation of one sealed retained interruption."""

    runtime.safe_path(root, directory, "quarantine directory", must_exist=True)
    observed = os.lstat(directory)
    if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
        raise CoordinatorError("quarantine entry is not a real directory")
    seal, seal_identity = runtime.read_json(root, directory / "quarantine.seal.json")
    receipt, receipt_identity = runtime.read_json(root, directory / "quarantine.receipt.json")
    evidence, evidence_identity = runtime.read_json(root, directory / "quarantine.evidence.json")
    if (
        set(seal) != {
            "schema_version", "analysis_id", "locus_entry_id", "category",
            "evidence_manifest", "receipt",
        }
        or seal.get("schema_version")
        != "sleep-atlas-track-b-finemapping-interruption-quarantine.1"
        or seal.get("analysis_id") != ANALYSIS_ID
        or seal.get("locus_entry_id") != locus_entry_id
        or seal.get("category") != category
        or not runtime.identity_matches(seal.get("receipt"), receipt_identity)
        or not runtime.identity_matches(seal.get("evidence_manifest"), evidence_identity)
        or set(receipt) != {
            "schema_version", "analysis_id", "locus_entry_id", "category", "reason",
            "moved_utc", "source_path", "quarantine_path", "scientific_output_deleted",
            "recovered_after_interruption", "evidence_manifest",
        }
        or receipt.get("schema_version")
        != "sleep-atlas-track-b-finemapping-interruption-quarantine.1"
        or receipt.get("analysis_id") != ANALYSIS_ID
        or receipt.get("locus_entry_id") != locus_entry_id
        or receipt.get("category") != category
        or not isinstance(receipt.get("reason"), str) or not receipt["reason"]
        or receipt.get("quarantine_path") != str(directory.relative_to(root))
        or receipt.get("scientific_output_deleted") is not False
        or not isinstance(receipt.get("recovered_after_interruption"), bool)
        or not runtime.identity_matches(receipt.get("evidence_manifest"), evidence_identity)
        or set(evidence) != {"schema_version", "analysis_id", "records"}
        or evidence.get("schema_version")
        != "sleep-atlas-track-b-finemapping-quarantine-evidence.1"
        or evidence.get("analysis_id") != ANALYSIS_ID
        or not isinstance(evidence.get("records"), dict)
        or evidence["records"] != _quarantine_evidence_records(runtime, root, directory)
    ):
        raise CoordinatorError("sealed quarantine evidence or receipt drifted")
    return {
        "seal_identity": seal_identity, "receipt_identity": receipt_identity,
        "evidence_identity": evidence_identity, "receipt": receipt,
    }


def _complete_recovered_quarantine_directory(
    runtime: Any, root: Path, directory: Path, *, category: str, locus_entry_id: str,
) -> dict[str, object]:
    """Finish sealing a directory moved just before a coordinator crash."""

    seal_path = directory / "quarantine.seal.json"
    if seal_path.exists() or seal_path.is_symlink():
        return verify_quarantine_directory(
            runtime, root, directory, category=category, locus_entry_id=locus_entry_id,
        )
    evidence_path = directory / "quarantine.evidence.json"
    receipt_path = directory / "quarantine.receipt.json"
    if evidence_path.is_symlink() or receipt_path.is_symlink():
        raise CoordinatorError("recovered quarantine metadata is unsafe")
    if evidence_path.exists():
        evidence, evidence_identity = runtime.read_json(root, evidence_path)
        if (
            evidence.get("schema_version")
            != "sleep-atlas-track-b-finemapping-quarantine-evidence.1"
            or evidence.get("analysis_id") != ANALYSIS_ID
            or evidence.get("records") != _quarantine_evidence_records(runtime, root, directory)
        ):
            raise CoordinatorError("interrupted quarantine evidence manifest drifted")
    else:
        evidence = {
            "schema_version": "sleep-atlas-track-b-finemapping-quarantine-evidence.1",
            "analysis_id": ANALYSIS_ID,
            "records": _quarantine_evidence_records(runtime, root, directory),
        }
        evidence_identity = _publish_quarantine_metadata(
            runtime, root, evidence_path,
            json.dumps(evidence, indent=2, sort_keys=True).encode() + b"\n",
        )
    if receipt_path.exists():
        receipt, receipt_identity = runtime.read_json(root, receipt_path)
        if (
            receipt.get("analysis_id") != ANALYSIS_ID
            or receipt.get("locus_entry_id") != locus_entry_id
            or receipt.get("category") != category
            or receipt.get("scientific_output_deleted") is not False
            or not runtime.identity_matches(receipt.get("evidence_manifest"), evidence_identity)
        ):
            raise CoordinatorError("interrupted quarantine receipt drifted")
    else:
        receipt = {
            "schema_version": "sleep-atlas-track-b-finemapping-interruption-quarantine.1",
            "analysis_id": ANALYSIS_ID, "locus_entry_id": locus_entry_id,
            "category": category,
            "reason": "RECOVERED_AFTER_ATOMIC_MOVE_INTERRUPTED_BEFORE_QUARANTINE_SEAL",
            "moved_utc": _utc_now(),
            "source_path": "UNKNOWN_AFTER_ATOMIC_QUARANTINE_MOVE",
            "quarantine_path": str(directory.relative_to(root)),
            "scientific_output_deleted": False,
            "recovered_after_interruption": True,
            "evidence_manifest": evidence_identity,
        }
        receipt_identity = _publish_quarantine_metadata(
            runtime, root, receipt_path,
            json.dumps(receipt, indent=2, sort_keys=True).encode() + b"\n",
        )
    seal = {
        "schema_version": "sleep-atlas-track-b-finemapping-interruption-quarantine.1",
        "analysis_id": ANALYSIS_ID, "locus_entry_id": locus_entry_id,
        "category": category, "evidence_manifest": evidence_identity,
        "receipt": receipt_identity,
    }
    _publish_quarantine_metadata(
        runtime, root, seal_path, json.dumps(seal, indent=2, sort_keys=True).encode() + b"\n",
    )
    return verify_quarantine_directory(
        runtime, root, directory, category=category, locus_entry_id=locus_entry_id,
    )


def _quarantine_materialization_input_stages(
    runtime: Any, root: Path, locus_entry_id: str, *, reason: str,
) -> list[dict[str, object]]:
    """Under the pair lease, retain/reconcile every hidden input-building attempt."""

    category = "materialization_input_attempt"
    input_root = runtime.safe_path(root, runtime.INPUT_ROOT_REL, "materialization input root")
    if input_root.exists():
        if not stat.S_ISDIR(os.lstat(input_root).st_mode) or input_root.is_symlink():
            raise CoordinatorError("materialization input root is unsafe")
        for child in sorted(input_root.glob(f".{locus_entry_id}.*.building")):
            _quarantine_directory(
                runtime, root, child, category=category,
                locus_entry_id=locus_entry_id, reason=reason,
            )
    category_root = runtime.safe_path(
        root, QUARANTINE_ROOT_REL / category, "materialization quarantine root",
    )
    if not category_root.exists():
        return []
    if not stat.S_ISDIR(os.lstat(category_root).st_mode) or category_root.is_symlink():
        raise CoordinatorError("materialization quarantine root is unsafe")
    output: list[dict[str, object]] = []
    for directory in sorted(category_root.glob(f"{locus_entry_id}.*")):
        output.append(_complete_recovered_quarantine_directory(
            runtime, root, directory, category=category, locus_entry_id=locus_entry_id,
        ))
    return output


def _verify_materialization_quarantine_seals(
    runtime: Any, root: Path, values: object, *, locus_entry_id: str,
) -> list[dict[str, object]]:
    """Revalidate receipt-bound materialization quarantine seals without writes."""

    if not isinstance(values, list):
        raise CoordinatorError("materialization quarantine seal family is malformed")
    output: list[dict[str, object]] = []
    observed_paths: set[str] = set()
    expected_parent = QUARANTINE_ROOT_REL / "materialization_input_attempt"
    for expected in values:
        if not isinstance(expected, Mapping) or set(expected) != {"path", "bytes", "sha256"}:
            raise CoordinatorError("materialization quarantine identity is malformed")
        path = Path(str(expected["path"]))
        if (
            path.name != "quarantine.seal.json" or path.parent.parent != expected_parent
            or not path.parent.name.startswith(f"{locus_entry_id}.")
            or str(path) in observed_paths
        ):
            raise CoordinatorError("materialization quarantine identity path drifted")
        observed_paths.add(str(path))
        verified = verify_quarantine_directory(
            runtime, root, runtime.safe_path(root, path.parent, "materialization quarantine", must_exist=True),
            category="materialization_input_attempt", locus_entry_id=locus_entry_id,
        )
        if not runtime.identity_matches(verified["seal_identity"], expected):
            raise CoordinatorError("materialization quarantine seal identity drifted")
        output.append(verified)
    return output


def _publish_measurement_attempt(
    runtime: Any, root: Path, receipt: Mapping[str, Any], log: bytes, *, success: bool,
) -> dict[str, object]:
    locus = str(receipt["locus_entry_id"])
    fingerprint = str(receipt["materialization_fingerprint"])
    if success:
        destination_rel = MATERIALIZATION_METRICS_REL / locus
    else:
        attempt = str(receipt["attempt_id"])
        destination_rel = MATERIALIZATION_ATTEMPTS_REL / f"{locus}.{attempt}"
    if (root / destination_rel).exists() or (root / destination_rel).is_symlink():
        raise CoordinatorError(f"materialization receipt destination already exists: {destination_rel}")
    stage = _new_stage(runtime, root, STAGE_ROOT_REL / "measurement_receipts", locus)
    try:
        _write_file(runtime, stage / "process.log", log)
        log_identity = runtime.stable_identity(root, stage / "process.log", allow_empty=True)
        finalized = dict(receipt)
        finalized["process_log"] = {
            "path": "process.log", "bytes": log_identity["bytes"],
            "sha256": log_identity["sha256"],
        }
        receipt_content = json.dumps(finalized, indent=2, sort_keys=True).encode() + b"\n"
        _write_file(runtime, stage / "receipt.json", receipt_content)
        seal = {
            "schema_version": MEASUREMENT_SCHEMA,
            "materialization_fingerprint": fingerprint,
            "receipt_sha256": hashlib.sha256(receipt_content).hexdigest(),
            "process_log_sha256": log_identity["sha256"],
        }
        _write_file(runtime, stage / "measurement.seal.json", json.dumps(
            seal, indent=2, sort_keys=True,
        ).encode() + b"\n")
        destination_parent = runtime.ensure_directory(root, destination_rel.parent)
        publish_directory_no_replace(runtime, root, stage, destination_parent / destination_rel.name)
        return runtime.stable_identity(root, destination_rel / "measurement.seal.json")
    except BaseException:
        # Owned incomplete staging is retained for explicit reconciliation.
        raise


def _quarantine_measurement_stages(
    runtime: Any, root: Path, locus_entry_id: str,
) -> list[dict[str, object]]:
    parent = runtime.safe_path(
        root, STAGE_ROOT_REL / "measurement_receipts", "measurement receipt staging",
    )
    if not parent.exists():
        return []
    if not stat.S_ISDIR(os.lstat(parent).st_mode) or parent.is_symlink():
        raise CoordinatorError("measurement receipt staging root is unsafe")
    output: list[dict[str, object]] = []
    for child in sorted(parent.glob(f".{locus_entry_id}.*.building")):
        output.append(_quarantine_directory(
            runtime, root, child, category="interrupted_materialization_receipt",
            locus_entry_id=locus_entry_id,
            reason="UNPUBLISHED_MEASUREMENT_STAGE_AFTER_PAIR_LEASE_BECAME_AVAILABLE",
        ))
    return output


def verify_materialization_measurement(
    root: Path, locus_entry_id: str, *, family: Mapping[str, Any] | None = None,
    runtime: Any | None = None,
) -> dict[str, Any]:
    """Deeply revalidate one successful materialization and PGID benchmark."""

    runtime = runtime or _load_runtime(root)
    family = family or runtime.validate_pre_materialization_family(root)
    matches = [row for row in family["rows"] if row["locus_entry_id"] == locus_entry_id]
    if len(matches) != 1:
        raise CoordinatorError("measurement locus is outside the frozen pre-family")
    directory_rel = MATERIALIZATION_METRICS_REL / locus_entry_id
    directory = runtime.safe_path(root, directory_rel, "materialization measurement", must_exist=True)
    if not stat.S_ISDIR(os.lstat(directory).st_mode) or directory.is_symlink():
        raise CoordinatorError("materialization measurement is not a real directory")
    observed_names = {item.name for item in directory.iterdir()}
    if observed_names != {"process.log", "receipt.json", "measurement.seal.json"}:
        raise CoordinatorError("materialization measurement file family drifted")
    receipt, receipt_identity = runtime.read_json(root, directory_rel / "receipt.json")
    seal, seal_identity = runtime.read_json(root, directory_rel / "measurement.seal.json")
    log_identity = runtime.stable_identity(root, directory_rel / "process.log", allow_empty=True)
    expected_keys = {
        "schema_version", "analysis_id", "locus_entry_id", "pair_id", "attempt_id",
        "materialization_fingerprint", "terminal_status", "pre_family_lock", "policy",
        "ram_measurement_amendment", "host_physical_memory_bytes",
        "available_memory_before_bytes", "configured_envelope_bytes",
        "effective_envelope_bytes", "larger_host_continuation",
        "maximum_owned_pgid_rss_bytes",
        "minimum_host_reserve_bytes", "command", "started_utc", "finished_utc",
        "elapsed_monotonic_seconds", "peak_rss_bytes", "rss_measurement_backend",
        "positive_pgid_samples", "exit_code", "guard_failure", "bundle_lock",
        "prelaunch_quarantine_seals", "attempt_quarantine_seals", "process_log",
    }
    expected_fingerprint = materialization_fingerprint(runtime, family, matches[0], root)
    configured_raw = receipt.get("configured_envelope_bytes")
    configured = (
        int(configured_raw)
        if isinstance(configured_raw, int) and not isinstance(configured_raw, bool)
        else -1
    )
    continuation = receipt.get("larger_host_continuation") is True
    expected_command = [
        sys.executable, str(root / RUNTIME_REL), "--root", str(root),
        "--materialize-one", locus_entry_id, "--execute",
        "--memory-envelope-bytes", str(configured),
    ]
    if continuation:
        expected_command.append("--larger-host-continuation")
    if (
        set(receipt) != expected_keys or receipt.get("schema_version") != MEASUREMENT_SCHEMA
        or receipt.get("analysis_id") != ANALYSIS_ID
        or receipt.get("locus_entry_id") != locus_entry_id
        or receipt.get("pair_id") != matches[0]["pair_id"]
        or receipt.get("materialization_fingerprint") != expected_fingerprint
        or receipt.get("terminal_status") != "COMPLETE_MATERIALIZED"
        or receipt.get("pre_family_lock") != family["lock_identity"]
        or receipt.get("policy") != family["policy_identity"]
        or receipt.get("ram_measurement_amendment") != ram_measurement_amendment(runtime, root)
        or receipt.get("command") != expected_command
        or (configured > PHYSICAL_ENVELOPE_BYTES and not continuation)
        or not isinstance(receipt.get("configured_envelope_bytes"), int)
        or int(receipt["configured_envelope_bytes"]) < PHYSICAL_ENVELOPE_BYTES
        or not isinstance(receipt.get("larger_host_continuation"), bool)
        or receipt.get("effective_envelope_bytes") != min(
            int(receipt.get("host_physical_memory_bytes", 0)),
            int(receipt["configured_envelope_bytes"]),
        )
        or receipt.get("maximum_owned_pgid_rss_bytes")
        != int(receipt["effective_envelope_bytes"]) - RESERVE_BYTES
        or receipt.get("rss_measurement_backend") != PGID_BACKEND
        or not isinstance(receipt.get("positive_pgid_samples"), int)
        or int(receipt["positive_pgid_samples"]) <= 0
        or not isinstance(receipt.get("peak_rss_bytes"), int)
        or int(receipt["peak_rss_bytes"]) <= 0
        or receipt.get("exit_code") != 0 or receipt.get("guard_failure") is not None
        or receipt.get("process_log") != {
            "path": "process.log", "bytes": log_identity["bytes"],
            "sha256": log_identity["sha256"],
        }
    ):
        raise CoordinatorError("successful materialization measurement drifted")
    if seal != {
        "schema_version": MEASUREMENT_SCHEMA,
        "materialization_fingerprint": expected_fingerprint,
        "receipt_sha256": receipt_identity["sha256"],
        "process_log_sha256": log_identity["sha256"],
    }:
        raise CoordinatorError("materialization measurement seal drifted")
    prelaunch_quarantines = _verify_materialization_quarantine_seals(
        runtime, root, receipt["prelaunch_quarantine_seals"], locus_entry_id=locus_entry_id,
    )
    attempt_quarantines = _verify_materialization_quarantine_seals(
        runtime, root, receipt["attempt_quarantine_seals"], locus_entry_id=locus_entry_id,
    )
    if {
        str(item["seal_identity"]["path"]) for item in prelaunch_quarantines
    } & {str(item["seal_identity"]["path"]) for item in attempt_quarantines}:
        raise CoordinatorError("materialization quarantine lineage is duplicated")
    if continuation:
        failures = _failed_materialization_attempts(runtime, root, family, locus_entry_id)
        if (
            not failures
            or any(
                attempt["terminal_status"] not in {"BLOCKED_BY_COMPUTE", "FAILED_RESOURCE_OOM"}
                for attempt in failures
            )
            or int(receipt["host_physical_memory_bytes"]) <= max(
                int(attempt["receipt"]["host_physical_memory_bytes"])
                for attempt in failures
            )
            or int(receipt["effective_envelope_bytes"]) <= max(
                int(attempt["receipt"]["effective_envelope_bytes"])
                for attempt in failures
            )
        ):
            raise CoordinatorError("successful larger-host continuation lacks a valid larger retry")
    bundle = runtime.verify_materialized_bundle(
        root, runtime.INPUT_ROOT_REL / locus_entry_id, family=family,
    )
    if receipt.get("bundle_lock") != bundle["bundle_lock_identity"]:
        raise CoordinatorError("measurement does not bind the current materialized bundle")
    return {
        "receipt": receipt, "receipt_identity": receipt_identity,
        "seal_identity": seal_identity, "bundle": bundle, "row": matches[0],
    }


def verify_failed_materialization_attempt(
    root: Path, attempt_rel: Path | str, *, family: Mapping[str, Any] | None = None,
    runtime: Any | None = None,
) -> dict[str, Any]:
    """Side-effect-free validation of one retained failed materialization attempt."""

    runtime = runtime or _load_runtime(root)
    family = family or runtime.validate_pre_materialization_family(root)
    attempt_rel = Path(attempt_rel)
    directory = runtime.safe_path(root, attempt_rel, "failed materialization attempt", must_exist=True)
    if not stat.S_ISDIR(os.lstat(directory).st_mode) or directory.is_symlink():
        raise CoordinatorError("failed materialization attempt is not a real directory")
    if {item.name for item in directory.iterdir()} != {
        "process.log", "receipt.json", "measurement.seal.json",
    }:
        raise CoordinatorError("failed materialization attempt file family drifted")
    receipt, receipt_identity = runtime.read_json(root, attempt_rel / "receipt.json")
    seal, seal_identity = runtime.read_json(root, attempt_rel / "measurement.seal.json")
    log_identity = runtime.stable_identity(root, attempt_rel / "process.log", allow_empty=True)
    locus = receipt.get("locus_entry_id")
    matches = [row for row in family["rows"] if row["locus_entry_id"] == locus]
    if len(matches) != 1:
        raise CoordinatorError("failed materialization locus is outside frozen family")
    expected_keys = {
        "schema_version", "analysis_id", "locus_entry_id", "pair_id", "attempt_id",
        "materialization_fingerprint", "terminal_status", "pre_family_lock", "policy",
        "ram_measurement_amendment", "host_physical_memory_bytes",
        "available_memory_before_bytes", "configured_envelope_bytes",
        "effective_envelope_bytes", "larger_host_continuation",
        "maximum_owned_pgid_rss_bytes", "minimum_host_reserve_bytes", "command",
        "started_utc", "finished_utc", "elapsed_monotonic_seconds", "peak_rss_bytes",
        "rss_measurement_backend", "positive_pgid_samples", "exit_code", "guard_failure",
        "bundle_lock", "prelaunch_quarantine_seals", "attempt_quarantine_seals",
        "process_log",
    }
    configured_raw = receipt.get("configured_envelope_bytes")
    configured = (
        int(configured_raw)
        if isinstance(configured_raw, int) and not isinstance(configured_raw, bool)
        else -1
    )
    continuation = receipt.get("larger_host_continuation") is True
    expected_command = [
        sys.executable, str(root / RUNTIME_REL), "--root", str(root),
        "--materialize-one", str(locus), "--execute",
        "--memory-envelope-bytes", str(configured),
    ]
    if continuation:
        expected_command.append("--larger-host-continuation")
    fingerprint = materialization_fingerprint(runtime, family, matches[0], root)
    terminal = receipt.get("terminal_status")
    if (
        set(receipt) != expected_keys or receipt.get("schema_version") != MEASUREMENT_SCHEMA
        or receipt.get("analysis_id") != ANALYSIS_ID
        or receipt.get("pair_id") != matches[0]["pair_id"]
        or receipt.get("materialization_fingerprint") != fingerprint
        or terminal not in {"BLOCKED_BY_COMPUTE", "FAILED_RESOURCE_OOM", "FAILED_VALIDATION"}
        or receipt.get("pre_family_lock") != family["lock_identity"]
        or receipt.get("policy") != family["policy_identity"]
        or receipt.get("ram_measurement_amendment") != ram_measurement_amendment(runtime, root)
        or configured < PHYSICAL_ENVELOPE_BYTES
        or (configured > PHYSICAL_ENVELOPE_BYTES and not continuation)
        or not isinstance(receipt.get("host_physical_memory_bytes"), int)
        or not isinstance(receipt.get("available_memory_before_bytes"), int)
        or receipt.get("effective_envelope_bytes") != min(
            int(receipt.get("host_physical_memory_bytes", 0)), configured,
        )
        or receipt.get("maximum_owned_pgid_rss_bytes")
        != int(receipt.get("effective_envelope_bytes", 0)) - RESERVE_BYTES
        or receipt.get("minimum_host_reserve_bytes") != RESERVE_BYTES
        or receipt.get("command") != expected_command
        or receipt.get("bundle_lock") is not None
        or receipt.get("process_log") != {
            "path": "process.log", "bytes": log_identity["bytes"],
            "sha256": log_identity["sha256"],
        }
        or attempt_rel.name != f"{locus}.{receipt.get('attempt_id')}"
    ):
        raise CoordinatorError("failed materialization receipt drifted")
    resource_blocked = receipt.get("exit_code") == 75 and receipt.get("guard_failure") is None
    resource_oom = (
        str(receipt.get("guard_failure") or "").startswith(
            "PHYSICAL_RAM_CAP_OR_HOST_RESERVE_EXCEEDED"
        )
        or receipt.get("exit_code") in {-9, 137}
    )
    if (
        (terminal == "BLOCKED_BY_COMPUTE" and not resource_blocked)
        or (terminal == "FAILED_RESOURCE_OOM" and not resource_oom)
        or (terminal == "FAILED_VALIDATION" and (resource_blocked or resource_oom))
    ):
        raise CoordinatorError("failed materialization terminal status contradicts evidence")
    if seal != {
        "schema_version": MEASUREMENT_SCHEMA,
        "materialization_fingerprint": fingerprint,
        "receipt_sha256": receipt_identity["sha256"],
        "process_log_sha256": log_identity["sha256"],
    }:
        raise CoordinatorError("failed materialization measurement seal drifted")
    prelaunch_quarantines = _verify_materialization_quarantine_seals(
        runtime, root, receipt["prelaunch_quarantine_seals"], locus_entry_id=str(locus),
    )
    attempt_quarantines = _verify_materialization_quarantine_seals(
        runtime, root, receipt["attempt_quarantine_seals"], locus_entry_id=str(locus),
    )
    if {
        str(item["seal_identity"]["path"]) for item in prelaunch_quarantines
    } & {str(item["seal_identity"]["path"]) for item in attempt_quarantines}:
        raise CoordinatorError("failed materialization quarantine lineage is duplicated")
    return {
        "terminal_status": terminal, "receipt": receipt,
        "receipt_identity": receipt_identity, "seal_identity": seal_identity,
        "row": matches[0],
    }


def _failed_materialization_attempts(
    runtime: Any, root: Path, family: Mapping[str, Any], locus: str,
) -> list[dict[str, Any]]:
    parent = runtime.safe_path(root, MATERIALIZATION_ATTEMPTS_REL, "materialization attempts")
    if not parent.exists():
        return []
    if not stat.S_ISDIR(os.lstat(parent).st_mode) or parent.is_symlink():
        raise CoordinatorError("materialization-attempt root is unsafe")
    output: list[dict[str, Any]] = []
    for child in sorted(parent.iterdir(), key=lambda item: item.name):
        if child.name.startswith(f"{locus}."):
            output.append(verify_failed_materialization_attempt(
                root, child.relative_to(root), family=family, runtime=runtime,
            ))
    return output


def _attempt_id() -> str:
    return dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + "." + os.urandom(6).hex()


def validate_larger_host_materialization_continuation(
    prior_failures: Sequence[Mapping[str, Any]], *, current_physical_bytes: int,
    configured_envelope_bytes: int,
) -> int:
    """Authorize only a strictly larger unchanged-locus resource continuation."""

    if not prior_failures or any(
        attempt.get("terminal_status") not in {"BLOCKED_BY_COMPUTE", "FAILED_RESOURCE_OOM"}
        for attempt in prior_failures
    ):
        raise CoordinatorError("materialization continuation lacks only retained resource failures")
    previous_physical = max(
        int(attempt["receipt"]["host_physical_memory_bytes"])
        for attempt in prior_failures
    )
    previous_effective = max(
        int(attempt["receipt"]["effective_envelope_bytes"])
        for attempt in prior_failures
    )
    current_effective = min(int(current_physical_bytes), int(configured_envelope_bytes))
    if (
        int(current_physical_bytes) <= previous_physical
        or current_effective <= previous_effective
        or int(configured_envelope_bytes) <= PHYSICAL_ENVELOPE_BYTES
    ):
        raise CoordinatorError(
            "materialization retry requires a strictly larger host and effective envelope"
        )
    return current_effective


def benchmark_materialization(
    runtime: Any, root: Path, family: Mapping[str, Any], pre_row: Mapping[str, str],
    *, configured_envelope_bytes: int = PHYSICAL_ENVELOPE_BYTES,
    larger_host_continuation: bool = False,
) -> dict[str, Any]:
    """Materialize one unchanged whole locus in a measured fresh process."""

    locus = pre_row["locus_entry_id"]
    configured_envelope_bytes = int(configured_envelope_bytes)
    if configured_envelope_bytes < PHYSICAL_ENVELOPE_BYTES:
        raise CoordinatorError("configured materialization envelope is below 8 GiB")
    if configured_envelope_bytes > PHYSICAL_ENVELOPE_BYTES and not larger_host_continuation:
        raise CoordinatorError("larger materialization envelope requires explicit continuation")
    prior_failures = _failed_materialization_attempts(runtime, root, family, locus)
    success_path = root / MATERIALIZATION_METRICS_REL / locus
    if success_path.exists() or success_path.is_symlink():
        return {"terminal_status": "COMPLETE_MATERIALIZED", **verify_materialization_measurement(
            root, locus, family=family, runtime=runtime,
        )}
    if larger_host_continuation and not prior_failures:
        raise CoordinatorError(
            "larger-host materialization continuation requires a retained resource failure"
        )
    if prior_failures:
        nonresource = [
            attempt for attempt in prior_failures
            if attempt["terminal_status"] not in {"BLOCKED_BY_COMPUTE", "FAILED_RESOURCE_OOM"}
        ]
        if nonresource:
            return nonresource[-1]
        if not larger_host_continuation:
            return prior_failures[-1]
        current_physical = runtime.host_physical_memory_bytes()
        validate_larger_host_materialization_continuation(
            prior_failures, current_physical_bytes=current_physical,
            configured_envelope_bytes=configured_envelope_bytes,
        )
    fingerprint = materialization_fingerprint(runtime, family, pre_row, root)
    bundle_path = root / runtime.INPUT_ROOT_REL / locus
    with worker_lease(runtime, root, pre_row["pair_id"], fingerprint) as lease_fd:
        _quarantine_measurement_stages(runtime, root, locus)
        prelaunch_quarantines = _quarantine_materialization_input_stages(
            runtime, root, locus,
            reason="ORPHAN_INPUT_BUILDING_STAGE_FOUND_AFTER_PAIR_LEASE_BECAME_AVAILABLE",
        )
        prelaunch_quarantine_seals = sorted(
            (item["seal_identity"] for item in prelaunch_quarantines),
            key=lambda item: str(item["path"]),
        )
        # An input bundle without a successful benchmark cannot be used. Before
        # the materialized-family lock exists it is safely retained in
        # quarantine and rematerialized unchanged so both phases are measured.
        if bundle_path.exists() or bundle_path.is_symlink():
            if (root / runtime.MATERIALIZED_LOCK_REL).exists():
                raise CoordinatorError("locked materialized family lacks its required benchmark")
            _quarantine_directory(
                runtime, root, bundle_path, category="unmeasured_materialization",
                locus_entry_id=locus, reason="SEALED_OR_PARTIAL_BUNDLE_WITHOUT_SUCCESS_MEASUREMENT",
            )
        physical = runtime.host_physical_memory_bytes()
        available = runtime._memory_available_bytes()
        effective_envelope = min(physical, configured_envelope_bytes)
        live_limit = effective_envelope - RESERVE_BYTES
        if live_limit <= 0:
            raise CoordinatorError("materialization host cannot retain the locked reserve")
        attempt = _attempt_id()
        command = [
            sys.executable, str(root / RUNTIME_REL), "--root", str(root),
            "--materialize-one", locus, "--execute",
            "--memory-envelope-bytes", str(configured_envelope_bytes),
        ]
        if larger_host_continuation:
            command.append("--larger-host-continuation")
        measured = monitored_command(
            command, root=root, inherited_fds=(lease_fd,),
            maximum_owned_pgid_rss_bytes=live_limit,
            extra_environment={
                "TRACK_B_FINEMAP_WORKER_LEASE_FD": str(lease_fd),
                "TRACK_B_FINEMAP_MEASURED_PHASE": "MATERIALIZATION",
            },
            memory_available=runtime._memory_available_bytes,
        )
        terminal = "FAILED_VALIDATION"
        bundle_lock: dict[str, object] | None = None
        if measured["exit_code"] == 75 and measured["guard_failure"] is None:
            terminal = "BLOCKED_BY_COMPUTE"
        elif (
            str(measured["guard_failure"] or "").startswith(
                "PHYSICAL_RAM_CAP_OR_HOST_RESERVE_EXCEEDED"
            )
            or measured["exit_code"] in {-9, 137}
        ):
            terminal = "FAILED_RESOURCE_OOM"
        elif measured["exit_code"] == 0 and measured["success_measurement"]:
            verified = runtime.verify_materialized_bundle(
                root, runtime.INPUT_ROOT_REL / locus, family=family,
            )
            bundle_lock = verified["bundle_lock_identity"]
            terminal = "COMPLETE_MATERIALIZED"
        post_attempt_quarantines = _quarantine_materialization_input_stages(
            runtime, root, locus,
            reason=f"RETAINED_INPUT_BUILDING_STAGE_AFTER_{terminal}",
        )
        prelaunch_paths = {str(item["path"]) for item in prelaunch_quarantine_seals}
        attempt_quarantine_seals = sorted(
            (
                item["seal_identity"] for item in post_attempt_quarantines
                if str(item["seal_identity"]["path"]) not in prelaunch_paths
            ),
            key=lambda item: str(item["path"]),
        )
        receipt = {
            "schema_version": MEASUREMENT_SCHEMA, "analysis_id": ANALYSIS_ID,
            "locus_entry_id": locus, "pair_id": pre_row["pair_id"], "attempt_id": attempt,
            "materialization_fingerprint": fingerprint, "terminal_status": terminal,
            "pre_family_lock": family["lock_identity"], "policy": family["policy_identity"],
            "ram_measurement_amendment": ram_measurement_amendment(runtime, root),
            "host_physical_memory_bytes": physical,
            "available_memory_before_bytes": available,
            "configured_envelope_bytes": configured_envelope_bytes,
            "effective_envelope_bytes": effective_envelope,
            "larger_host_continuation": larger_host_continuation,
            "maximum_owned_pgid_rss_bytes": live_limit,
            "minimum_host_reserve_bytes": RESERVE_BYTES,
            "command": command, "started_utc": measured["started_utc"],
            "finished_utc": measured["finished_utc"],
            "elapsed_monotonic_seconds": measured["elapsed_monotonic_seconds"],
            "peak_rss_bytes": measured["peak_rss_bytes"],
            "rss_measurement_backend": measured["rss_measurement_backend"],
            "positive_pgid_samples": measured["positive_pgid_samples"],
            "exit_code": measured["exit_code"], "guard_failure": measured["guard_failure"],
            "bundle_lock": bundle_lock,
            "prelaunch_quarantine_seals": prelaunch_quarantine_seals,
            "attempt_quarantine_seals": attempt_quarantine_seals,
        }
        seal_identity = _publish_measurement_attempt(
            runtime, root, receipt, measured["log"], success=terminal == "COMPLETE_MATERIALIZED",
        )
        if terminal != "COMPLETE_MATERIALIZED":
            return verify_failed_materialization_attempt(
                root, MATERIALIZATION_ATTEMPTS_REL / f"{locus}.{attempt}",
                family=family, runtime=runtime,
            )
    return {"terminal_status": "COMPLETE_MATERIALIZED", **verify_materialization_measurement(
        root, locus, family=family, runtime=runtime,
    )}


def run_fingerprint(
    runtime: Any, materialized: Mapping[str, Any], manifest_row: Mapping[str, str],
    task_validation: Mapping[str, Any], root: Path = ROOT,
) -> str:
    """Fingerprint the unchanged scientific task, never host capacity/results."""

    scripts = runtime._runtime_script_identities(root)
    task = task_validation["task"]
    return _digest({
        "schema_version": RUN_SCHEMA, "analysis_id": ANALYSIS_ID,
        "manifest_lock": materialized["lock_identity"],
        "manifest_row": dict(manifest_row),
        "task_identity": task_validation["task_identity"],
        "input_identities": task_validation["identities"],
        "policy": materialized["pre_family"]["policy_identity"],
        "engine": runtime.stable_identity(root, R_ENGINE_REL),
        "scripts": scripts, "scalar_n_convention": runtime.SCALAR_N_CONVENTION,
        "execution_amendment": runtime.EXECUTION_AMENDMENT,
        "execution_amendment_sha256": runtime.EXECUTION_AMENDMENT_SHA256,
        "runtime_package_lock": task_validation["identities"]["runtime_package_lock"],
        "trait1_scalar_N": task["trait1_scalar_N"],
        "trait1_scalar_N_rule": task["trait1_scalar_N_rule"],
        "trait1_per_snp_N_sha256": task["trait1_per_snp_N_sha256"],
        "trait2_scalar_N": task["trait2_scalar_N"],
        "trait2_scalar_N_rule": task["trait2_scalar_N_rule"],
        "trait2_per_snp_N_sha256": task["trait2_per_snp_N_sha256"],
        "trait1_N_dispersion": {
            key: value for key, value in task.items()
            if key.startswith("trait1_per_snp_N_") or key == "trait1_N_dispersion_status"
        },
        "trait2_N_dispersion": {
            key: value for key, value in task.items()
            if key.startswith("trait2_per_snp_N_") or key == "trait2_N_dispersion_status"
        },
        "model": {
            "method": "susieR::susie_rss_via_coloc::runsusie",
            "susieR": "0.14.2", "L": 10, "coverage": 0.95,
            "min_abs_corr": 0.5, "maxit": 1000, "prior": "FLAT",
            "estimate_residual_variance": False,
            "coloc": "5.2.3::coloc.susie", "p1": 1e-4, "p2": 1e-4,
            "p12_grid": [1e-6, 5e-6, 1e-5, 5e-5], "abf_fallback": False,
        },
    })


def _clean_error(value: object) -> str:
    text = str(value or "NA").replace("\t", " ").replace("\r", " ").replace("\n", " ")
    return text[:10000] or "NA"


def _output_manifest(
    runtime: Any, root: Path, stage: Path, *, complete_engine_bundle: bool,
    terminal_status: str,
) -> tuple[dict[str, Any], bytes]:
    records: dict[str, dict[str, object]] = {}
    for path in sorted(stage.iterdir(), key=lambda item: item.name):
        observed = os.lstat(path)
        if not stat.S_ISREG(observed.st_mode) or stat.S_ISLNK(observed.st_mode) or observed.st_nlink != 1:
            raise CoordinatorError(f"run staging contains an unsafe entry: {path.name}")
        runtime.fsync_file(path)
        identity = runtime.stable_identity(root, path, allow_empty=(path.name == "process.log"))
        records[path.name] = {"bytes": identity["bytes"], "sha256": identity["sha256"]}
    payload = {
        "schema_version": "sleep-atlas-track-b-finemapping-run-output-manifest.1",
        "analysis_id": ANALYSIS_ID, "terminal_status": terminal_status,
        "complete_engine_bundle": complete_engine_bundle, "artifacts": records,
    }
    content = json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n"
    return payload, content


def _resource_row(
    runtime: Any, root: Path, task_validation: Mapping[str, Any], fingerprint: str, attempt: str,
    materialization: Mapping[str, Any], engine_measurement: Mapping[str, Any],
    *, physical: int, available: int, estimated: int, admission_status: str,
    terminal_status: str, output_manifest_sha256: str, error: str,
) -> dict[str, str]:
    task = task_validation["task"]
    material_receipt = materialization["receipt"]
    material_peak = int(material_receipt["peak_rss_bytes"])
    engine_peak = int(engine_measurement.get("peak_rss_bytes", 0))
    material_runtime = float(material_receipt["elapsed_monotonic_seconds"])
    engine_runtime = float(engine_measurement.get("elapsed_monotonic_seconds", 0.0))
    successful_backend = engine_measurement.get("rss_measurement_backend")
    if terminal_status == "COMPLETE_SEALED" and successful_backend != PGID_BACKEND:
        raise CoordinatorError("successful R run lacks authoritative PGID aggregate RSS")
    return {
        "analysis_id": ANALYSIS_ID, "run_fingerprint": fingerprint,
        "locus_entry_id": task["locus_entry_id"], "pair_id": task["pair_id"],
        "attempt_id": attempt,
        "worker_environment_id": f"{os.uname().nodename}:{sys.platform}",
        "process_isolation": "TWO_FRESH_OS_SESSIONS_MATERIALIZATION_THEN_R;NO_CONCURRENT_LOCI",
        "variant_count": task["variant_count"], "ld_dimension": task["variant_count"],
        "estimated_peak_rss_bytes": str(estimated),
        "host_physical_memory_bytes": str(physical),
        "available_memory_before_bytes": str(available),
        "admission_status": admission_status,
        "started_utc": str(material_receipt["started_utc"]),
        "finished_utc": str(engine_measurement.get("finished_utc", _utc_now())),
        "elapsed_monotonic_seconds": format(material_runtime + engine_runtime, ".17g"),
        "peak_rss_bytes": str(max(material_peak, engine_peak)),
        "rss_measurement_backend": PGID_BACKEND if terminal_status == "COMPLETE_SEALED" else str(
            successful_backend or material_receipt["rss_measurement_backend"]
        ),
        "exit_code": str(engine_measurement.get("exit_code", "NA")),
        "terminal_status": terminal_status,
        "task_sha256": task_validation["task_identity"]["sha256"],
        "trait1_summary_sha256": task["summary1_sha256"],
        "trait2_summary_sha256": task["summary2_sha256"],
        "variant_order_sha256": task["variant_order_sha256"],
        "signed_ld_sha256": task["signed_ld_sha256"],
        "engine_sha256": runtime.stable_identity(root, R_ENGINE_REL)["sha256"],
        "policy_sha256": task_validation["identities"]["policy"]["sha256"],
        "output_manifest_sha256": output_manifest_sha256,
        "error": _clean_error(error),
    }


def _write_run_metadata_and_publish(
    runtime: Any, root: Path, stage: Path, destination_rel: Path,
    *, materialized: Mapping[str, Any], manifest_row: Mapping[str, str],
    task_validation: Mapping[str, Any], fingerprint: str, attempt: str,
    materialization: Mapping[str, Any], engine_measurement: Mapping[str, Any],
    physical: int, available: int, estimated: int, admission_status: str,
    terminal_status: str, error: str, complete_engine_bundle: bool,
    engine_continuation: Mapping[str, object],
) -> None:
    effective_envelope = min(physical, int(engine_measurement["configured_envelope_bytes"]))
    if (
        not isinstance(engine_continuation.get("larger_host_continuation"), bool)
        or not isinstance(engine_continuation.get("predecessors"), list)
        or engine_continuation.get("current_effective_envelope_bytes") != effective_envelope
    ):
        raise CoordinatorError("R-engine continuation provenance is malformed")
    if not (stage / "process.log").exists():
        _write_file(runtime, stage / "process.log", bytes(engine_measurement.get("log", b"")))
    output_manifest, output_content = _output_manifest(
        runtime, root, stage, complete_engine_bundle=complete_engine_bundle,
        terminal_status=terminal_status,
    )
    _write_file(runtime, stage / "output.manifest.json", output_content)
    output_identity = runtime.stable_identity(root, stage / "output.manifest.json")
    fields = materialized["policy"]["ram_aware_execution_contract"]["resource_metrics_schema"]
    resource_row = _resource_row(
        runtime, root, task_validation, fingerprint, attempt, materialization,
        engine_measurement, physical=physical, available=available, estimated=estimated,
        admission_status=admission_status, terminal_status=terminal_status,
        output_manifest_sha256=output_identity["sha256"], error=error,
    )
    _write_file(runtime, stage / "resource_metrics.tsv", runtime.tsv_bytes(fields, [resource_row]))
    metrics_identity = runtime.stable_identity(root, stage / "resource_metrics.tsv")
    bound_output_identity = dict(output_identity)
    bound_output_identity["path"] = str(destination_rel / "output.manifest.json")
    bound_metrics_identity = dict(metrics_identity)
    bound_metrics_identity["path"] = str(destination_rel / "resource_metrics.tsv")
    measurement_without_log = {
        key: value for key, value in engine_measurement.items() if key != "log"
    }
    provenance = {
        "schema_version": RUN_SCHEMA, "analysis_id": ANALYSIS_ID,
        "artifact_role": "IMMUTABLE_PER_LOCUS_TERMINAL_RUN",
        "run_fingerprint": fingerprint, "attempt_id": attempt,
        "locus_entry_id": manifest_row["locus_entry_id"],
        "pair_id": manifest_row["pair_id"], "terminal_status": terminal_status,
        "scientific_task_unchanged": True,
        "materialized_manifest_lock": materialized["lock_identity"],
        "materialized_manifest_row_sha256": _digest(dict(manifest_row)),
        "task": task_validation["task_identity"],
        "input_identities": task_validation["identities"],
        "materialization_measurement_seal": materialization["seal_identity"],
        "materialization_measurement_receipt": materialization["receipt_identity"],
        "r_engine_measurement": measurement_without_log,
        "ram_measurement_amendment": ram_measurement_amendment(runtime, root),
        "memory_envelope_bytes": effective_envelope,
        "larger_host_continuation": engine_continuation["larger_host_continuation"],
        "engine_continuation_predecessors": engine_continuation["predecessors"],
        "minimum_host_reserve_bytes": RESERVE_BYTES,
        "output_manifest": bound_output_identity, "resource_metrics": bound_metrics_identity,
        "scalar_n_convention": runtime.SCALAR_N_CONVENTION,
        "execution_amendment": runtime.EXECUTION_AMENDMENT,
        "execution_amendment_sha256": runtime.EXECUTION_AMENDMENT_SHA256,
        "script_identities": runtime._runtime_script_identities(root),
        "model_or_prior_retry": False, "whole_locus_subsetting": False,
        "error": _clean_error(error),
    }
    provenance_content = json.dumps(provenance, indent=2, sort_keys=True).encode() + b"\n"
    _write_file(runtime, stage / "run.provenance.json", provenance_content)
    records: dict[str, dict[str, object]] = {}
    for path in sorted(stage.iterdir(), key=lambda item: item.name):
        identity = runtime.stable_identity(root, path, allow_empty=(path.name == "process.log"))
        records[path.name] = {"bytes": identity["bytes"], "sha256": identity["sha256"]}
    seal = {
        "schema_version": RUN_SCHEMA, "analysis_id": ANALYSIS_ID,
        "run_fingerprint": fingerprint, "attempt_id": attempt,
        "locus_entry_id": manifest_row["locus_entry_id"],
        "terminal_status": terminal_status, "records": records,
    }
    _write_file(runtime, stage / "run.seal.json", json.dumps(
        seal, indent=2, sort_keys=True,
    ).encode() + b"\n")
    destination_parent = runtime.ensure_directory(root, destination_rel.parent)
    publish_directory_no_replace(runtime, root, stage, destination_parent / destination_rel.name)


def _iter_run_directories(runtime: Any, root: Path, locus: str) -> list[Path]:
    parent = runtime.safe_path(root, RUN_ROOT_REL / locus, "locus run family")
    if not parent.exists():
        return []
    observed = os.lstat(parent)
    if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
        raise CoordinatorError("locus run family is not a real directory")
    output: list[Path] = []
    for child in sorted(parent.iterdir(), key=lambda item: item.name):
        status = os.lstat(child)
        if not stat.S_ISDIR(status.st_mode) or stat.S_ISLNK(status.st_mode):
            raise CoordinatorError("locus run family contains a non-directory entry")
        output.append(child)
    return output


def _bound_engine_predecessors(
    runtime: Any, root: Path, materialized: Mapping[str, Any], run_rel: Path,
    *, locus: str, fingerprint: str,
) -> list[dict[str, object]]:
    """Read and bind exactly the resource terminals preceding one attempt ID."""

    current_name = run_rel.name
    fields = materialized["policy"]["ram_aware_execution_contract"]["resource_metrics_schema"]
    output: list[dict[str, object]] = []
    for directory in _iter_run_directories(runtime, root, locus):
        if directory.name >= current_name:
            continue
        relative = directory.relative_to(root)
        seal, seal_identity = runtime.read_json(root, relative / "run.seal.json")
        provenance, provenance_identity = runtime.read_json(
            root, relative / "run.provenance.json",
        )
        _, metric_rows, resource_identity = runtime.read_tsv(
            root, relative / "resource_metrics.tsv", fields,
        )
        records = seal.get("records")
        terminal = seal.get("terminal_status")
        if (
            seal.get("schema_version") != RUN_SCHEMA or seal.get("analysis_id") != ANALYSIS_ID
            or seal.get("locus_entry_id") != locus or seal.get("attempt_id") != directory.name
            or seal.get("run_fingerprint") != fingerprint
            or terminal not in {"BLOCKED_BY_COMPUTE", "FAILED_RESOURCE_OOM"}
            or not isinstance(records, Mapping)
            or len(metric_rows) != 1
        ):
            raise CoordinatorError("R-engine continuation predecessor seal is invalid")
        for name, identity in (
            ("run.provenance.json", provenance_identity),
            ("resource_metrics.tsv", resource_identity),
        ):
            expected = records.get(name)
            if (
                not isinstance(expected, Mapping)
                or expected.get("bytes") != identity["bytes"]
                or expected.get("sha256") != identity["sha256"]
            ):
                raise CoordinatorError("R-engine predecessor evidence differs from its seal")
        metric = metric_rows[0]
        measurement = provenance.get("r_engine_measurement")
        try:
            host_physical = int(metric["host_physical_memory_bytes"])
            effective = int(provenance["memory_envelope_bytes"])
            configured = int(measurement["configured_envelope_bytes"])
        except (KeyError, TypeError, ValueError) as error:
            raise CoordinatorError("R-engine predecessor memory receipt is malformed") from error
        if (
            provenance.get("run_fingerprint") != fingerprint
            or provenance.get("attempt_id") != directory.name
            or provenance.get("locus_entry_id") != locus
            or provenance.get("terminal_status") != terminal
            or metric.get("run_fingerprint") != fingerprint
            or metric.get("attempt_id") != directory.name
            or metric.get("terminal_status") != terminal
            or effective != min(host_physical, configured)
        ):
            raise CoordinatorError("R-engine predecessor provenance/resource receipt drifted")
        output.append({
            "attempt_id": directory.name, "terminal_status": str(terminal),
            "run_directory": str(relative), "run_seal": seal_identity,
            "run_provenance": provenance_identity, "resource_metrics": resource_identity,
            "host_physical_memory_bytes": host_physical,
            "effective_envelope_bytes": effective,
        })
    return output


def _verify_engine_continuation_provenance(
    runtime: Any, root: Path, materialized: Mapping[str, Any], run_rel: Path,
    *, locus: str, fingerprint: str, provenance: Mapping[str, Any],
    measurement: Mapping[str, Any], metric: Mapping[str, str],
) -> None:
    """Prove initial 8-GiB execution or a genuine strictly larger continuation."""

    predecessors = _bound_engine_predecessors(
        runtime, root, materialized, run_rel, locus=locus, fingerprint=fingerprint,
    )
    requested = provenance.get("larger_host_continuation")
    if not isinstance(requested, bool) or provenance.get(
        "engine_continuation_predecessors"
    ) != predecessors:
        raise CoordinatorError("R-engine continuation lineage differs from predecessor seals")
    try:
        physical = int(metric["host_physical_memory_bytes"])
        configured = int(measurement["configured_envelope_bytes"])
        effective = int(provenance["memory_envelope_bytes"])
    except (KeyError, TypeError, ValueError) as error:
        raise CoordinatorError("R-engine current memory receipt is malformed") from error
    if effective != min(physical, configured):
        raise CoordinatorError("R-engine effective envelope differs from its exact receipt")
    if not predecessors:
        if requested or configured != PHYSICAL_ENVELOPE_BYTES:
            raise CoordinatorError("first R-engine attempt must use the baseline 8-GiB envelope")
        return
    if (
        not requested or configured <= PHYSICAL_ENVELOPE_BYTES
        or physical <= max(int(item["host_physical_memory_bytes"]) for item in predecessors)
        or effective <= max(int(item["effective_envelope_bytes"]) for item in predecessors)
    ):
        raise CoordinatorError(
            "R-engine continuation is not strictly larger than every resource predecessor"
        )


def verify_sealed_run(
    root: Path, run_rel: Path | str, *, materialized: Mapping[str, Any] | None = None,
    runtime: Any | None = None,
) -> dict[str, Any]:
    """Exhaustively revalidate one immutable terminal per-locus attempt."""

    runtime = runtime or _load_runtime(root)
    materialized = materialized or runtime.validate_materialized_manifest(root)
    run_rel = Path(run_rel)
    directory = runtime.safe_path(root, run_rel, "sealed run", must_exist=True)
    if not stat.S_ISDIR(os.lstat(directory).st_mode) or directory.is_symlink():
        raise CoordinatorError("sealed run is not a real directory")
    seal, seal_identity = runtime.read_json(root, run_rel / "run.seal.json")
    required_seal = {
        "schema_version", "analysis_id", "run_fingerprint", "attempt_id",
        "locus_entry_id", "terminal_status", "records",
    }
    if (
        set(seal) != required_seal or seal.get("schema_version") != RUN_SCHEMA
        or seal.get("analysis_id") != ANALYSIS_ID or seal.get("attempt_id") != run_rel.name
    ):
        raise CoordinatorError("run seal schema/identity drifted")
    terminal = seal.get("terminal_status")
    allowed = set(materialized["policy"]["ram_aware_execution_contract"]["terminal_statuses"])
    if terminal not in allowed:
        raise CoordinatorError("run terminal status is outside policy")
    records = seal.get("records")
    if not isinstance(records, dict) or not records:
        raise CoordinatorError("run seal lacks file records")
    names = {path.name for path in directory.iterdir()}
    if names != set(records) | {"run.seal.json"}:
        raise CoordinatorError("sealed run file family drifted")
    identities: dict[str, dict[str, object]] = {}
    for name, expected in records.items():
        if not isinstance(expected, dict) or set(expected) != {"bytes", "sha256"}:
            raise CoordinatorError("run record shape drifted")
        identity = runtime.stable_identity(
            root, run_rel / name, allow_empty=(name == "process.log"),
        )
        if identity["bytes"] != expected["bytes"] or identity["sha256"] != expected["sha256"]:
            raise CoordinatorError(f"sealed run artifact drifted: {name}")
        identities[name] = identity
    provenance, provenance_identity = runtime.read_json(root, run_rel / "run.provenance.json")
    locus = str(seal["locus_entry_id"])
    matches = [row for row in materialized["rows"] if row["locus_entry_id"] == locus]
    if len(matches) != 1:
        raise CoordinatorError("sealed run locus is outside materialized family")
    row = matches[0]
    task_validation = runtime.deep_validate_locus_task(
        root, row["task_path"], expected_manifest_row=row,
    )
    fingerprint = run_fingerprint(runtime, materialized, row, task_validation, root)
    if seal.get("run_fingerprint") != fingerprint:
        raise CoordinatorError("sealed run fingerprint differs from unchanged task")
    materialization = verify_materialization_measurement(
        root, locus, family=materialized["pre_family"], runtime=runtime,
    )
    expected_provenance = {
        "schema_version": RUN_SCHEMA, "analysis_id": ANALYSIS_ID,
        "artifact_role": "IMMUTABLE_PER_LOCUS_TERMINAL_RUN",
        "run_fingerprint": fingerprint, "attempt_id": seal["attempt_id"],
        "locus_entry_id": locus, "pair_id": row["pair_id"],
        "terminal_status": terminal, "scientific_task_unchanged": True,
        "materialized_manifest_lock": materialized["lock_identity"],
        "materialized_manifest_row_sha256": _digest(dict(row)),
        "task": task_validation["task_identity"],
        "input_identities": task_validation["identities"],
        "materialization_measurement_seal": materialization["seal_identity"],
        "materialization_measurement_receipt": materialization["receipt_identity"],
        "r_engine_measurement": provenance.get("r_engine_measurement"),
        "ram_measurement_amendment": ram_measurement_amendment(runtime, root),
        "memory_envelope_bytes": provenance.get("memory_envelope_bytes"),
        "larger_host_continuation": provenance.get("larger_host_continuation"),
        "engine_continuation_predecessors": provenance.get(
            "engine_continuation_predecessors"
        ),
        "minimum_host_reserve_bytes": RESERVE_BYTES,
        "output_manifest": identities["output.manifest.json"],
        "resource_metrics": identities["resource_metrics.tsv"],
        "scalar_n_convention": runtime.SCALAR_N_CONVENTION,
        "execution_amendment": runtime.EXECUTION_AMENDMENT,
        "execution_amendment_sha256": runtime.EXECUTION_AMENDMENT_SHA256,
        "script_identities": runtime._runtime_script_identities(root),
        "model_or_prior_retry": False, "whole_locus_subsetting": False,
        "error": provenance.get("error"),
    }
    if provenance != expected_provenance:
        raise CoordinatorError("run provenance differs from its exact bound family")
    measurement = provenance["r_engine_measurement"]
    if not isinstance(measurement, dict) or "log" in measurement:
        raise CoordinatorError("R measurement provenance is malformed")
    output_manifest, _ = runtime.read_json(root, run_rel / "output.manifest.json")
    if (
        output_manifest.get("terminal_status") != terminal
        or output_manifest.get("analysis_id") != ANALYSIS_ID
        or not isinstance(output_manifest.get("artifacts"), dict)
    ):
        raise CoordinatorError("run output manifest drifted")
    for name, expected in output_manifest["artifacts"].items():
        identity = identities.get(name)
        if identity is None or identity["bytes"] != expected.get("bytes") or identity["sha256"] != expected.get("sha256"):
            raise CoordinatorError("run output-manifest member drifted")
    metric_fields = materialized["policy"]["ram_aware_execution_contract"]["resource_metrics_schema"]
    _, metric_rows, resource_identity = runtime.read_tsv(
        root, run_rel / "resource_metrics.tsv", metric_fields,
    )
    if len(metric_rows) != 1:
        raise CoordinatorError("sealed run must contain one resource row")
    metric = metric_rows[0]
    if (
        metric["analysis_id"] != ANALYSIS_ID or metric["run_fingerprint"] != fingerprint
        or metric["locus_entry_id"] != locus or metric["pair_id"] != row["pair_id"]
        or metric["attempt_id"] != str(seal["attempt_id"])
        or metric["terminal_status"] != terminal
        or metric["task_sha256"] != task_validation["task_identity"]["sha256"]
        or metric["output_manifest_sha256"] != identities["output.manifest.json"]["sha256"]
    ):
        raise CoordinatorError("sealed resource metric identity drifted")
    _verify_engine_continuation_provenance(
        runtime, root, materialized, run_rel, locus=locus, fingerprint=fingerprint,
        provenance=provenance, measurement=measurement, metric=metric,
    )
    engine_bundle: dict[str, Any] | None = None
    if terminal == "COMPLETE_SEALED":
        if (
            output_manifest.get("complete_engine_bundle") is not True
            or measurement.get("rss_measurement_backend") != PGID_BACKEND
            or int(measurement.get("positive_pgid_samples", 0)) <= 0
            or int(measurement.get("peak_rss_bytes", 0)) <= 0
            or measurement.get("exit_code") != 0 or measurement.get("guard_failure") is not None
            or metric["rss_measurement_backend"] != PGID_BACKEND
        ):
            raise CoordinatorError("successful run lacks exact positive PGID measurement")
        engine_bundle = runtime.validate_engine_bundle(root, run_rel, task_validation)
    elif output_manifest.get("complete_engine_bundle") is not False:
        raise CoordinatorError("failed run falsely claims a complete engine bundle")
    return {
        "seal": seal, "seal_identity": seal_identity, "provenance": provenance,
        "provenance_identity": provenance_identity, "resource_identity": resource_identity,
        "terminal_status": terminal, "row": row, "task": task_validation,
        "materialization": materialization, "resource_row": metric,
        "output_manifest": output_manifest, "engine_bundle": engine_bundle,
        "run_rel": str(run_rel),
    }


def _quarantine_interrupted_run_stages(
    runtime: Any, root: Path, locus: str,
) -> list[dict[str, object]]:
    parent = runtime.safe_path(root, STAGE_ROOT_REL / locus, "run staging family")
    if not parent.exists():
        return []
    if not stat.S_ISDIR(os.lstat(parent).st_mode) or parent.is_symlink():
        raise CoordinatorError("run staging family is unsafe")
    receipts: list[dict[str, object]] = []
    for child in sorted(parent.iterdir(), key=lambda item: item.name):
        receipts.append(_quarantine_directory(
            runtime, root, child, category="interrupted_r_engine_attempt",
            locus_entry_id=locus, reason="UNSEALED_STAGE_AFTER_PAIR_LEASE_BECAME_AVAILABLE",
        ))
    return receipts


def _existing_runs(
    runtime: Any, root: Path, materialized: Mapping[str, Any], locus: str,
) -> list[dict[str, Any]]:
    return [
        verify_sealed_run(
            root, path.relative_to(root), materialized=materialized, runtime=runtime,
        )
        for path in _iter_run_directories(runtime, root, locus)
    ]


def _engine_predecessor_records(
    prior_runs: Sequence[Mapping[str, Any]],
) -> list[dict[str, object]]:
    """Bind every prior resource terminal by attempt and immutable file identities."""

    records: list[dict[str, object]] = []
    for run in prior_runs:
        if run.get("terminal_status") not in {"BLOCKED_BY_COMPUTE", "FAILED_RESOURCE_OOM"}:
            raise CoordinatorError("engine continuation predecessor is not a resource terminal")
        provenance = run["provenance"]
        metric = run["resource_row"]
        records.append({
            "attempt_id": str(run["seal"]["attempt_id"]),
            "terminal_status": str(run["terminal_status"]),
            "run_directory": str(run["run_rel"]),
            "run_seal": run["seal_identity"],
            "run_provenance": run["provenance_identity"],
            "resource_metrics": run["resource_identity"],
            "host_physical_memory_bytes": int(metric["host_physical_memory_bytes"]),
            "effective_envelope_bytes": int(provenance["memory_envelope_bytes"]),
        })
    return records


def validate_larger_host_engine_continuation(
    prior_runs: Sequence[Mapping[str, Any]], *, current_physical_bytes: int,
    configured_envelope_bytes: int,
) -> dict[str, object]:
    """Require a genuine strictly larger retry of retained resource-only terminals."""

    if not prior_runs:
        raise CoordinatorError("R-engine larger-host continuation requires a retained resource failure")
    predecessors = _engine_predecessor_records(prior_runs)
    current_physical = int(current_physical_bytes)
    configured = int(configured_envelope_bytes)
    current_effective = min(current_physical, configured)
    if (
        configured <= PHYSICAL_ENVELOPE_BYTES
        or current_physical <= max(int(item["host_physical_memory_bytes"]) for item in predecessors)
        or current_effective <= max(int(item["effective_envelope_bytes"]) for item in predecessors)
    ):
        raise CoordinatorError(
            "R-engine retry requires a strictly larger host and effective envelope"
        )
    return {
        "larger_host_continuation": True,
        "predecessors": predecessors,
        "current_effective_envelope_bytes": current_effective,
    }


def execute_locus(
    runtime: Any, root: Path, materialized: Mapping[str, Any],
    manifest_row: Mapping[str, str], *, configured_envelope_bytes: int = PHYSICAL_ENVELOPE_BYTES,
    larger_host_continuation: bool = False,
) -> dict[str, Any]:
    """Execute or resume one locus while preserving every prior attempt."""

    locus = manifest_row["locus_entry_id"]
    prior = _existing_runs(runtime, root, materialized, locus)
    complete = [run for run in prior if run["terminal_status"] == "COMPLETE_SEALED"]
    if len(complete) > 1:
        raise CoordinatorError("more than one successful run exists for one frozen locus")
    if complete:
        return complete[0]
    if larger_host_continuation and not prior:
        raise CoordinatorError(
            "R-engine larger-host continuation requires a retained resource failure"
        )
    if prior:
        nonresource = [
            run for run in prior
            if run["terminal_status"] not in {"BLOCKED_BY_COMPUTE", "FAILED_RESOURCE_OOM"}
        ]
        if nonresource:
            # Invalid/engine failures are retained and never retried with a
            # parameter alteration or silently interpreted as null evidence.
            return nonresource[-1]
        if not larger_host_continuation:
            return prior[-1]
    task_validation = runtime.deep_validate_locus_task(
        root, manifest_row["task_path"], expected_manifest_row=manifest_row,
    )
    fingerprint = run_fingerprint(runtime, materialized, manifest_row, task_validation, root)
    materialization = verify_materialization_measurement(
        root, locus, family=materialized["pre_family"], runtime=runtime,
    )
    attempt = _attempt_id()
    with worker_lease(runtime, root, manifest_row["pair_id"], fingerprint) as lease_fd:
        locked_prior = _existing_runs(runtime, root, materialized, locus)
        if [run["seal_identity"] for run in locked_prior] != [
            run["seal_identity"] for run in prior
        ]:
            raise CoordinatorError("R-engine predecessor family changed while acquiring pair lease")
        prior = locked_prior
        _quarantine_interrupted_run_stages(runtime, root, locus)
        # Revalidate all immutable inputs after acquiring the pair-stable lease
        # and immediately before any result process is launched.
        task_validation = runtime.deep_validate_locus_task(
            root, manifest_row["task_path"], expected_manifest_row=manifest_row,
        )
        if run_fingerprint(runtime, materialized, manifest_row, task_validation, root) != fingerprint:
            raise CoordinatorError("task fingerprint changed while acquiring worker lease")
        physical = runtime.host_physical_memory_bytes()
        available = runtime._memory_available_bytes()
        effective_envelope = min(physical, int(configured_envelope_bytes))
        if configured_envelope_bytes < PHYSICAL_ENVELOPE_BYTES:
            raise CoordinatorError("configured memory envelope is below the frozen 8-GiB envelope")
        if configured_envelope_bytes > PHYSICAL_ENVELOPE_BYTES and not larger_host_continuation:
            raise CoordinatorError("larger memory envelope requires explicit larger-host continuation")
        continuation = (
            validate_larger_host_engine_continuation(
                prior, current_physical_bytes=physical,
                configured_envelope_bytes=configured_envelope_bytes,
            )
            if larger_host_continuation
            else {
                "larger_host_continuation": False, "predecessors": [],
                "current_effective_envelope_bytes": effective_envelope,
            }
        )
        estimated = runtime.estimate_peak_rss_bytes(int(manifest_row["variant_count"]), materialized["policy"])
        admitted = available >= estimated + RESERVE_BYTES and estimated + RESERVE_BYTES <= effective_envelope
        stage = _new_stage(runtime, root, STAGE_ROOT_REL / locus, attempt)
        destination_rel = RUN_ROOT_REL / locus / attempt
        if not admitted:
            measurement = {
                "command": [], "started_utc": _utc_now(), "finished_utc": _utc_now(),
                "elapsed_monotonic_seconds": 0.0, "peak_rss_bytes": 0,
                "rss_measurement_backend": materialization["receipt"]["rss_measurement_backend"],
                "positive_pgid_samples": 0, "exit_code": "NA", "guard_failure": None,
                "success_measurement": False, "log": b"",
                "configured_envelope_bytes": int(configured_envelope_bytes),
            }
            error = (
                "UNCHANGED_WHOLE_LOCUS_NOT_ADMITTED: available memory or envelope cannot retain "
                "the locked 1-GiB host reserve"
            )
            _write_run_metadata_and_publish(
                runtime, root, stage, destination_rel, materialized=materialized,
                manifest_row=manifest_row, task_validation=task_validation,
                fingerprint=fingerprint, attempt=attempt, materialization=materialization,
                engine_measurement=measurement, physical=physical, available=available,
                estimated=estimated, admission_status="BLOCKED_BY_COMPUTE",
                terminal_status="BLOCKED_BY_COMPUTE", error=error,
                complete_engine_bundle=False, engine_continuation=continuation,
            )
            return verify_sealed_run(
                root, destination_rel, materialized=materialized, runtime=runtime,
            )

        command = [
            str(root / RSCRIPT_REL), "--vanilla", str(root / R_ENGINE_REL),
            str(root / manifest_row["task_path"]), str(stage),
        ]
        measurement = monitored_command(
            command, root=root, inherited_fds=(lease_fd,),
            maximum_owned_pgid_rss_bytes=effective_envelope - RESERVE_BYTES,
            extra_environment={"TRACK_B_FINEMAP_WORKER_LEASE_FD": str(lease_fd)},
            memory_available=runtime._memory_available_bytes,
        )
        measurement["configured_envelope_bytes"] = int(configured_envelope_bytes)
        _write_file(runtime, stage / "process.log", measurement["log"])
        error = measurement["guard_failure"] or "NA"
        terminal = "FAILED_ENGINE"
        complete_bundle = False
        if (
            str(measurement["guard_failure"] or "").startswith(
                "PHYSICAL_RAM_CAP_OR_HOST_RESERVE_EXCEEDED"
            )
            or measurement["exit_code"] in {-9, 137}
        ):
            terminal = "FAILED_RESOURCE_OOM"
            error = measurement["guard_failure"] or f"R process exit {measurement['exit_code']}"
        elif measurement["exit_code"] != 0:
            terminal = "FAILED_ENGINE"
            error = measurement["guard_failure"] or f"R process exit {measurement['exit_code']}"
        elif not measurement["success_measurement"]:
            terminal = "FAILED_VALIDATION"
            error = "successful R exit lacks a positive authoritative live PGID RSS measurement"
        else:
            try:
                runtime.adapt_engine_outputs(root, stage.relative_to(root), task_validation)
                runtime.validate_engine_bundle(root, stage.relative_to(root), task_validation)
            except BaseException as validation_error:
                terminal = "FAILED_VALIDATION"
                error = f"{type(validation_error).__name__}: {validation_error}"
            else:
                terminal = "COMPLETE_SEALED"
                complete_bundle = True
                error = "NA"
        _write_run_metadata_and_publish(
            runtime, root, stage, destination_rel, materialized=materialized,
            manifest_row=manifest_row, task_validation=task_validation,
            fingerprint=fingerprint, attempt=attempt, materialization=materialization,
            engine_measurement=measurement, physical=physical, available=available,
            estimated=estimated, admission_status="ADMITTED_WHOLE_LOCUS",
            terminal_status=terminal, error=error, complete_engine_bundle=complete_bundle,
            engine_continuation=continuation,
        )
    return verify_sealed_run(root, destination_rel, materialized=materialized, runtime=runtime)


def _publish_or_verify_bytes(
    runtime: Any, root: Path, relative: Path, content: bytes,
) -> dict[str, object]:
    path = runtime.safe_path(root, relative, "versioned immutable artifact")
    if path.exists() or path.is_symlink():
        identity = runtime.stable_identity(root, relative, allow_empty=(len(content) == 0))
        if identity["bytes"] != len(content) or identity["sha256"] != hashlib.sha256(content).hexdigest():
            raise CoordinatorError(f"versioned artifact conflicts with existing content: {relative}")
        return identity
    return runtime.publish_bytes_no_replace(root, relative, content)


def component_report_paths(report_digest: str) -> dict[str, Path]:
    if not HEX64.fullmatch(report_digest):
        raise CoordinatorError("invalid component report digest")
    stem = f"RAM_BY_LOCUS.{report_digest[:20]}"
    return {
        "report": Path("results/track_b/finemapping") / f"{stem}.tsv",
        "namespace": Path("results/track_b/finemapping") / f"{stem}.namespace.json",
        "provenance": Path("results/track_b/finemapping") / f"{stem}.provenance.json",
    }


def build_component_ram_rows(
    runtime: Any, family: Mapping[str, Any], materialization_outcomes: Mapping[str, Mapping[str, Any]],
    run_outcomes: Mapping[str, Mapping[str, Any]] | None = None,
) -> list[dict[str, str]]:
    """Build exactly one authoritative row for every mapped/unavailable entry."""

    runs = run_outcomes or {}
    output: list[dict[str, str]] = []
    for pre_row in family["rows"]:
        locus = pre_row["locus_entry_id"]
        outcome = materialization_outcomes.get(locus)
        if outcome is None:
            raise CoordinatorError("component report omits a frozen mapped locus")
        run = runs.get(locus)
        if run is not None:
            metric = run["resource_row"]
            n_snps = metric["variant_count"]
            peak = int(metric["peak_rss_bytes"])
            runtime_seconds = float(metric["elapsed_monotonic_seconds"])
            terminal = metric["terminal_status"]
            output_hash = metric["output_manifest_sha256"]
        else:
            receipt = outcome["receipt"]
            bundle = receipt.get("bundle_lock")
            n_snps = (
                str(outcome["bundle"]["manifest_row"]["variant_count"])
                if outcome.get("bundle") is not None else "NA"
            )
            peak = int(receipt.get("peak_rss_bytes", 0))
            runtime_seconds = float(receipt.get("elapsed_monotonic_seconds", 0.0))
            terminal = str(outcome["terminal_status"])
            output_hash = str(bundle.get("sha256")) if isinstance(bundle, Mapping) else str(
                outcome["seal_identity"]["sha256"]
            )
        output.append({
            "analysis": "TRACK_B_FINEMAPPING_TRAIT_COLOC_V1",
            "pair": pre_row["pair_id"], "locus": locus,
            "chromosome": pre_row["CHR"], "n_snps": n_snps,
            "peak_ram_gb": format(peak / 1024**3, ".12g"),
            "runtime_sec": format(runtime_seconds, ".12g"),
            "exit_status": terminal, "output_hash": output_hash,
        })
    for row in family["unavailable_rows"]:
        output.append({
            "analysis": "TRACK_B_FINEMAPPING_TRAIT_COLOC_V1",
            "pair": row["pair_id"], "locus": row["unavailable_entry_id"],
            "chromosome": row["CHR"], "n_snps": "0", "peak_ram_gb": "0",
            "runtime_sec": "0",
            "exit_status": "BLOCKED_BY_DATA:NO_PREDECLARED_SIGNED_LD_BLOCK",
            "output_hash": family["unavailable_identity"]["sha256"],
        })
    keys = [(row["analysis"], row["pair"], row["locus"], row["chromosome"]) for row in output]
    if len(keys) != len(set(keys)):
        raise CoordinatorError("component RAM collision key is duplicated")
    return output


def publish_component_ram_report(
    runtime: Any, root: Path, family: Mapping[str, Any],
    materialization_outcomes: Mapping[str, Mapping[str, Any]],
    run_outcomes: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    rows = build_component_ram_rows(runtime, family, materialization_outcomes, run_outcomes)
    report_content = runtime.tsv_bytes(COMPONENT_RAM_FIELDS, rows)
    report_digest = _digest({
        "pre_family_lock": family["lock_identity"], "rows_sha256": hashlib.sha256(report_content).hexdigest(),
        "schema": COMPONENT_RAM_FIELDS,
    })
    paths = component_report_paths(report_digest)
    report_expected = {
        "path": str(paths["report"]), "bytes": len(report_content),
        "sha256": hashlib.sha256(report_content).hexdigest(),
    }
    namespace = {
        "schema_version": "sleep-atlas-component-ram-namespace.1",
        "component": "TRACK_B_FINEMAPPING_TRAIT_COLOC_V1",
        "analysis_id": ANALYSIS_ID, "table": str(paths["report"]),
        "exact_fields": COMPONENT_RAM_FIELDS,
        "collision_key": ["analysis", "pair", "locus", "chromosome"],
        "authoritative_peak_rule": (
            "MAX_MATERIALIZATION_AND_R_ENGINE_DEDUPLICATED_COORDINATOR_PARENT_PLUS_"
            "ALL_OWNED_PGID_MEMBER_RSS_FROM_ONE_SNAPSHOT"
        ),
        "runtime_rule": "MATERIALIZATION_MONOTONIC_SECONDS_PLUS_R_ENGINE_MONOTONIC_SECONDS",
        "successful_measurement_backend": PGID_BACKEND,
        "pre_family_lock": family["lock_identity"],
    }
    namespace_content = json.dumps(namespace, indent=2, sort_keys=True).encode() + b"\n"
    namespace_expected = {
        "path": str(paths["namespace"]), "bytes": len(namespace_content),
        "sha256": hashlib.sha256(namespace_content).hexdigest(),
    }
    provenance = {
        "schema_version": "sleep-atlas-track-b-finemapping-component-ram.1",
        "analysis_id": ANALYSIS_ID, "artifact_role": "VERSIONED_COMPONENT_RAM_FEDERATION_INPUT",
        "report_digest": report_digest, "report": report_expected,
        "namespace": namespace_expected, "row_count": len(rows),
        "collision_key_unique": True,
        "mapped_locus_count": len(family["rows"]),
        "unavailable_locus_count": len(family["unavailable_rows"]),
        "pre_family_lock": family["lock_identity"],
        "ram_measurement_amendment": ram_measurement_amendment(runtime, root),
        "terminal_status_counts": dict(sorted(
            __import__("collections").Counter(row["exit_status"] for row in rows).items()
        )),
    }
    provenance_content = json.dumps(provenance, indent=2, sort_keys=True).encode() + b"\n"
    report_identity = _publish_or_verify_bytes(runtime, root, paths["report"], report_content)
    namespace_identity = _publish_or_verify_bytes(runtime, root, paths["namespace"], namespace_content)
    provenance_identity = _publish_or_verify_bytes(runtime, root, paths["provenance"], provenance_content)
    return {
        "paths": paths, "rows": rows, "report_identity": report_identity,
        "namespace_identity": namespace_identity, "provenance_identity": provenance_identity,
    }


def verify_component_ram_report(
    root: Path = ROOT, report_rel: Path | str | None = None, *, runtime: Any | None = None,
) -> dict[str, Any]:
    """Side-effect-free deep verifier for the nine-column federation table."""

    runtime = runtime or _load_runtime(root)
    if report_rel is None:
        directory = runtime.safe_path(root, Path("results/track_b/finemapping"), "component report root", must_exist=True)
        candidates = sorted(directory.glob("RAM_BY_LOCUS.*.tsv"))
        if len(candidates) != 1:
            raise CoordinatorError("specify one component report when zero or multiple versions exist")
        report_rel = candidates[0].relative_to(root)
    report_rel = Path(report_rel)
    if report_rel.is_absolute():
        try:
            report_rel = report_rel.relative_to(root.resolve())
        except ValueError as exc:
            raise CoordinatorError("component report is outside the repository root") from exc
    match = re.fullmatch(r"RAM_BY_LOCUS\.([0-9a-f]{20})\.tsv", report_rel.name)
    if match is None:
        raise CoordinatorError("component report filename is not versioned by its digest")
    stem = report_rel.name[:-4]
    namespace_rel = report_rel.parent / f"{stem}.namespace.json"
    provenance_rel = report_rel.parent / f"{stem}.provenance.json"
    _, rows, report_identity = runtime.read_tsv(root, report_rel, COMPONENT_RAM_FIELDS)
    keys = [(row["analysis"], row["pair"], row["locus"], row["chromosome"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise CoordinatorError("component report collision key is duplicated")
    for row in rows:
        peak = float(row["peak_ram_gb"])
        elapsed = float(row["runtime_sec"])
        if not math.isfinite(peak) or peak < 0 or not math.isfinite(elapsed) or elapsed < 0:
            raise CoordinatorError("component report has invalid resource values")
        if row["n_snps"] != "NA" and (not row["n_snps"].isdigit() or int(row["n_snps"]) < 0):
            raise CoordinatorError("component report has invalid SNP count")
        if not HEX64.fullmatch(row["output_hash"]):
            raise CoordinatorError("component report output hash is invalid")
    namespace, namespace_identity = runtime.read_json(root, namespace_rel)
    provenance, provenance_identity = runtime.read_json(root, provenance_rel)
    expected_namespace = {
        "schema_version": "sleep-atlas-component-ram-namespace.1",
        "component": "TRACK_B_FINEMAPPING_TRAIT_COLOC_V1",
        "analysis_id": ANALYSIS_ID, "table": str(report_rel),
        "exact_fields": COMPONENT_RAM_FIELDS,
        "collision_key": ["analysis", "pair", "locus", "chromosome"],
        "authoritative_peak_rule": (
            "MAX_MATERIALIZATION_AND_R_ENGINE_DEDUPLICATED_COORDINATOR_PARENT_PLUS_"
            "ALL_OWNED_PGID_MEMBER_RSS_FROM_ONE_SNAPSHOT"
        ),
        "runtime_rule": "MATERIALIZATION_MONOTONIC_SECONDS_PLUS_R_ENGINE_MONOTONIC_SECONDS",
        "successful_measurement_backend": PGID_BACKEND,
        "pre_family_lock": provenance.get("pre_family_lock"),
    }
    if namespace != expected_namespace:
        raise CoordinatorError("component RAM namespace drifted")
    report_digest = provenance.get("report_digest")
    if not isinstance(report_digest, str) or not HEX64.fullmatch(report_digest) or report_digest[:20] != match.group(1):
        raise CoordinatorError("component report digest/name drifted")
    if report_digest != _digest({
        "pre_family_lock": provenance.get("pre_family_lock"),
        "rows_sha256": report_identity["sha256"], "schema": COMPONENT_RAM_FIELDS,
    }):
        raise CoordinatorError("component report content digest drifted")
    if (
        provenance.get("schema_version") != "sleep-atlas-track-b-finemapping-component-ram.1"
        or provenance.get("analysis_id") != ANALYSIS_ID
        or provenance.get("report") != report_identity
        or provenance.get("namespace") != namespace_identity
        or provenance.get("row_count") != len(rows)
        or provenance.get("collision_key_unique") is not True
        or provenance.get("ram_measurement_amendment") != ram_measurement_amendment(runtime, root)
    ):
        raise CoordinatorError("component RAM provenance drifted")
    expected_counts = dict(sorted(__import__("collections").Counter(
        row["exit_status"] for row in rows
    ).items()))
    if provenance.get("terminal_status_counts") != expected_counts:
        raise CoordinatorError("component report terminal counts drifted")
    return {
        "state": "VERIFIED_COMPONENT_RAM_REPORT",
        "rows": rows, "report_identity": report_identity,
        "namespace_identity": namespace_identity, "provenance_identity": provenance_identity,
        "collision_key": ["analysis", "pair", "locus", "chromosome"],
        "paths": {
            "report": str(report_rel), "namespace": str(namespace_rel),
            "provenance": str(provenance_rel),
        },
    }


def _canonical_targets(policy: Mapping[str, Any]) -> dict[str, Path]:
    return {
        "10_finemap_trait1.tsv": Path(policy["required_science_outputs"]["output_10"]["path"]),
        "11_finemap_trait2.tsv": Path(policy["required_science_outputs"]["output_11"]["path"]),
        "12_trait_trait_coloc.tsv": Path(policy["required_science_outputs"]["output_12"]["path"]),
        "credible_sets.tsv": Path("results/track_b/finemapping/credible_sets.tsv"),
        "finemap_qc.tsv": Path("results/track_b/finemapping/finemap_qc.tsv"),
        "diagnostics.tsv": Path("results/track_b/finemapping/diagnostics.tsv"),
        "resource_metrics.tsv": RESOURCE_METRICS_REL,
        "run_index.tsv": RUN_INDEX_REL,
        "results.provenance.json": RESULT_PROVENANCE_REL,
    }


def _complete_runs(
    runtime: Any, root: Path, materialized: Mapping[str, Any],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in materialized["rows"]:
        attempts = _existing_runs(runtime, root, materialized, row["locus_entry_id"])
        successful = [attempt for attempt in attempts if attempt["terminal_status"] == "COMPLETE_SEALED"]
        if len(successful) != 1:
            raise CoordinatorError(
                f"canonical family requires exactly one complete sealed run for {row['locus_entry_id']}"
            )
        output.append(successful[0])
    if not output:
        raise CoordinatorError("nonzero family cannot publish header-only science")
    return output


def _canonical_rows(
    runtime: Any, materialized: Mapping[str, Any], runs: Sequence[Mapping[str, Any]],
) -> dict[str, list[dict[str, str]]]:
    output: dict[str, list[dict[str, str]]] = {
        "10_finemap_trait1.tsv": [], "11_finemap_trait2.tsv": [],
        "12_trait_trait_coloc.tsv": [], "credible_sets.tsv": [],
        "finemap_qc.tsv": [], "diagnostics.tsv": [], "resource_metrics.tsv": [],
        "run_index.tsv": [],
    }
    if len(runs) != len(materialized["rows"]):
        raise CoordinatorError("canonical run family count differs from materialized loci")
    for manifest_row, run in zip(materialized["rows"], runs, strict=True):
        if run["row"] != manifest_row or run["engine_bundle"] is None:
            raise CoordinatorError("canonical run order/family drifted")
        rows = run["engine_bundle"]["rows"]
        output["10_finemap_trait1.tsv"].extend(rows["10_finemap_trait1.fragment.tsv"])
        output["11_finemap_trait2.tsv"].extend(rows["11_finemap_trait2.fragment.tsv"])
        output["12_trait_trait_coloc.tsv"].extend(rows["12_trait_trait_coloc.fragment.tsv"])
        output["credible_sets.tsv"].extend(rows["credible_sets.tsv"])
        output["finemap_qc.tsv"].extend(rows["finemap_qc.tsv"])
        output["diagnostics.tsv"].extend(rows["diagnostics.tsv"])
        output["resource_metrics.tsv"].append(run["resource_row"])
        output["run_index.tsv"].append({
            "analysis_id": ANALYSIS_ID,
            "run_fingerprint": run["seal"]["run_fingerprint"],
            "locus_entry_id": manifest_row["locus_entry_id"],
            "pair_id": manifest_row["pair_id"], "CHR": manifest_row["CHR"],
            "ld_block_id": manifest_row["ld_block_id"],
            "variant_count": manifest_row["variant_count"],
            "terminal_status": run["terminal_status"],
            "run_directory": run["run_rel"],
            "run_seal_sha256": run["seal_identity"]["sha256"],
            "output_manifest_sha256": run["resource_row"]["output_manifest_sha256"],
            "error": run["resource_row"]["error"],
        })
    expected_loci = [row["locus_entry_id"] for row in materialized["rows"]]
    if [row["locus_entry_id"] for row in output["run_index.tsv"]] != expected_loci:
        raise CoordinatorError("run index lost deterministic complete-family order")
    for name in ("10_finemap_trait1.tsv", "11_finemap_trait2.tsv", "12_trait_trait_coloc.tsv"):
        if not output[name]:
            raise CoordinatorError(f"nonzero family would publish header-only science: {name}")
    for role_name in ("10_finemap_trait1.tsv", "11_finemap_trait2.tsv"):
        observed = {row["locus_entry_id"] for row in output[role_name]}
        if observed != set(expected_loci):
            raise CoordinatorError("fine-map science output omits a frozen locus")
    if {row["locus_entry_id"] for row in output["12_trait_trait_coloc.tsv"]} != set(expected_loci):
        raise CoordinatorError("coloc science output omits a frozen locus")
    if len(output["diagnostics.tsv"]) != 2 * len(expected_loci) or len(output["finemap_qc.tsv"]) != 2 * len(expected_loci):
        raise CoordinatorError("diagnostic/QC family does not retain both traits per locus")
    return output


def _canonical_contents(
    runtime: Any, materialized: Mapping[str, Any], rows: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict[str, bytes]:
    policy = materialized["policy"]
    schemas = {
        "10_finemap_trait1.tsv": policy["required_science_outputs"]["output_10"]["schema"],
        "11_finemap_trait2.tsv": policy["required_science_outputs"]["output_11"]["schema"],
        "12_trait_trait_coloc.tsv": policy["required_science_outputs"]["output_12"]["schema"],
        "credible_sets.tsv": runtime.CREDIBLE_SET_FIELDS,
        "finemap_qc.tsv": runtime.FINEMAP_QC_FIELDS,
        "diagnostics.tsv": policy["diagnostic_contract"]["required_schema"],
        "resource_metrics.tsv": policy["ram_aware_execution_contract"]["resource_metrics_schema"],
        "run_index.tsv": RUN_INDEX_FIELDS,
    }
    return {name: runtime.tsv_bytes(fields, rows[name]) for name, fields in schemas.items()}


def _canonical_expected(
    runtime: Any, root: Path, materialized: Mapping[str, Any], runs: Sequence[Mapping[str, Any]],
    component_report: Mapping[str, Any],
) -> tuple[dict[str, bytes], dict[str, Any]]:
    rows = _canonical_rows(runtime, materialized, runs)
    contents = _canonical_contents(runtime, materialized, rows)
    targets = _canonical_targets(materialized["policy"])
    identities = {
        name: {
            "path": str(targets[name]), "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
        }
        for name, content in contents.items()
    }
    classifications: dict[str, int] = {}
    for row in rows["12_trait_trait_coloc.tsv"]:
        classifications[row["classification"]] = classifications.get(row["classification"], 0) + 1
    fine_statuses: dict[str, int] = {}
    for name in ("10_finemap_trait1.tsv", "11_finemap_trait2.tsv"):
        for row in rows[name]:
            fine_statuses[row["analysis_status"]] = fine_statuses.get(row["analysis_status"], 0) + 1
    provenance = {
        "schema_version": PUBLICATION_SCHEMA, "analysis_id": ANALYSIS_ID,
        "artifact_role": "CANONICAL_COMPLETE_TRACK_B_FINE_MAPPING_AND_TRAIT_COLOC_FAMILY",
        "materialized_manifest": materialized["manifest_identity"],
        "materialized_manifest_lock": materialized["lock_identity"],
        "pre_family_manifest": materialized["pre_family"]["manifest_identity"],
        "pre_family_lock": materialized["pre_family"]["lock_identity"],
        "policy": materialized["pre_family"]["policy_identity"],
        "locus_count": len(materialized["rows"]),
        "locus_entry_ids_in_order": [row["locus_entry_id"] for row in materialized["rows"]],
        "unavailable_locus_count": len(materialized["pre_family"]["unavailable_rows"]),
        "unavailable_manifest": materialized["pre_family"]["unavailable_identity"],
        "unavailable_terminal_entries": [
            {
                "unavailable_entry_id": row["unavailable_entry_id"],
                "pair_id": row["pair_id"], "CHR": row["CHR"], "BP": row["BP"],
                "terminal_status": "BLOCKED_BY_DATA",
                "blocker_reason": "NO_PREDECLARED_SIGNED_LD_BLOCK",
                "evidence_ids": row["evidence_ids"],
            }
            for row in materialized["pre_family"]["unavailable_rows"]
        ],
        "upstream_evidence_denominator": materialized["pre_family"]["lock"]["upstream_evidence_count"],
        "mapped_evidence_count": materialized["pre_family"]["lock"]["evidence_count"],
        "unavailable_evidence_count": materialized["pre_family"]["lock"]["unavailable_evidence_count"],
        "run_seals": {run["row"]["locus_entry_id"]: run["seal_identity"] for run in runs},
        "outputs": identities,
        "component_ram_report": component_report["report_identity"],
        "component_ram_namespace": component_report["namespace_identity"],
        "component_ram_provenance": component_report["provenance_identity"],
        "trait_coloc_classification_counts": dict(sorted(classifications.items())),
        "fine_mapping_status_counts": dict(sorted(fine_statuses.items())),
        "scalar_n_convention": runtime.SCALAR_N_CONVENTION,
        "rss_ld_s_policy": {
            "hard_invalidity_threshold": None,
            "warning_trigger": "s>0.10",
            "warning_label": "LD_UNCERTAINTY",
            "strong_claim_allowed_when_warning_present": False,
            "result_row_retained": True,
        },
        "kriging_outlier_rule": "logLR>2 AND abs(z)>2",
        "ram_measurement_amendment": ram_measurement_amendment(runtime, root),
        "complete_family_validated_before_publication": True,
        "header_only_science_outputs": False,
        "publication_commit_marker": str(RESULT_PROVENANCE_REL),
        "script_identities": runtime._runtime_script_identities(root),
    }
    provenance_content = json.dumps(provenance, indent=2, sort_keys=True).encode() + b"\n"
    contents["results.provenance.json"] = provenance_content
    return contents, provenance


def _package_manifest(
    runtime: Any, contents: Mapping[str, bytes], targets: Mapping[str, Path],
    fingerprint: str,
) -> tuple[dict[str, Any], bytes]:
    artifacts = {
        name: {
            "canonical_path": str(targets[name]), "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
        }
        for name, content in sorted(contents.items())
    }
    payload = {
        "schema_version": PUBLICATION_SCHEMA, "analysis_id": ANALYSIS_ID,
        "publication_fingerprint": fingerprint, "artifacts": artifacts,
        "commit_marker": "results.provenance.json",
        "exclusive_no_replace": True,
    }
    return payload, json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n"


def _verify_publication_package(
    runtime: Any, root: Path, package_rel: Path, contents: Mapping[str, bytes],
    targets: Mapping[str, Path], fingerprint: str,
) -> dict[str, Any]:
    package = runtime.safe_path(root, package_rel, "publication package", must_exist=True)
    if not stat.S_ISDIR(os.lstat(package).st_mode) or package.is_symlink():
        raise CoordinatorError("publication package is unsafe")
    manifest, manifest_identity = runtime.read_json(root, package_rel / "package.manifest.json")
    expected_manifest, _ = _package_manifest(runtime, contents, targets, fingerprint)
    if manifest != expected_manifest:
        raise CoordinatorError("publication package manifest drifted")
    names = {path.name for path in package.iterdir()}
    if names != set(contents) | {"package.manifest.json"}:
        raise CoordinatorError("publication package file family drifted")
    identities: dict[str, dict[str, object]] = {}
    for name, content in contents.items():
        identity = runtime.stable_identity(root, package_rel / name)
        if identity["bytes"] != len(content) or identity["sha256"] != hashlib.sha256(content).hexdigest():
            raise CoordinatorError(f"publication package artifact drifted: {name}")
        identities[name] = identity
    return {"manifest": manifest, "manifest_identity": manifest_identity, "identities": identities}


def publish_canonical_family(
    runtime: Any, root: Path, materialized: Mapping[str, Any], runs: Sequence[Mapping[str, Any]],
    component_report: Mapping[str, Any],
) -> dict[str, Any]:
    """Stage, atomically seal, and exclusively promote the complete family."""

    contents, provenance = _canonical_expected(runtime, root, materialized, runs, component_report)
    targets = _canonical_targets(materialized["policy"])
    fingerprint = _digest({
        "schema_version": PUBLICATION_SCHEMA,
        "materialized_lock": materialized["lock_identity"],
        "run_seals": [run["seal_identity"] for run in runs],
        "output_hashes": {name: hashlib.sha256(content).hexdigest() for name, content in contents.items()},
    })
    package_rel = PUBLICATION_ROOT_REL / fingerprint
    package = root / package_rel
    if not package.exists() and not package.is_symlink():
        stage = _new_stage(runtime, root, PUBLICATION_STAGE_REL, fingerprint)
        try:
            for name, content in contents.items():
                _write_file(runtime, stage / name, content)
            _, manifest_content = _package_manifest(runtime, contents, targets, fingerprint)
            _write_file(runtime, stage / "package.manifest.json", manifest_content)
            destination_parent = runtime.ensure_directory(root, package_rel.parent)
            publish_directory_no_replace(runtime, root, stage, destination_parent / fingerprint)
        except BaseException:
            raise
    package_verified = _verify_publication_package(
        runtime, root, package_rel, contents, targets, fingerprint,
    )
    # The provenance is a commit marker. Existing exact-prefix members from a
    # SIGKILL are verified byte-for-byte and the deterministic promotion
    # resumes; no conflicting artifact is replaced.
    promotion_order = [name for name in contents if name != "results.provenance.json"]
    promotion_order.append("results.provenance.json")
    promoted: dict[str, dict[str, object]] = {}
    for name in promotion_order:
        expected = package_verified["manifest"]["artifacts"][name]
        promoted[name] = _stream_publish_no_replace(
            runtime, root, package / name, targets[name], expected,
        )
    verified = validate_canonical_family(
        root, runtime=runtime, materialized=materialized,
        component_report_rel=component_report["paths"]["report"],
    )
    if verified["provenance"] != provenance:
        raise CoordinatorError("canonical verifier returned different provenance")
    return {"publication_fingerprint": fingerprint, "package": package_verified, **verified}


def validate_canonical_family(
    root: Path = ROOT, *, runtime: Any | None = None,
    materialized: Mapping[str, Any] | None = None,
    component_report_rel: Path | str | None = None,
) -> dict[str, Any]:
    """Side-effect-free exhaustive canonical family verifier."""

    runtime = runtime or _load_runtime(root)
    materialized = materialized or runtime.validate_materialized_manifest(root)
    policy = materialized.get("policy") or materialized.get("pre_family", {}).get("policy")
    if materialized["state"] == runtime.ZERO_STATE:
        for key in ("output_10", "output_11", "output_12"):
            if runtime.safe_path(root, Path(policy["required_science_outputs"][key]["path"]), "zero-family output").exists():
                raise CoordinatorError("zero family has a forbidden science placeholder")
        if runtime.safe_path(root, RESULT_PROVENANCE_REL, "zero-family result provenance").exists():
            raise CoordinatorError("zero family has nonzero canonical result provenance")
        return {"state": runtime.ZERO_STATE, "rows": [], "provenance": None}
    runs = _complete_runs(runtime, root, materialized)
    if component_report_rel is None:
        provenance_observed, _ = runtime.read_json(root, RESULT_PROVENANCE_REL)
        component_report_rel = provenance_observed.get("component_ram_report", {}).get("path")
    component = verify_component_ram_report(root, component_report_rel, runtime=runtime)
    component_shape = {
        "report_identity": component["report_identity"],
        "namespace_identity": component["namespace_identity"],
        "provenance_identity": component["provenance_identity"],
    }
    contents, expected_provenance = _canonical_expected(
        runtime, root, materialized, runs, component_shape,
    )
    targets = _canonical_targets(materialized["policy"])
    identities: dict[str, dict[str, object]] = {}
    for name, expected_content in contents.items():
        identity = runtime.stable_identity(root, targets[name])
        if identity["bytes"] != len(expected_content) or identity["sha256"] != hashlib.sha256(expected_content).hexdigest():
            raise CoordinatorError(f"canonical artifact differs from complete sealed family: {targets[name]}")
        identities[name] = identity
    observed_provenance, provenance_identity = runtime.read_json(root, RESULT_PROVENANCE_REL)
    if observed_provenance != expected_provenance:
        raise CoordinatorError("canonical result provenance drifted")
    return {
        "state": "COMPLETE_CANONICAL_FAMILY", "runs": runs,
        "identities": identities, "provenance": observed_provenance,
        "provenance_identity": provenance_identity,
        "component_report": component,
    }


def execute_family(
    root: Path = ROOT, *, configured_envelope_bytes: int = PHYSICAL_ENVELOPE_BYTES,
    larger_host_continuation: bool = False,
) -> dict[str, Any]:
    """Materialize, freeze, execute, retain, and publish the frozen family."""

    runtime = _load_runtime(root)
    with coordinator_lease(runtime, root):
        family = runtime.validate_pre_materialization_family(root)
        if family["state"] == runtime.ZERO_STATE:
            frozen = runtime.freeze_materialized_manifest(root)
            verified = runtime.validate_materialized_manifest(root)
            if frozen["state"] != runtime.ZERO_STATE or verified["state"] != runtime.ZERO_STATE:
                raise CoordinatorError("zero family failed terminal provenance verification")
            return {
                "state": runtime.ZERO_STATE, "mapped_loci": 0, "unavailable_loci": 0,
                "science_outputs_created": False, "zero_provenance": verified["zero_provenance"],
            }

        materialization_outcomes: dict[str, dict[str, Any]] = {}
        for pre_row in family["rows"]:
            outcome = benchmark_materialization(
                runtime, root, family, pre_row,
                configured_envelope_bytes=configured_envelope_bytes,
                larger_host_continuation=larger_host_continuation,
            )
            materialization_outcomes[pre_row["locus_entry_id"]] = outcome
        failed_materialization = {
            locus: outcome for locus, outcome in materialization_outcomes.items()
            if outcome["terminal_status"] != "COMPLETE_MATERIALIZED"
        }
        if failed_materialization:
            component = publish_component_ram_report(
                runtime, root, family, materialization_outcomes,
            )
            return {
                "state": "INCOMPLETE_MATERIALIZATION_RETAINED",
                "mapped_loci": len(family["rows"]),
                "unavailable_loci": len(family["unavailable_rows"]),
                "terminal_statuses": {
                    locus: outcome["terminal_status"] for locus, outcome in failed_materialization.items()
                },
                "component_report": component["report_identity"],
                "science_outputs_created": False,
            }

        # This lock is frozen only after every mapped input bundle and both
        # materialization measurements deeply revalidate. It precedes any R
        # engine launch/result access.
        materialized = runtime.freeze_materialized_manifest(root)
        if materialized["state"] != "READY_FRESH_PROCESS_SUSIE_COLOC":
            raise CoordinatorError("nonzero family did not freeze a runnable materialized manifest")
        materialized = runtime.validate_materialized_manifest(root)
        by_id = {row["locus_entry_id"]: row for row in materialized["rows"]}
        execution_order = runtime.representative_execution_order(materialized["rows"])
        if set(execution_order) != set(by_id) or len(execution_order) != len(by_id):
            raise CoordinatorError("representative-first execution order lost a locus")
        run_outcomes: dict[str, dict[str, Any]] = {}
        for locus in execution_order:
            run_outcomes[locus] = execute_locus(
                runtime, root, materialized, by_id[locus],
                configured_envelope_bytes=configured_envelope_bytes,
                larger_host_continuation=larger_host_continuation,
            )
        component = publish_component_ram_report(
            runtime, root, family, materialization_outcomes, run_outcomes,
        )
        incomplete = {
            locus: result["terminal_status"] for locus, result in run_outcomes.items()
            if result["terminal_status"] != "COMPLETE_SEALED"
        }
        if incomplete or (not materialized["rows"] and family["unavailable_rows"]):
            return {
                "state": (
                    "BLOCKED_BY_DATA_ONLY_FAMILY_RETAINED"
                    if not materialized["rows"] and family["unavailable_rows"]
                    else "INCOMPLETE_FAMILY_RETAINED"
                ),
                "mapped_loci": len(family["rows"]),
                "unavailable_loci": len(family["unavailable_rows"]),
                "terminal_statuses": incomplete,
                "blocked_by_data": [row["unavailable_entry_id"] for row in family["unavailable_rows"]],
                "component_report": component["report_identity"],
                "science_outputs_created": False,
            }
        ordered_runs = [run_outcomes[row["locus_entry_id"]] for row in materialized["rows"]]
        canonical = publish_canonical_family(
            runtime, root, materialized, ordered_runs, component,
        )
        return {
            "state": "COMPLETE_CANONICAL_FAMILY",
            "mapped_loci": len(materialized["rows"]),
            "unavailable_loci": len(family["unavailable_rows"]),
            "execution_order": execution_order,
            "component_report": component["report_identity"],
            "result_provenance": canonical["provenance_identity"],
            "science_outputs_created": True,
        }


def verify_current_state(root: Path = ROOT) -> dict[str, Any]:
    """Side-effect-free audit of the deepest currently published boundary."""

    runtime = _load_runtime(root)
    family = runtime.validate_pre_materialization_family(root)
    if family["state"] == runtime.ZERO_STATE:
        materialized = runtime.validate_materialized_manifest(root)
        canonical = validate_canonical_family(root, runtime=runtime, materialized=materialized)
        return {"state": runtime.ZERO_STATE, "family": family, "canonical": canonical}
    materialized_lock = runtime.safe_path(root, runtime.MATERIALIZED_LOCK_REL, "materialized lock")
    if not materialized_lock.exists():
        outcomes: dict[str, Any] = {}
        failed_attempts: dict[str, list[dict[str, Any]]] = {}
        for row in family["rows"]:
            failed_attempts[row["locus_entry_id"]] = _failed_materialization_attempts(
                runtime, root, family, row["locus_entry_id"],
            )
            measurement_path = root / MATERIALIZATION_METRICS_REL / row["locus_entry_id"]
            if measurement_path.exists():
                outcomes[row["locus_entry_id"]] = verify_materialization_measurement(
                    root, row["locus_entry_id"], family=family, runtime=runtime,
                )
        return {
            "state": "PRE_MATERIALIZATION_OR_INCOMPLETE_MATERIALIZATION",
            "verified_materializations": outcomes,
            "verified_failed_materialization_attempts": failed_attempts,
            "mapped_loci": len(family["rows"]),
            "unavailable_loci": len(family["unavailable_rows"]),
        }
    materialized = runtime.validate_materialized_manifest(root)
    attempts: dict[str, list[dict[str, Any]]] = {}
    for row in materialized["rows"]:
        attempts[row["locus_entry_id"]] = _existing_runs(
            runtime, root, materialized, row["locus_entry_id"],
        )
    if runtime.safe_path(root, RESULT_PROVENANCE_REL, "canonical result provenance").exists():
        canonical = validate_canonical_family(root, runtime=runtime, materialized=materialized)
        return {"state": "COMPLETE_CANONICAL_FAMILY", "canonical": canonical}
    return {
        "state": "MATERIALIZED_OR_PARTIAL_RUN_FAMILY", "materialized": materialized,
        "attempts": attempts,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--execute", action="store_true")
    actions.add_argument("--verify", action="store_true")
    actions.add_argument("--verify-component", metavar="REPORT_TSV")
    actions.add_argument("--print-interface", action="store_true")
    parser.add_argument("--memory-envelope-gib", type=float, default=8.0)
    parser.add_argument("--larger-host-continuation", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve(strict=True)
    if args.print_interface:
        runtime = _load_runtime(root)
        print(json.dumps({
            "analysis_id": ANALYSIS_ID,
            "successful_rss_backend": PGID_BACKEND,
            "execution_environment_requirement": (
                "UNSANDBOXED_PS_PROCESS_INVENTORY_PERMISSION_REQUIRED_FOR_PRODUCTION_EXECUTE"
            ),
            "component_ram_fields": COMPONENT_RAM_FIELDS,
            "component_collision_key": ["analysis", "pair", "locus", "chromosome"],
            "resource_metrics_fields": runtime.load_policy(root)[0]["ram_aware_execution_contract"]["resource_metrics_schema"],
            "public_verifiers": [
                "verify_materialization_measurement", "verify_failed_materialization_attempt",
                "verify_sealed_run",
                "verify_component_ram_report", "validate_canonical_family", "verify_current_state",
            ],
        }, indent=2, sort_keys=True))
        return 0
    if args.verify_component:
        result = verify_component_ram_report(root, args.verify_component)
        print(json.dumps({
            "status": "COMPONENT_RAM_REPORT_VERIFIED", "rows": len(result["rows"]),
            "report": result["report_identity"], "namespace": result["namespace_identity"],
            "provenance": result["provenance_identity"],
        }, sort_keys=True))
        return 0
    if args.verify:
        result = verify_current_state(root)
        print(json.dumps({"status": "TRACK_B_FINEMAPPING_STATE_VERIFIED", "state": result["state"]}, sort_keys=True))
        return 0
    envelope = int(args.memory_envelope_gib * 1024**3)
    if not math.isfinite(args.memory_envelope_gib) or envelope < PHYSICAL_ENVELOPE_BYTES:
        raise CoordinatorError("memory envelope must be finite and at least 8 GiB")
    result = execute_family(
        root, configured_envelope_bytes=envelope,
        larger_host_continuation=args.larger_host_continuation,
    )
    print(json.dumps(result, sort_keys=True, default=str))
    return 0 if result["state"] in {"COMPLETE_CANONICAL_FAMILY", _load_runtime(root).ZERO_STATE} else 3


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (CoordinatorError, OSError) as error:
        raise SystemExit(f"ERROR: {error}") from error
