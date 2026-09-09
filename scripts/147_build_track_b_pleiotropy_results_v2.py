#!/usr/bin/env python3
"""Build the real Track B post-PLACO result family.

This is an additive, fail-closed continuation of scripts 123--143.  It does
not estimate PLACO nuisance parameters and it never consumes LAVA rows as
scientific support.  The complete, already validated A/B/CONTROL PLACO
ledgers are streamed into table 07; eligible variants are exact-matched to the
sealed 1000 Genomes EUR BIM, clumped with pinned PLINK, independently checked
for pairwise r2 < 0.1 within 1 Mb, and only then mapped to the official LAVA
blocks as downstream work partitions.

The pre-result contract must be sealed before this module opens a PLACO result
ledger.  All production mutations require ``--execute``.  Publication is a
no-replace, crash-resumable hard-link commit whose result provenance is linked
last.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import errno
import fcntl
import gzip
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
import sqlite3
import stat
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Iterable, Iterator, Mapping, NamedTuple, Sequence, TextIO


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path("scripts/147_build_track_b_pleiotropy_results_v2.py")
POLICY = Path("config/track_b_pleiotropy_policy.json")
PLEIOTROPY_CONTRACT_SCRIPT = Path("scripts/123_track_b_pleiotropy_contract.py")
PLEIOTROPY_CONTRACT_LOCK = Path("results/track_b/pleiotropy/contract.lock.json")
INPUT_GATE_LOCK = Path("results/track_b/pleiotropy/input_gate.lock.json")
LD_MATERIALIZER_SCRIPT = Path("scripts/127_prepare_track_b_pleiotropy_ld.py")
TERMINAL_GATE_SCRIPT = Path("scripts/139_build_track_b_placo_terminal_gate_v2.py")
TERMINAL_GATE_LOCK = Path(
    "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/terminal_gate.lock.json"
)
SEQUENTIAL_SCRIPT = Path("scripts/143_run_track_b_placo_sequential_v2.py")
SEQUENTIAL_CONTRACT = Path(
    "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/"
    "sequential_execution.contract.json"
)
LAVA_EXPORTER_SCRIPT = Path("scripts/142_export_track_b_lava_terminal_diagnostics.py")
ARCHIVE_GUARD_PATH = "scripts/146_manage_lava_reference_archives.py"
ARCHIVE_GUARD_SCRIPT = Path(ARCHIVE_GUARD_PATH)
OFFICIAL_BLOCKS = Path("ref/lava/blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")
LOCAL_POLICY = Path("config/track_b_local_analysis_policy.json")
LOCAL_INPUT_LOCK = Path("results/track_b/local_analysis_input.lock.json")

CONTRACT_ROOT = Path(
    "results/track_b/pleiotropy/continuations/post_placo_collation_v2"
)
COLLATION_CONTRACT = CONTRACT_ROOT / "collation.contract.json"
EXECUTION_MUTEX = CONTRACT_ROOT / ".execution.lock"
PACKAGE_ROOT = Path("results/track_b/pleiotropy/results/post_placo_collation_v2")
PACKAGE_STAGING_ROOT = Path("results/track_b/pleiotropy/.post_placo_collation_staging")

PLAC0_07_PATH = Path("results/track_b/07_placo_plus_variants.tsv")
SHARED_LOCI_08_PATH = Path("results/track_b/08_shared_loci.tsv")
COMPARISON_09_PATH = Path("results/track_b/09_pleiotropy_comparison.tsv")
RESULT_PROVENANCE_PATH = Path("results/track_b/pleiotropy/results/results.provenance.json")
LD_COVERAGE_PATH = Path("results/track_b/pleiotropy/results/placo_ld_coverage.tsv")
ZERO_FAMILY_PATH = Path("results/track_b/pleiotropy/results/zero_family.provenance.json")
BLOCKED_COVERAGE_PATH = Path(
    "results/track_b/pleiotropy/results/ld_coverage_blocked_family.provenance.json"
)
CONTROL_ROOT = Path("results/track_b/control_insomnia_frailty")
CONTROL_08_PATH = CONTROL_ROOT / "08_shared_loci.tsv"
CONTROL_09_PATH = CONTROL_ROOT / "09_pleiotropy_comparison.tsv"
CONTROL_COVERAGE_PATH = CONTROL_ROOT / "ld_reference_coverage.tsv"
RAM_COMPONENT_PATH = Path(
    "results/track_b/pleiotropy/results/RAM_BENCHMARK.post_placo_collation_v2.tsv"
)
RAM_NAMESPACE_PATH = Path(
    "results/track_b/pleiotropy/results/RAM_BENCHMARK.post_placo_collation_v2.namespace.json"
)
RAM_REPORT_PATH = Path(
    "results/track_b/pleiotropy/results/RAM_AWARE_EXECUTION_REPORT.post_placo_collation_v2.md"
)

# Literal constants and an executable call are intentionally present for
# scripts/146_manage_lava_reference_archives.py --additional-consumer.
ARCHIVES_PRESENT = "ARCHIVES_PRESENT_FULLY_VERIFIED"
ARCHIVES_EVICTED = (
    "ARCHIVES_INTENTIONALLY_EVICTED_WITH_FINAL_RECEIPT_AND_ALL_44_"
    "EXTRACTED_PAYLOADS_REHASHED"
)
ACCEPTED_ARCHIVE_STATES = {ARCHIVES_PRESENT, ARCHIVES_EVICTED}

SCHEMA = "sleep-atlas-track-b-post-placo-collation.2"
CONTRACT_SCHEMA = "sleep-atlas-track-b-post-placo-collation-contract.2"
ZERO_SCHEMA = "sleep-atlas-track-b-post-placo-zero-family.2"
VERIFIER_SCHEMA = "sleep-atlas-track-b-post-placo-verified-family.2"
OFFICIAL_BLOCK_SHA256 = "462e81bce9ca85c9f11c0feb4d1bda24dd2c0fe095bfcc83fd32f77dff9a0882"
PAIR_ORDER = ("A", "B", "CONTROL")
PAIR_IDENTITIES = {
    "A": ("snoring", "parental_lifespan", "PRIMARY_DISCOVERY"),
    "B": ("insomnia", "adhd", "PRIMARY_DISCOVERY"),
    "CONTROL": ("insomnia", "frailty", "POSITIVE_CONTROL"),
}
COMPLETE_STATES = {"COMPLETE_WITH_HITS", "TESTED_NO_HIT"}
CONJFDR_BLOCKED = "NOT_TESTED_METHOD_BLOCKED"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
RSID_RE = re.compile(r"^rs[0-9]+$", re.IGNORECASE)
BASES = {"A", "C", "G", "T"}
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}

PLACO_LEDGER_FIELDS = [
    "analysis_id", "pair_id", "SNP", "CHR", "BP", "A1", "A2", "Z1", "Z2", "P1", "P2",
    "T_PLACO_PLUS", "P_PLACO_PLUS", "PLACO_BH_Q", "within_pair_family_n",
    "analysis_status", "numerical_error",
]
PLAC0_07_FIELDS = ["analysis_id", "pair_id", "family_role", *PLACO_LEDGER_FIELDS[2:]]
SHARED_LOCI_08_PREFIX = [
    "pair_id", "lead_variant", "chr", "position", "PLACO_P", "FDR",
    "trait1_P", "trait2_P", "locus_start", "locus_end",
    "independent_signal", "annotations",
]
SHARED_LOCI_08_EXTRA_FIELDS = [
    "analysis_id", "family_role", "trait1", "trait2", "lead_A1", "lead_A2",
    "lead_Z1", "lead_Z2", "evidence_labels", "primary_headline",
    "eligibility_basis", "clump_id", "clump_evidence_ids", "lava_block_id",
    "lava_block_start", "lava_block_end", "pair_block_evidence_ids",
    "ld_reference_status", "plink_clump_status", "pairwise_r2_status",
    "max_pairwise_r2_within_1mb", "lava_terminal_state",
    "lava_failed_qc_accounting", "conditional_status", "method_availability",
    "claim_limit",
]
SHARED_LOCI_08_FIELDS = SHARED_LOCI_08_PREFIX + SHARED_LOCI_08_EXTRA_FIELDS
COMPARISON_09_FIELDS = [
    "analysis_id", "pair_id", "family_role", "union_locus_id", "CHR", "START", "STOP",
    "placo_lead_snp", "placo_P", "placo_BH_Q", "conjfdr_lead_snp", "conjfdr",
    "evidence_labels", "method_availability", "comparison_label", "claim_limit",
]
LD_COVERAGE_FIELDS = [
    "analysis_id", "pair_id", "family_role", "evidence_id", "SNP", "CHR", "BP",
    "A1", "A2", "P_PLACO_PLUS", "PLACO_BH_Q", "eligibility_basis",
    "reference_status", "reference_CHR", "reference_BP", "reference_A1", "reference_A2",
    "clump_status", "failure_reason",
]
RAM_FIELDS = [
    "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
    "runtime_sec", "exit_status", "output_hash",
]

RESULT_PATHS = {
    "complete_placo_07": PLAC0_07_PATH,
    "primary_shared_loci_08": SHARED_LOCI_08_PATH,
    "primary_comparison_09": COMPARISON_09_PATH,
    "control_shared_loci_08": CONTROL_08_PATH,
    "control_comparison_09": CONTROL_09_PATH,
    "ld_coverage": LD_COVERAGE_PATH,
    "control_ld_coverage": CONTROL_COVERAGE_PATH,
    "result_provenance": RESULT_PROVENANCE_PATH,
}
RESULT_SCHEMAS = {
    "complete_placo_07": PLAC0_07_FIELDS,
    "shared_loci_08": SHARED_LOCI_08_FIELDS,
    "pleiotropy_comparison_09": COMPARISON_09_FIELDS,
    "ld_coverage": LD_COVERAGE_FIELDS,
    "ram_benchmark": RAM_FIELDS,
}

CANONICAL_SCIENCE_PATHS = {
    PLAC0_07_PATH, SHARED_LOCI_08_PATH, COMPARISON_09_PATH,
    RESULT_PROVENANCE_PATH, ZERO_FAMILY_PATH, BLOCKED_COVERAGE_PATH,
    CONTROL_08_PATH, CONTROL_09_PATH, LD_COVERAGE_PATH, CONTROL_COVERAGE_PATH,
    RAM_COMPONENT_PATH, RAM_NAMESPACE_PATH, RAM_REPORT_PATH,
}
COMPETING_PATHS = {
    Path("results/track_b/finemapping/locus_manifest.tsv"),
    Path("results/track_b/finemapping/locus_manifest.lock.json"),
    Path("results/track_b/finemapping/zero_family.provenance.json"),
}
PUBLICATION_ORDER_NONZERO = [
    PLAC0_07_PATH, LD_COVERAGE_PATH, SHARED_LOCI_08_PATH, COMPARISON_09_PATH,
    CONTROL_COVERAGE_PATH, CONTROL_08_PATH, CONTROL_09_PATH,
    RAM_COMPONENT_PATH, RAM_NAMESPACE_PATH, RAM_REPORT_PATH, RESULT_PROVENANCE_PATH,
]
PUBLICATION_ORDER_ZERO = [
    PLAC0_07_PATH, LD_COVERAGE_PATH, CONTROL_COVERAGE_PATH,
    RAM_COMPONENT_PATH, RAM_NAMESPACE_PATH, RAM_REPORT_PATH,
    ZERO_FAMILY_PATH, RESULT_PROVENANCE_PATH,
]
PUBLICATION_ORDER_BLOCKED = [
    PLAC0_07_PATH, LD_COVERAGE_PATH, CONTROL_COVERAGE_PATH,
    RAM_COMPONENT_PATH, RAM_NAMESPACE_PATH, RAM_REPORT_PATH,
    BLOCKED_COVERAGE_PATH, RESULT_PROVENANCE_PATH,
]


class CollationError(RuntimeError):
    """A fail-closed contract, scientific, or publication violation."""


class Block(NamedTuple):
    locus: int
    chromosome: int
    start: int
    stop: int


class PairInput(NamedTuple):
    pair_id: str
    ledger: Path
    provenance: Path
    terminal_state: str
    expected_rows: int
    ledger_identity: Mapping[str, Any]
    provenance_identity: Mapping[str, Any]


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise CollationError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _canonical_json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _digest_json(value: object) -> str:
    return hashlib.sha256(_canonical_json(value)).hexdigest()


def _safe_path(root: Path, relative: Path | str, label: str) -> Path:
    relative = Path(relative)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise CollationError(f"unsafe {label} path: {relative}")
    try:
        root_real = root.resolve(strict=True)
    except OSError as error:
        raise CollationError(f"repository root is unavailable: {error}") from error
    current = root
    for offset, part in enumerate(relative.parts):
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise CollationError(f"could not inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise CollationError(f"{label} contains a symbolic link: {current}")
        if offset < len(relative.parts) - 1 and not stat.S_ISDIR(observed.st_mode):
            raise CollationError(f"{label} has a non-directory ancestor: {current}")
    try:
        (root / relative).resolve(strict=False).relative_to(root_real)
    except (OSError, RuntimeError, ValueError) as error:
        raise CollationError(f"{label} escapes repository: {relative}") from error
    return root / relative


def _relative(root: Path, path: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError as error:
        raise CollationError(f"path escapes repository: {path}") from error


def _stat_tuple(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def stable_identity(
    root: Path, relative: Path | str, *, allow_empty: bool = False,
) -> dict[str, object]:
    path = _safe_path(root, relative, "artifact")
    digest = hashlib.sha256()
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        named_before = os.lstat(path)
        descriptor = os.open(path, flags)
        try:
            opened = os.fstat(descriptor)
            if (
                not stat.S_ISREG(named_before.st_mode)
                or not stat.S_ISREG(opened.st_mode)
                or (opened.st_size == 0 and not allow_empty)
                or (named_before.st_dev, named_before.st_ino)
                != (opened.st_dev, opened.st_ino)
            ):
                raise CollationError(f"artifact is not one stable regular file: {path}")
            while True:
                chunk = os.read(descriptor, 4 * 1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        named_after = os.lstat(path)
    except CollationError:
        raise
    except OSError as error:
        raise CollationError(f"could not hash artifact {path}: {error}") from error
    if (
        len({_stat_tuple(item) for item in (named_before, opened, after, named_after)}) != 1
        or not stat.S_ISREG(named_after.st_mode)
    ):
        raise CollationError(f"artifact changed while hashing: {path}")
    return {
        "path": _relative(root, path), "device": int(after.st_dev), "inode": int(after.st_ino),
        "bytes": int(after.st_size), "mtime_ns": int(after.st_mtime_ns),
        "ctime_ns": int(after.st_ctime_ns), "sha256": digest.hexdigest(),
    }


def portable_identity(root: Path, relative: Path | str, *, allow_empty: bool = False) -> dict[str, object]:
    observed = stable_identity(root, relative, allow_empty=allow_empty)
    return {key: observed[key] for key in ("path", "bytes", "sha256")}


def identity_matches(root: Path, record: Mapping[str, Any]) -> bool:
    try:
        observed = portable_identity(root, str(record.get("path", "")), allow_empty=record.get("bytes") == 0)
    except CollationError:
        return False
    return all(observed.get(key) == record.get(key) for key in ("path", "bytes", "sha256"))


def read_json(root: Path, relative: Path | str) -> dict[str, Any]:
    path = _safe_path(root, relative, "JSON artifact")
    identity = stable_identity(root, relative)
    if int(identity["bytes"]) > 64 * 1024**2:
        raise CollationError(f"structured artifact exceeds 64-MiB bound: {path}")
    try:
        with stable_plain_text(
            root, relative,
            {key: identity[key] for key in ("path", "bytes", "sha256")},
        ) as handle:
            value = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CollationError(f"unreadable JSON artifact {path}: {error}") from error
    if not isinstance(value, dict):
        raise CollationError(f"JSON artifact must contain one object: {path}")
    return value


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _ensure_directories(root: Path, relative: Path) -> Path:
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise CollationError(f"unsafe directory path: {relative}")
    current = root
    for part in relative.parts:
        candidate = current / part
        try:
            os.mkdir(candidate)
        except FileExistsError:
            observed = os.lstat(candidate)
            if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
                raise CollationError(f"directory component is not a real directory: {candidate}")
        else:
            _fsync_directory(current)
        current = candidate
    return _safe_path(root, relative, "directory")


def publish_bytes_no_replace(root: Path, relative: Path, content: bytes) -> dict[str, object]:
    parent = _ensure_directories(root, relative.parent)
    path = _safe_path(root, relative, "publication")
    if path.exists() or path.is_symlink():
        raise CollationError(f"refusing to replace existing artifact: {relative}")
    digest = hashlib.sha256(content).hexdigest()
    temporary = parent / f".{path.name}.{digest}.{os.getpid()}.tmp"
    if temporary.exists() or temporary.is_symlink():
        raise CollationError(f"publication temporary already exists: {temporary}")
    linked: tuple[int, int] | None = None
    try:
        with temporary.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        source = os.lstat(temporary)
        if not stat.S_ISREG(source.st_mode):
            raise CollationError("publication temporary is not regular")
        try:
            os.link(temporary, path, follow_symlinks=False)
        except FileExistsError as error:
            raise CollationError(f"publication destination appeared concurrently: {relative}") from error
        linked = source.st_dev, source.st_ino
        destination = os.lstat(path)
        if not stat.S_ISREG(destination.st_mode) or linked != (destination.st_dev, destination.st_ino):
            raise CollationError("no-replace publication did not link the exact source inode")
        _fsync_directory(parent)
    except BaseException:
        if linked is not None:
            try:
                current = os.lstat(path)
                if (current.st_dev, current.st_ino) == linked:
                    path.unlink()
                    _fsync_directory(parent)
            except FileNotFoundError:
                pass
        raise
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return portable_identity(root, relative, allow_empty=len(content) == 0)


@contextmanager
def execution_lock(root: Path) -> Iterator[None]:
    parent = _ensure_directories(root, EXECUTION_MUTEX.parent)
    path = parent / EXECUTION_MUTEX.name
    flags = os.O_RDWR | os.O_CREAT
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
            or (opened.st_dev, opened.st_ino) != (named.st_dev, named.st_ino)
            or opened.st_nlink != 1 or named.st_nlink != 1
        ):
            raise CollationError("collation mutex is not one private regular file")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise CollationError("another post-PLACO collation action holds the mutex") from error
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def _artifact_present(path: Path) -> bool:
    return path.exists() or path.is_symlink()


def assert_no_competing_publication(root: Path) -> None:
    offenders = [
        relative for relative in sorted(CANONICAL_SCIENCE_PATHS | COMPETING_PATHS, key=str)
        if _artifact_present(root / relative)
    ]
    package = root / PACKAGE_ROOT
    if _artifact_present(package):
        if package.is_symlink() or not package.is_dir() or any(package.iterdir()):
            offenders.append(PACKAGE_ROOT)
    control = root / CONTROL_ROOT
    if _artifact_present(control):
        if control.is_symlink() or not control.is_dir() or any(control.iterdir()):
            offenders.append(CONTROL_ROOT)
    if offenders:
        raise CollationError(
            "pre-result collation contract conflicts with canonical/competing evidence: "
            + ", ".join(str(path) for path in offenders[:8])
        )


def assert_no_competing_finemapping(root: Path) -> None:
    offenders = [
        relative for relative in sorted(COMPETING_PATHS, key=str)
        if _artifact_present(root / relative)
    ]
    if offenders:
        raise CollationError(
            "competing fine-mapping evidence exists before post-PLACO execution: "
            + ", ".join(str(path) for path in offenders)
        )


def verify_archive_state(guard: Any, root: Path) -> dict[str, Any]:
    state = guard.verify_reference_state(root)
    if not isinstance(state, dict) or state.get("state") not in ACCEPTED_ARCHIVE_STATES:
        raise CollationError("LAVA archive guard returned a state outside its two-state interface")
    if (
        state.get("extracted_payload_count") != 44
        or not SHA256_RE.fullmatch(str(state.get("archive_family_sha256", "")))
        or not SHA256_RE.fullmatch(str(state.get("extracted_family_sha256", "")))
        or not isinstance(state.get("reference_provenance_identity"), Mapping)
    ):
        raise CollationError("LAVA archive guard returned incomplete reference identity")
    return state


def archive_state_invariant(root: Path, state: Mapping[str, Any]) -> dict[str, Any]:
    """Freeze scientific identity while allowing only the guard's two lifecycle states."""
    return {
        "accepted_states": [ARCHIVES_PRESENT, ARCHIVES_EVICTED],
        "guard_implementation": portable_identity(root, ARCHIVE_GUARD_SCRIPT),
        "reference_provenance": state["reference_provenance_identity"],
        "archive_family_sha256": state["archive_family_sha256"],
        "extracted_payload_count": state["extracted_payload_count"],
        "extracted_family_sha256": state["extracted_family_sha256"],
        "verification": (
            "CURRENT_STATE_REVALIDATED_EVERY_CALL;MUTABLE_PRESENT_OR_EVICTED_TOKEN_"
            "INTENTIONALLY_EXCLUDED_FROM_IMMUTABLE_COLLATION_CONTRACT"
        ),
    }


def _validate_official_blocks(root: Path) -> dict[str, Any]:
    identity = portable_identity(root, OFFICIAL_BLOCKS)
    local_policy = read_json(root, LOCAL_POLICY)
    local_lock = read_json(root, LOCAL_INPUT_LOCK)
    locus = local_policy.get("locus_definition", {})
    if (
        identity["sha256"] != OFFICIAL_BLOCK_SHA256
        or locus.get("path") != str(OFFICIAL_BLOCKS)
        or locus.get("sha256") != OFFICIAL_BLOCK_SHA256
        or local_policy.get("expected_loci") != 2495
        or local_lock.get("locus_file_sha256") != OFFICIAL_BLOCK_SHA256
    ):
        raise CollationError("official 2,495-block LAVA partition identity drifted")
    load_official_blocks(root, identity)
    return identity | {"block_count": 2495, "role": "DOWNSTREAM_WORK_PARTITION_ONLY"}


def _load_production_modules(root: Path) -> dict[str, Any]:
    return {
        "contract": _load_module("track_b_pleiotropy_contract_for_147", root / PLEIOTROPY_CONTRACT_SCRIPT),
        "ld": _load_module("track_b_pleiotropy_ld_for_147", root / LD_MATERIALIZER_SCRIPT),
        "terminal": _load_module("track_b_placo_terminal_for_147", root / TERMINAL_GATE_SCRIPT),
        "sequential": _load_module("track_b_placo_sequential_for_147", root / SEQUENTIAL_SCRIPT),
        "lava_exporter": _load_module("track_b_lava_exporter_for_147", root / LAVA_EXPORTER_SCRIPT),
        "archive_guard": _load_module("track_b_lava_archive_guard_for_147", root / ARCHIVE_GUARD_SCRIPT),
    }


def validate_lava_accounting_live(
    root: Path, terminal_gate: Mapping[str, Any], exporter: Any,
) -> dict[str, Any]:
    """Deeply verify the upstream 142 failed-QC diagnostic family directly.

    This deliberately has no dependency on script 145: fine-mapping is a
    downstream consumer of this module, never an input to its freeze.
    """
    evidence = terminal_gate.get("lava_scientific_evidence", {})
    terminal = terminal_gate.get("lava_terminal", {})
    payload = terminal.get("terminal_payload", {})
    if (
        terminal_gate.get("schema_version") != "sleep-atlas-track-b-pleiotropy-terminal-gate.2"
        or evidence.get("status") != "FAILED_QC_NOT_CONSUMED"
        or evidence.get("scientific_validation_passed") is not False
        or evidence.get("passed_lava_results_used_by_placo") is not False
        or evidence.get("consumption") != "ORDERING_TERMINALITY_ATTESTATION_ONLY"
        or payload.get("state") != "TERMINAL_FAILED_QC"
        or payload.get("full_family_complete") is not True
        or payload.get("scientific_validation_passed") is not False
        or payload.get("canonical_publication_allowed") is not False
    ):
        raise CollationError("LAVA terminal gate is not complete FAILED_QC accounting-only evidence")
    source = terminal.get("source_discovery_fingerprint")
    continuation = terminal.get("continuation_execution_fingerprint")
    if not isinstance(source, str) or not SHA256_RE.fullmatch(source):
        raise CollationError("invalid LAVA source fingerprint")
    if not isinstance(continuation, str) or not SHA256_RE.fullmatch(continuation):
        raise CollationError("invalid LAVA continuation fingerprint")
    try:
        sealed = exporter.deep_validate_terminal(source, continuation, progress=None)
        family = exporter.validate_family(sealed)
        outputs, derived = exporter.render_outputs(family, sealed)
        export_path = exporter.verify_outputs(sealed, outputs, derived)
    except BaseException as error:
        raise CollationError(f"versioned 142 LAVA diagnostic family is invalid: {error}") from error
    export_relative = export_path.relative_to(root)
    export_lock = read_json(root, export_relative)
    output_records = export_lock.get("outputs")
    expected_paths = {
        "results/track_b/lava/local_univariate_results.tsv",
        "results/track_b/lava/local_bivariate_results.tsv",
        "results/track_b/lava/local_bivariate_fdr.tsv",
        "results/track_b/lava/local_failures.tsv",
        "results/track_b/lava/local_family_accounting.tsv",
        "results/track_b/LAVA_DISCOVERY_SUMMARY.md",
    }
    if (
        export_lock.get("schema_version") != "track-b-lava-terminal-diagnostic-export.1"
        or export_lock.get("scientific_state") != "TERMINAL_FAILED_QC"
        or export_lock.get("scientific_validation_passed") is not False
        or export_lock.get("canonical_publication_allowed") is not False
        or export_lock.get("source_discovery_fingerprint") != source
        or export_lock.get("continuation_execution_fingerprint") != continuation
        or export_lock.get("conditional_checkpoint_family", {}).get("count") != 2495
        or export_lock.get("derived_accounting", {}).get("official_loci") != 2495
        or not isinstance(output_records, list)
        or {record.get("path") for record in output_records if isinstance(record, Mapping)}
        != expected_paths
    ):
        raise CollationError("142 LAVA diagnostic export lock lost complete failed-QC accounting")
    for record in output_records:
        if not isinstance(record, Mapping) or not identity_matches(root, record):
            raise CollationError("142 LAVA diagnostic output identity drifted")
    return {
        "scientific_state": "TERMINAL_FAILED_QC",
        "consumption": "ACCOUNTING_ONLY_NO_SELECTION_SUPPORT_OR_RANKING",
        "source_discovery_fingerprint": source,
        "continuation_execution_fingerprint": continuation,
        "terminal_receipt": terminal.get("receipt"),
        "terminal_gate": portable_identity(root, TERMINAL_GATE_LOCK),
        "diagnostic_export_lock": portable_identity(root, export_relative),
        "diagnostic_outputs": output_records,
        "derived_accounting": export_lock.get("derived_accounting"),
        "all_lava_loci_excluded_from_scientific_support": True,
    }


def build_pre_result_contract(root: Path = ROOT) -> dict[str, Any]:
    """Build the additive contract without opening any PLACO result ledger."""
    modules = _load_production_modules(root)
    try:
        policy = modules["contract"].validate_policy(root)
        original_contract = modules["contract"].verify_contract(root)
        terminal = modules["terminal"].verify_gate()
        sequential = modules["sequential"].verify_execution_contract()
        policy_obj, ld_policy, policy_path, contract_sha = modules["ld"].load_policy_context(
            root, POLICY, test_fixture=False,
        )
        ld_provenance = modules["ld"].verify_reference(
            root, policy_obj, ld_policy, policy_path, contract_sha,
            test_fixture=False, quiet=True,
        )
        archive = verify_archive_state(modules["archive_guard"], root)
        lava = validate_lava_accounting_live(root, terminal, modules["lava_exporter"])
    except CollationError:
        raise
    except BaseException as error:
        raise CollationError(f"pre-result dependency validation failed: {error}") from error
    if policy_obj != policy or original_contract.get("analysis_id") != policy["analysis_id"]:
        raise CollationError("pleiotropy policy/contract loaders disagree")
    future_paths = policy.get("future_result_paths", {})
    future_schemas = policy.get("required_future_result_schemas", {})
    if (
        future_paths.get("placo_plus_variants_07") != str(PLACO_07_PATH)
        or future_paths.get("shared_loci_08") != str(SHARED_LOCI_08_PATH)
        or future_paths.get("pleiotropy_comparison_09") != str(COMPARISON_09_PATH)
        or future_paths.get("result_provenance") != str(RESULT_PROVENANCE_PATH)
        or future_schemas.get("placo_full_p_ledger") != PLACO_LEDGER_FIELDS
        or future_schemas.get("shared_loci_08_required_prefix") != SHARED_LOCI_08_PREFIX
        or policy.get("pair_order") != list(PAIR_ORDER)
        or policy.get("primary_pair_family") != ["A", "B"]
        or policy.get("control_pair_family") != ["CONTROL"]
    ):
        raise CollationError("collation paths/schemas/pair roles differ from frozen policy")
    plink = ld_provenance.get("plink", {})
    if (
        plink.get("path") != ld_policy["plink_path"]
        or plink.get("sha256") != ld_policy["plink_sha256"]
        or not str(plink.get("version", "")).startswith(ld_policy["plink_version"])
    ):
        raise CollationError("verified LD reference did not bind the pinned PLINK identity")
    block = _validate_official_blocks(root)
    identities = {
        str(path): portable_identity(root, path)
        for path in (
            POLICY, PLEIOTROPY_CONTRACT_LOCK, INPUT_GATE_LOCK, TERMINAL_GATE_LOCK,
            SEQUENTIAL_CONTRACT, PLEIOTROPY_CONTRACT_SCRIPT, LD_MATERIALIZER_SCRIPT,
            TERMINAL_GATE_SCRIPT, SEQUENTIAL_SCRIPT, LAVA_EXPORTER_SCRIPT,
            ARCHIVE_GUARD_SCRIPT, SCRIPT,
        )
    }
    return {
        "schema_version": CONTRACT_SCHEMA,
        "analysis_id": "track-b-v1.0-pleiotropy",
        "selection_timing": "BEFORE_POST_PLACO_RESULT_ACCESS",
        "placo_result_ledgers_opened_while_building_contract": False,
        "pair_order": list(PAIR_ORDER),
        "primary_pair_family": ["A", "B"],
        "control_pair_family": ["CONTROL"],
        "thresholds": {
            "eligible_union": "P_PLACO_PLUS<=5e-8 OR COMPLETE_WITHIN_PAIR_BH_Q<=0.05",
            "primary_headline_A_B_only": 2.5e-8,
            "clump_window_kb": 1000,
            "clump_r2": 0.1,
            "lead_independence": "STRICT_R2_LT_0.1_WITHIN_1_MB",
        },
        "conjfdr_state_for_this_implementation": CONJFDR_BLOCKED,
        "conjfdr_valid_attestation_policy": (
            "FAIL_CLOSED_IF_ATTESTED_FAMILY_EXISTS_UNTIL_SIGNAL_LEVEL_R2_GE_0.6_"
            "CROSS_METHOD_RECONCILIATION_IS_SEPARATELY_IMPLEMENTED"
        ),
        "upstream": identities,
        "terminal_gate_snapshot_sha256": _digest_json(terminal),
        "sequential_contract_snapshot_sha256": _digest_json(sequential),
        "ld_reference": {
            "provenance": portable_identity(root, ld_policy["materialized_provenance"]),
            "manifest": portable_identity(root, ld_policy["materialized_manifest"]),
            "bed": portable_identity(root, f"{ld_policy['materialized_prefix']}.bed"),
            "bim": portable_identity(root, f"{ld_policy['materialized_prefix']}.bim"),
            "fam": portable_identity(root, f"{ld_policy['materialized_prefix']}.fam"),
            "source_release": ld_policy["source_release"],
            "ancestry": ld_policy["ancestry"], "build": ld_policy["build"],
            "variant_count": ld_policy["variant_count"],
            "autosomal_variant_count": ld_policy["expected_autosomal_variant_count"],
            "chromosome_lengths": ld_policy["grch37_chromosome_lengths"],
            "plink": plink,
        },
        "official_lava_partition": block,
        "lava_failed_qc_accounting": lava,
        "archive_reference_state_guard": archive_state_invariant(root, archive),
        "output_schemas": {
            str(PLAC0_07_PATH): PLAC0_07_FIELDS,
            str(SHARED_LOCI_08_PATH): SHARED_LOCI_08_FIELDS,
            str(COMPARISON_09_PATH): COMPARISON_09_FIELDS,
            str(LD_COVERAGE_PATH): LD_COVERAGE_FIELDS,
            str(RAM_COMPONENT_PATH): RAM_FIELDS,
        },
        "publication": {
            "production_mutations_require_execute": True,
            "package_then_hardlink_no_replace": True,
            "result_provenance_is_last_commit_marker": True,
            "partial_exact_prefix_is_restartable": True,
            "canonical_overwrite": "FORBIDDEN",
        },
    }


def seal_pre_result_contract(root: Path = ROOT) -> dict[str, Any]:
    with execution_lock(root):
        assert_no_competing_publication(root)
        payload = build_pre_result_contract(root)
        path = root / COLLATION_CONTRACT
        if _artifact_present(path):
            if read_json(root, COLLATION_CONTRACT) != payload:
                raise CollationError("pre-result collation contract already exists with different bytes")
            return payload
        assert_no_competing_publication(root)
        publish_bytes_no_replace(
            root, COLLATION_CONTRACT,
            json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n",
        )
        if read_json(root, COLLATION_CONTRACT) != payload:
            raise CollationError("sealed pre-result collation contract failed verification")
        return payload


def verify_pre_result_contract(root: Path = ROOT) -> dict[str, Any]:
    observed = read_json(root, COLLATION_CONTRACT)
    expected = build_pre_result_contract(root)
    if observed != expected:
        raise CollationError("pre-result collation contract differs from live frozen dependencies")
    return observed


def materialize_ld_reference(root: Path = ROOT) -> dict[str, Any]:
    assert_no_competing_publication(root)
    module = _load_module("track_b_pleiotropy_ld_materialize_for_147", root / LD_MATERIALIZER_SCRIPT)
    policy, ld, policy_path, contract_hash = module.load_policy_context(root, POLICY, test_fixture=False)
    return module.materialize_reference(
        root, policy, ld, policy_path, contract_hash, test_fixture=False,
    )


# Public aliases used by downstream consumers; retain the misspelled internal
# name only for backward compatibility with early review copies of this file.
PLACO_07_PATH = PLAC0_07_PATH
PLACO_07_FIELDS = PLAC0_07_FIELDS


def parse_blocks(
    text: str, *, expected_count: int | None = None,
) -> tuple[list[Block], dict[int, tuple[list[int], list[Block]]]]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or lines[0].split() != ["LOC", "CHR", "START", "STOP"]:
        raise CollationError("official block file schema drifted")
    blocks: list[Block] = []
    for row_number, line in enumerate(lines[1:], start=2):
        fields = line.split()
        if len(fields) != 4:
            raise CollationError(f"malformed official block row {row_number}")
        try:
            locus, chromosome, start, stop = map(int, fields)
        except ValueError as error:
            raise CollationError(f"nonnumeric official block row {row_number}") from error
        if locus < 1 or chromosome not in range(1, 23) or start < 1 or stop < start:
            raise CollationError(f"invalid official block row {row_number}")
        blocks.append(Block(locus, chromosome, start, stop))
    if expected_count is not None and len(blocks) != expected_count:
        raise CollationError(
            f"official block family has {len(blocks)} rows, expected {expected_count}"
        )
    if [block.locus for block in blocks] != list(range(1, len(blocks) + 1)):
        raise CollationError("official block IDs are not the ordered 1..N family")
    if blocks != sorted(blocks, key=lambda item: (item.chromosome, item.start, item.stop)):
        raise CollationError("official block family is not in genomic order")
    grouped: dict[int, list[Block]] = defaultdict(list)
    for block in blocks:
        grouped[block.chromosome].append(block)
    if expected_count == 2495 and set(grouped) != set(range(1, 23)):
        raise CollationError("official block family does not cover all autosomes")
    index: dict[int, tuple[list[int], list[Block]]] = {}
    for chromosome, family in grouped.items():
        previous_stop = 0
        for block in family:
            if block.start <= previous_stop:
                raise CollationError(f"official blocks overlap on chromosome {chromosome}")
            previous_stop = block.stop
        index[chromosome] = ([block.start for block in family], family)
    return blocks, index


def load_official_blocks(
    root: Path, expected_identity: Mapping[str, Any],
) -> tuple[list[Block], dict[int, tuple[list[int], list[Block]]]]:
    with stable_plain_text(root, OFFICIAL_BLOCKS, expected_identity) as handle:
        return parse_blocks(handle.read(), expected_count=2495)


def map_to_block(
    chromosome: int, position: int,
    index: Mapping[int, tuple[list[int], list[Block]]],
) -> Block | None:
    family = index.get(chromosome)
    if family is None:
        return None
    starts, blocks = family
    offset = bisect.bisect_right(starts, position) - 1
    if offset < 0 or position > blocks[offset].stop:
        return None
    return blocks[offset]


def _number(value: Any, label: str) -> float:
    try:
        observed = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise CollationError(f"invalid numeric {label}: {value!r}") from error
    if not math.isfinite(observed):
        raise CollationError(f"non-finite numeric {label}: {value!r}")
    return observed


def _probability(value: Any, label: str) -> float:
    observed = _number(value, label)
    if not 0 <= observed <= 1:
        raise CollationError(f"{label} outside [0,1]: {value!r}")
    return observed


def evidence_labels(pair_id: str, p_value: float, q_value: float) -> list[str]:
    """Return every applicable policy label in frozen order."""
    if pair_id not in PAIR_ORDER:
        raise CollationError(f"unknown Track B pair: {pair_id}")
    labels: list[str] = []
    if pair_id == "CONTROL":
        if p_value <= 5e-8:
            labels.append("CONTROL_PLACO_GWS_RECOVERED")
        if q_value <= 0.05:
            labels.append("CONTROL_PLACO_WITHIN_PAIR_FDR")
    else:
        if p_value <= 2.5e-8:
            labels.append("PLACO_PRIMARY_HEADLINE")
        if p_value <= 5e-8:
            labels.append("PLACO_PAIRWISE_GWS")
        if q_value <= 0.05:
            labels.append("PLACO_WITHIN_PAIR_FDR")
    return labels


def eligibility_basis(pair_id: str, p_value: float, q_value: float) -> str | None:
    """The eligible union is P<=5e-8 OR complete-family BH q<=0.05."""
    del pair_id  # CONTROL has the same union threshold but a different role/label family.
    values: list[str] = []
    if p_value <= 5e-8:
        values.append("PLACO_P_LE_5E_8")
    if q_value <= 0.05:
        values.append("PLACO_COMPLETE_WITHIN_PAIR_BH_Q_LE_0_05")
    return ";".join(values) if values else None


def _encode_tsv_row(fields: Sequence[str], row: Mapping[str, Any]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
        extrasaction="raise",
    )
    writer.writerow({field: row[field] for field in fields})
    return buffer.getvalue().encode("utf-8")


def _tsv_header(fields: Sequence[str]) -> bytes:
    return ("\t".join(fields) + "\n").encode("utf-8")


@contextmanager
def stable_gzip_text(
    root: Path, relative: Path | str, expected: Mapping[str, Any] | None = None,
) -> Iterator[TextIO]:
    """Hash and parse one gzip through a stable open inode.

    Hashing happens before decompression on the same descriptor, then path,
    inode, size, mtime, and ctime are rechecked after the complete gzip read.
    This catches same-size mutation and named-inode replacement attacks.
    """
    path = _safe_path(root, relative, "compressed PLACO ledger")
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    text_handle: TextIO | None = None
    binary_handle: Any | None = None
    try:
        named_before = os.lstat(path)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(named_before.st_mode) or not stat.S_ISREG(opened.st_mode)
            or opened.st_size <= 0
            or (named_before.st_dev, named_before.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise CollationError(f"PLACO ledger is not one nonempty regular file: {path}")
        digest = hashlib.sha256()
        while True:
            chunk = os.read(descriptor, 4 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
        record = {
            "path": _relative(root, path), "bytes": int(opened.st_size),
            "sha256": digest.hexdigest(),
        }
        if expected is not None and any(record.get(key) != expected.get(key) for key in record):
            raise CollationError(f"PLACO ledger differs from its sealed identity: {path}")
        os.lseek(descriptor, 0, os.SEEK_SET)
        binary_handle = os.fdopen(os.dup(descriptor), "rb")
        gzip_handle = gzip.GzipFile(fileobj=binary_handle, mode="rb")
        text_handle = io.TextIOWrapper(gzip_handle, encoding="utf-8", newline="")
        yield text_handle
        # Force CRC/trailer validation even if a caller stopped at a boundary.
        while text_handle.read(1024 * 1024):
            pass
        text_handle.close()
        text_handle = None
        binary_handle.close()
        binary_handle = None
        after = os.fstat(descriptor)
        named_after = os.lstat(path)
        if (
            len({_stat_tuple(item) for item in (named_before, opened, after, named_after)}) != 1
            or not stat.S_ISREG(named_after.st_mode)
        ):
            raise CollationError(f"PLACO ledger changed while streaming: {path}")
    except (OSError, EOFError, UnicodeError, gzip.BadGzipFile) as error:
        raise CollationError(f"could not completely stream PLACO ledger {path}: {error}") from error
    finally:
        if text_handle is not None:
            text_handle.close()
        if binary_handle is not None:
            binary_handle.close()
        os.close(descriptor)


def _validate_pair_provenance(pair: PairInput, payload: Mapping[str, Any]) -> None:
    if (
        payload.get("schema_version") != "sleep-atlas-track-b-placo-result.1"
        or payload.get("analysis_id") != "track-b-v1.0-pleiotropy"
        or payload.get("pair_id") != pair.pair_id
        or payload.get("terminal_result_state") != pair.terminal_state
        or payload.get("terminal_result_state") not in COMPLETE_STATES
        or payload.get("qc_status") != "PASS"
        or payload.get("output_rows") != pair.expected_rows
        or payload.get("within_pair_bh_family_n") != pair.expected_rows
        or payload.get("output_sha256") != pair.ledger_identity.get("sha256")
        or payload.get("output_bytes") != pair.ledger_identity.get("bytes")
        or payload.get("exact_schema") != ",".join(PLACO_LEDGER_FIELDS)
        or payload.get("reference_sha256")
        != "NOT_APPLICABLE_TO_PLACO_FULL_P_SCAN_LD_REQUIRED_FOR_LATER_CLUMPING"
    ):
        raise CollationError(f"PLACO {pair.pair_id} provenance is not a complete PASS family")
    counts = payload.get("complete_family_counts")
    nuisance = payload.get("nuisance_estimation")
    if (
        not isinstance(counts, Mapping) or counts.get("rows") != pair.expected_rows
        or not isinstance(nuisance, Mapping) or nuisance.get("input_rows") != pair.expected_rows
        or nuisance.get("scope")
        != "ALL_VALID_ALIGNED_GENOME_WIDE_VARIANTS_BEFORE_ANY_SINGLE_VARIANT_SHARD"
    ):
        raise CollationError(f"PLACO {pair.pair_id} denominator or global nuisance scope drifted")


def resolve_production_pair_inputs(root: Path, sequential: Any) -> list[PairInput]:
    """Invoke the complete 143 family verifier, then resolve canonical inputs."""
    try:
        outcomes = sequential.verify_family()
    except BaseException as error:
        raise CollationError(f"complete deeply verified 143 PLACO family is unavailable: {error}") from error
    wanted = [{"pair": pair, "outcome": "VERIFIED_COMPLETE_CLEAN"} for pair in PAIR_ORDER]
    if outcomes != wanted:
        raise CollationError("143 verifier returned an unexpected A/B/CONTROL family")
    pairs: list[PairInput] = []
    for pair_id in PAIR_ORDER:
        try:
            context, _ = sequential.BRIDGE.configure_context(pair_id)
            paths = sequential.ENGINE.materialized_paths(context)
        except BaseException as error:
            raise CollationError(f"could not resolve verified PLACO pair {pair_id}: {error}") from error
        ledger_identity = portable_identity(root, paths["canonical_ledger"].relative_to(root))
        provenance_identity = portable_identity(root, paths["canonical_provenance"].relative_to(root))
        provenance = read_json(root, provenance_identity["path"])
        pair = PairInput(
            pair_id=pair_id,
            ledger=root / str(ledger_identity["path"]),
            provenance=root / str(provenance_identity["path"]),
            terminal_state=str(provenance.get("terminal_result_state", "")),
            expected_rows=int(provenance.get("output_rows", -1)),
            ledger_identity=ledger_identity,
            provenance_identity=provenance_identity,
        )
        _validate_pair_provenance(pair, provenance)
        pairs.append(pair)
    return pairs


def create_index_database(path: Path) -> sqlite3.Connection:
    if path.exists() or path.is_symlink():
        raise CollationError(f"collation index path already exists: {path}")
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA journal_mode=DELETE")
    connection.execute("PRAGMA synchronous=FULL")
    connection.executescript(
        """
        CREATE TABLE eligible (
          row_id INTEGER PRIMARY KEY,
          pair_id TEXT NOT NULL,
          family_role TEXT NOT NULL,
          source_row_index INTEGER NOT NULL,
          evidence_id TEXT NOT NULL,
          snp TEXT NOT NULL,
          chr INTEGER NOT NULL,
          bp INTEGER NOT NULL,
          a1 TEXT NOT NULL,
          a2 TEXT NOT NULL,
          z1 TEXT NOT NULL,
          z2 TEXT NOT NULL,
          p1 TEXT NOT NULL,
          p2 TEXT NOT NULL,
          p TEXT NOT NULL,
          q TEXT NOT NULL,
          eligibility_basis TEXT NOT NULL,
          evidence_labels TEXT NOT NULL,
          reference_status TEXT,
          reference_chr TEXT,
          reference_bp TEXT,
          reference_a1 TEXT,
          reference_a2 TEXT,
          clump_status TEXT,
          failure_reason TEXT,
          UNIQUE(pair_id, evidence_id),
          UNIQUE(pair_id, snp)
        );
        CREATE INDEX eligible_snp_index ON eligible(snp);
        CREATE TABLE bim_hits (
          snp TEXT NOT NULL, chr INTEGER NOT NULL, bp INTEGER NOT NULL,
          a1 TEXT NOT NULL, a2 TEXT NOT NULL
        );
        CREATE INDEX bim_hits_snp_index ON bim_hits(snp);
        """
    )
    return connection


def _stream_pair_rows(
    root: Path, pair: PairInput,
    callback: Callable[[int, Mapping[str, str], float, float], None] | None = None,
    output: Any | None = None,
) -> dict[str, Any]:
    role = PAIR_IDENTITIES[pair.pair_id][2]
    provenance = read_json(root, pair.provenance.relative_to(root))
    _validate_pair_provenance(pair, provenance)
    rows = failures = primary = pairwise = bh = 0
    minimum_p = 1.0
    with stable_gzip_text(root, pair.ledger.relative_to(root), pair.ledger_identity) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if list(reader.fieldnames or []) != PLACO_LEDGER_FIELDS:
            raise CollationError(f"PLACO {pair.pair_id} ledger schema drifted")
        for rows, row in enumerate(reader, start=1):
            if None in row or any(value is None for value in row.values()):
                raise CollationError(f"ragged PLACO {pair.pair_id} row {rows}")
            if row["analysis_id"] != "track-b-v1.0-pleiotropy" or row["pair_id"] != pair.pair_id:
                raise CollationError(f"PLACO row escaped frozen pair {pair.pair_id}")
            try:
                chromosome = int(row["CHR"])
                position = int(row["BP"])
                family_n = int(row["within_pair_family_n"])
            except ValueError as error:
                raise CollationError(f"invalid PLACO coordinate/family value at {pair.pair_id}:{rows}") from error
            p_value = _probability(row["P_PLACO_PLUS"], "PLACO P")
            q_value = _probability(row["PLACO_BH_Q"], "PLACO BH q")
            for field in ("Z1", "Z2"):
                _number(row[field], f"{pair.pair_id} {field}")
            for field in ("P1", "P2"):
                _probability(row[field], f"{pair.pair_id} {field}")
            snp, a1, a2 = row["SNP"], row["A1"].upper(), row["A2"].upper()
            if (
                chromosome not in range(1, 23) or position < 1 or family_n != pair.expected_rows
                or not RSID_RE.fullmatch(snp) or a1 not in BASES or a2 not in BASES
                or a1 == a2 or frozenset((a1, a2)) in PALINDROMIC
            ):
                raise CollationError(f"invalid PLACO variant identity at {pair.pair_id}:{rows}")
            status = row["analysis_status"]
            if status == "NUMERICAL_FAILURE_P_SET_TO_ONE":
                if (
                    p_value != 1 or q_value != 1
                    or row["T_PLACO_PLUS"] not in {"", "NA", "NaN", "nan"}
                    or not row["numerical_error"].strip() or row["numerical_error"] == "NA"
                ):
                    raise CollationError(f"invalid retained numerical failure at {pair.pair_id}:{rows}")
                failures += 1
            elif status == "TESTED":
                _number(row["T_PLACO_PLUS"], "tested PLACO statistic")
                if row["numerical_error"] not in {"", "NA"}:
                    raise CollationError(f"tested PLACO row has a numerical error at {pair.pair_id}:{rows}")
            else:
                raise CollationError(f"unknown PLACO row status at {pair.pair_id}:{rows}")
            minimum_p = min(minimum_p, p_value)
            primary += int(pair.pair_id != "CONTROL" and p_value <= 2.5e-8)
            pairwise += int(p_value <= 5e-8)
            bh += int(q_value <= 0.05)
            if callback is not None:
                callback(rows, row, p_value, q_value)
            if output is not None:
                output.write(_encode_tsv_row(
                    PLACO_07_FIELDS,
                    {"analysis_id": row["analysis_id"], "pair_id": pair.pair_id,
                     "family_role": role, **{field: row[field] for field in PLACO_LEDGER_FIELDS[2:]}},
                ))
    if rows != pair.expected_rows:
        raise CollationError(
            f"PLACO {pair.pair_id} ledger rows={rows} expected={pair.expected_rows}"
        )
    terminal = "COMPLETE_WITH_HITS" if primary or pairwise or bh else "TESTED_NO_HIT"
    counts = {
        "rows": rows, "failures": failures, "primary": primary, "pairwise": pairwise,
        "bh": bh, "minimum_p": minimum_p, "terminal_status": terminal,
    }
    expected_counts = provenance["complete_family_counts"]
    if any(expected_counts.get(key) != value for key, value in counts.items()):
        raise CollationError(f"PLACO {pair.pair_id} complete-family counts drifted")
    return counts


def index_complete_placo_family(
    root: Path, pairs: Sequence[PairInput], connection: sqlite3.Connection,
) -> tuple[dict[str, dict[str, Any]], int]:
    if [pair.pair_id for pair in pairs] != list(PAIR_ORDER):
        raise CollationError("PLACO input family is not exact ordered A/B/CONTROL")
    family_counts: dict[str, dict[str, Any]] = {}
    projected_07_bytes = len(_tsv_header(PLACO_07_FIELDS))
    insert_sql = (
        "INSERT INTO eligible (pair_id,family_role,source_row_index,evidence_id,snp,chr,bp,"
        "a1,a2,z1,z2,p1,p2,p,q,eligibility_basis,evidence_labels) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    )
    for pair in pairs:
        role = PAIR_IDENTITIES[pair.pair_id][2]
        pending: list[tuple[object, ...]] = []

        def consume(index: int, row: Mapping[str, str], p_value: float, q_value: float) -> None:
            nonlocal projected_07_bytes
            output_row = {
                "analysis_id": row["analysis_id"], "pair_id": pair.pair_id,
                "family_role": role, **{field: row[field] for field in PLACO_LEDGER_FIELDS[2:]},
            }
            projected_07_bytes += len(_encode_tsv_row(PLACO_07_FIELDS, output_row))
            basis = eligibility_basis(pair.pair_id, p_value, q_value)
            if row["analysis_status"] != "TESTED" or basis is None:
                return
            chromosome, position = int(row["CHR"]), int(row["BP"])
            evidence_id = f"PLACO:{pair.pair_id}:{row['SNP']}:{chromosome}:{position}"
            labels = evidence_labels(pair.pair_id, p_value, q_value)
            if not labels:
                raise CollationError("eligible PLACO signal lacks its method-specific label")
            pending.append((
                pair.pair_id, role, index, evidence_id, row["SNP"], chromosome, position,
                row["A1"].upper(), row["A2"].upper(), row["Z1"], row["Z2"],
                row["P1"], row["P2"], row["P_PLACO_PLUS"], row["PLACO_BH_Q"],
                basis, ";".join(labels),
            ))
            if len(pending) >= 10_000:
                try:
                    connection.executemany(insert_sql, pending)
                except sqlite3.IntegrityError as error:
                    raise CollationError(
                        f"eligible PLACO variant/evidence identity is duplicated for {pair.pair_id}"
                    ) from error
                pending.clear()

        family_counts[pair.pair_id] = _stream_pair_rows(root, pair, consume)
        if pending:
            try:
                connection.executemany(insert_sql, pending)
            except sqlite3.IntegrityError as error:
                raise CollationError(
                    f"eligible PLACO variant/evidence identity is duplicated for {pair.pair_id}"
                ) from error
        connection.commit()
    return family_counts, projected_07_bytes


def write_complete_07(
    root: Path, pairs: Sequence[PairInput], destination: Path,
    expected_bytes: int,
) -> dict[str, object]:
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    if destination.exists() or destination.is_symlink() or temporary.exists() or temporary.is_symlink():
        raise CollationError("07 staging destination or temporary already exists")
    digest = hashlib.sha256()

    class HashingWriter:
        def __init__(self, handle: Any):
            self.handle = handle
            self.bytes = 0

        def write(self, payload: bytes) -> None:
            self.handle.write(payload)
            digest.update(payload)
            self.bytes += len(payload)

    try:
        with temporary.open("xb") as raw:
            writer = HashingWriter(raw)
            writer.write(_tsv_header(PLACO_07_FIELDS))
            for pair in pairs:
                _stream_pair_rows(root, pair, output=writer)
            raw.flush()
            os.fsync(raw.fileno())
        if writer.bytes != expected_bytes:
            raise CollationError(
                f"exact 07 byte projection drifted: observed={writer.bytes} expected={expected_bytes}"
            )
        source = os.lstat(temporary)
        try:
            os.link(temporary, destination, follow_symlinks=False)
        except FileExistsError as error:
            raise CollationError("07 staging destination appeared concurrently") from error
        published = os.lstat(destination)
        if (
            not stat.S_ISREG(published.st_mode)
            or (source.st_dev, source.st_ino) != (published.st_dev, published.st_ino)
        ):
            raise CollationError("07 staging no-replace link changed the source inode")
        _fsync_directory(destination.parent)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    identity = {
        "path": destination.name, "bytes": destination.stat().st_size,
        "sha256": digest.hexdigest(),
    }
    if identity["bytes"] != expected_bytes:
        raise CollationError("07 staging artifact size drifted after atomic rename")
    return identity


class BloomFilter:
    """Small bounded accelerator for the disk-backed eligible-rsID table."""

    def __init__(self, expected_items: int):
        target_bits = max(1 << 20, max(1, expected_items) * 16)
        power = min(29, max(20, (target_bits - 1).bit_length()))
        self.bit_count = 1 << power
        self.mask = self.bit_count - 1
        self.bits = bytearray(self.bit_count // 8)

    @staticmethod
    def _hashes(value: str) -> tuple[int, int, int]:
        digest = hashlib.blake2b(value.encode("ascii"), digest_size=16).digest()
        first = int.from_bytes(digest[:8], "little")
        second = int.from_bytes(digest[8:], "little") | 1
        return first, first + second, first + 2 * second

    def add(self, value: str) -> None:
        for raw in self._hashes(value):
            bit = raw & self.mask
            self.bits[bit >> 3] |= 1 << (bit & 7)

    def __contains__(self, value: str) -> bool:
        return all(
            self.bits[(raw & self.mask) >> 3] & (1 << ((raw & self.mask) & 7))
            for raw in self._hashes(value)
        )


@contextmanager
def stable_plain_text(
    root: Path, relative: Path | str, expected: Mapping[str, Any] | None = None,
    *, allow_empty: bool = False,
) -> Iterator[TextIO]:
    path = _safe_path(root, relative, "text artifact")
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags)
    text_handle: TextIO | None = None
    try:
        named_before = os.lstat(path)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(named_before.st_mode) or not stat.S_ISREG(opened.st_mode)
            or (opened.st_size <= 0 and not allow_empty)
            or (named_before.st_dev, named_before.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise CollationError(f"text artifact is not one nonempty regular file: {path}")
        digest = hashlib.sha256()
        while True:
            chunk = os.read(descriptor, 4 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
        record = {"path": _relative(root, path), "bytes": opened.st_size, "sha256": digest.hexdigest()}
        if expected is not None and any(record.get(key) != expected.get(key) for key in record):
            raise CollationError(f"text artifact differs from its sealed identity: {path}")
        os.lseek(descriptor, 0, os.SEEK_SET)
        text_handle = io.TextIOWrapper(os.fdopen(os.dup(descriptor), "rb"), encoding="utf-8", newline="")
        yield text_handle
        while text_handle.read(1024 * 1024):
            pass
        text_handle.close()
        text_handle = None
        after = os.fstat(descriptor)
        named_after = os.lstat(path)
        if len({_stat_tuple(item) for item in (named_before, opened, after, named_after)}) != 1:
            raise CollationError(f"text artifact changed while streaming: {path}")
    except (OSError, UnicodeError) as error:
        raise CollationError(f"could not stream text artifact {path}: {error}") from error
    finally:
        if text_handle is not None:
            text_handle.close()
        os.close(descriptor)


def exact_reference_match(
    input_chr: int, input_bp: int, input_a1: str, input_a2: str,
    hits: Sequence[tuple[int, int, str, str]],
) -> tuple[str, str, str, str, str, str]:
    """Classify exact rsID+coordinate+unordered-nonpalindromic-allele coverage."""
    if not hits:
        return "REFERENCE_ABSENT", "NA", "NA", "NA", "NA", "RSID_ABSENT_FROM_SEALED_BIM"
    if len(hits) != 1:
        return (
            "REFERENCE_RSID_AMBIGUOUS", "MULTIPLE", "MULTIPLE", "MULTIPLE", "MULTIPLE",
            "RSID_HAS_MULTIPLE_RECORDS_IN_SEALED_BIM",
        )
    chromosome, position, allele1, allele2 = hits[0]
    reference_fields = str(chromosome), str(position), allele1, allele2
    if chromosome != input_chr or position != input_bp:
        return (
            "REFERENCE_COORDINATE_MISMATCH", *reference_fields,
            "RSID_PRESENT_BUT_CHR_OR_BP_DIFFERS_FROM_PLACO",
        )
    ref_pair = frozenset((allele1, allele2))
    input_pair = frozenset((input_a1, input_a2))
    if (
        allele1 not in BASES or allele2 not in BASES or allele1 == allele2
        or ref_pair in PALINDROMIC
    ):
        return (
            "REFERENCE_ALLELES_INVALID_OR_PALINDROMIC", *reference_fields,
            "SEALED_BIM_ALLELE_PAIR_IS_NOT_DISTINCT_NONPALINDROMIC_ACGT",
        )
    if ref_pair != input_pair:
        return (
            "REFERENCE_ALLELE_MISMATCH", *reference_fields,
            "UNORDERED_ALLELE_PAIR_DIFFERS_FROM_SEALED_BIM",
        )
    return "MATCHED", *reference_fields, "NONE"


def match_eligible_to_bim(
    root: Path, connection: sqlite3.Connection, bim_relative: Path | str,
    expected_identity: Mapping[str, Any] | None = None,
) -> dict[str, int]:
    eligible_count = int(connection.execute("SELECT COUNT(*) FROM eligible").fetchone()[0])
    bloom = BloomFilter(eligible_count)
    for (snp,) in connection.execute("SELECT DISTINCT snp FROM eligible"):
        bloom.add(str(snp))
    pending: list[tuple[str, int, int, str, str]] = []
    bim_rows = candidate_rows = 0
    with stable_plain_text(root, bim_relative, expected_identity) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) != 6:
                raise CollationError(f"BIM schema drift at line {line_number}")
            chromosome_raw, snp, cm_raw, position_raw, allele1, allele2 = fields
            try:
                chromosome, position = int(chromosome_raw), int(position_raw)
                genetic_distance = float(cm_raw)
            except ValueError as error:
                raise CollationError(f"BIM numeric field drift at line {line_number}") from error
            allele1, allele2 = allele1.upper(), allele2.upper()
            if (
                chromosome not in range(1, 24) or position < 1
                or not math.isfinite(genetic_distance) or genetic_distance < 0
                or not snp or snp == "." or allele1 not in BASES or allele2 not in BASES
                or allele1 == allele2
            ):
                raise CollationError(f"BIM identity field drift at line {line_number}")
            bim_rows += 1
            if snp in bloom:
                pending.append((snp, chromosome, position, allele1, allele2))
                candidate_rows += 1
                if len(pending) >= 20_000:
                    connection.executemany("INSERT INTO bim_hits VALUES (?,?,?,?,?)", pending)
                    pending.clear()
        if pending:
            connection.executemany("INSERT INTO bim_hits VALUES (?,?,?,?,?)", pending)
    updates: list[tuple[str, str, str, str, str, str, int]] = []
    counts: dict[str, int] = defaultdict(int)
    cursor = connection.execute(
        "SELECT row_id,chr,bp,a1,a2,snp FROM eligible ORDER BY row_id"
    )
    for row_id, chromosome, position, allele1, allele2, snp in cursor:
        hits = [
            (int(hit_chr), int(hit_bp), str(hit_a1), str(hit_a2))
            for hit_chr, hit_bp, hit_a1, hit_a2 in connection.execute(
                "SELECT chr,bp,a1,a2 FROM bim_hits WHERE snp=? ORDER BY chr,bp,a1,a2", (snp,),
            )
        ]
        result = exact_reference_match(
            int(chromosome), int(position), str(allele1), str(allele2), hits,
        )
        status, ref_chr, ref_bp, ref_a1, ref_a2, reason = result
        counts[status] += 1
        updates.append((status, ref_chr, ref_bp, ref_a1, ref_a2, reason, int(row_id)))
        if len(updates) >= 10_000:
            connection.executemany(
                "UPDATE eligible SET reference_status=?,reference_chr=?,reference_bp=?,"
                "reference_a1=?,reference_a2=?,failure_reason=? WHERE row_id=?", updates,
            )
            updates.clear()
    if updates:
        connection.executemany(
            "UPDATE eligible SET reference_status=?,reference_chr=?,reference_bp=?,"
            "reference_a1=?,reference_a2=?,failure_reason=? WHERE row_id=?", updates,
        )
    connection.commit()
    if sum(counts.values()) != eligible_count:
        raise CollationError("eligible/reference coverage denominator did not reconcile")
    return {"eligible_rows": eligible_count, "bim_rows_streamed": bim_rows,
            "bloom_candidate_rows": candidate_rows, **dict(counts)}


def disk_preflight(
    root: Path, exact_07_bytes: int, eligible_rows: int, *, free_bytes: int | None = None,
) -> dict[str, int | str]:
    """Require exact 07 bytes plus a conservative bounded auxiliary allowance."""
    if exact_07_bytes <= 0 or eligible_rows < 0:
        raise CollationError("invalid disk-preflight family dimensions")
    free = shutil.disk_usage(root).free if free_bytes is None else int(free_bytes)
    auxiliary = 256 * 1024**2 + eligible_rows * 4096
    required = exact_07_bytes + auxiliary
    if free < required:
        raise CollationError(
            "BLOCKED_BY_COMPUTE insufficient storage for complete plain 07 plus clumping/QC family: "
            f"free={free} required={required} exact_07_bytes={exact_07_bytes}"
        )
    return {
        "free_bytes_before_07": free,
        "exact_projected_07_bytes": exact_07_bytes,
        "auxiliary_allowance_bytes": auxiliary,
        "required_bytes": required,
        "headroom_bytes": free - required,
        "rule": "EXACT_PROJECTED_07_PLUS_256_MIB_PLUS_4096_BYTES_PER_ELIGIBLE_SIGNAL",
    }


def parse_plink_clumped_text(text: str) -> list[dict[str, Any]]:
    """Parse standard PLINK 1.9 .clumped output without inferring loci by distance."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise CollationError("PLINK clump output is empty")
    header = lines[0].split()
    required = {"CHR", "SNP", "BP", "P", "SP2"}
    if not required.issubset(header) or len(header) != len(set(header)):
        raise CollationError(f"PLINK clump schema drifted: {header}")
    records: list[dict[str, Any]] = []
    seen_leads: set[str] = set()
    assigned: set[str] = set()
    for line_number, line in enumerate(lines[1:], start=2):
        values = line.split()
        if len(values) != len(header):
            raise CollationError(f"ragged PLINK clump row {line_number}")
        raw = dict(zip(header, values, strict=True))
        try:
            chromosome, position = int(raw["CHR"]), int(raw["BP"])
            p_value = float(raw["P"])
        except ValueError as error:
            raise CollationError(f"invalid PLINK clump row {line_number}") from error
        lead = raw["SNP"]
        if (
            chromosome not in range(1, 23) or position < 1
            or not math.isfinite(p_value) or not 0 <= p_value <= 1
            or not lead or lead in seen_leads
        ):
            raise CollationError(f"invalid/duplicate PLINK clump lead at row {line_number}")
        members = [lead]
        if raw["SP2"] not in {"NONE", "NA", "."}:
            for token in raw["SP2"].split(","):
                match = re.fullmatch(r"(.+?)(?:\([0-9]+\))?", token)
                if match is None or not match.group(1):
                    raise CollationError(f"invalid PLINK SP2 token at row {line_number}: {token}")
                members.append(match.group(1))
        if len(members) != len(set(members)) or assigned.intersection(members):
            raise CollationError("PLINK clump memberships overlap or contain duplicates")
        seen_leads.add(lead)
        assigned.update(members)
        records.append({
            "lead": lead, "chr": chromosome, "bp": position, "plink_p": p_value,
            "members": members, "raw": raw,
        })
    return records


def parse_pairwise_r2_text(
    text: str,
    leads: Mapping[str, tuple[int, int]],
    *, window_bp: int = 1_000_000, strict_upper_bound: float = 0.1,
) -> dict[str, Any]:
    """Verify every within-window lead pair is present and has r2 strictly below 0.1."""
    if strict_upper_bound != 0.1 or window_bp != 1_000_000:
        raise CollationError("pairwise lead-independence contract must remain r2<0.1 within 1 Mb")
    lead_ids = set(leads)
    expected = {
        tuple(sorted((left, right)))
        for offset, left in enumerate(sorted(lead_ids))
        for right in sorted(lead_ids)[offset + 1:]
        if leads[left][0] == leads[right][0]
        and abs(leads[left][1] - leads[right][1]) <= window_bp
    }
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        if expected:
            raise CollationError("PLINK pairwise-r2 output is empty for within-window lead pairs")
        return {"expected_pairs": 0, "observed_pairs": 0, "maximum_r2": None,
                "status": "PASS_NO_WITHIN_WINDOW_LEAD_PAIR"}
    header = lines[0].split()
    required = {"CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B", "SNP_B", "R2"}
    if not required.issubset(header) or len(header) != len(set(header)):
        raise CollationError(f"PLINK pairwise-r2 schema drifted: {header}")
    observed: dict[tuple[str, str], float] = {}
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split()
        if len(fields) != len(header):
            raise CollationError(f"ragged PLINK pairwise-r2 row {line_number}")
        row = dict(zip(header, fields, strict=True))
        left, right = row["SNP_A"], row["SNP_B"]
        if left == right or left not in lead_ids or right not in lead_ids:
            continue
        try:
            r2 = float(row["R2"])
            chr_a, chr_b = int(row["CHR_A"]), int(row["CHR_B"])
            bp_a, bp_b = int(row["BP_A"]), int(row["BP_B"])
        except ValueError as error:
            raise CollationError(f"invalid PLINK pairwise-r2 row {line_number}") from error
        if not math.isfinite(r2) or not 0 <= r2 <= 1:
            raise CollationError(f"out-of-range PLINK r2 at row {line_number}")
        if (chr_a, bp_a) != leads[left] or (chr_b, bp_b) != leads[right]:
            raise CollationError("PLINK pairwise-r2 coordinates differ from sealed BIM/lead identities")
        key = tuple(sorted((left, right)))
        if key not in expected:
            continue
        if key in observed and not math.isclose(observed[key], r2, rel_tol=1e-12, abs_tol=1e-15):
            raise CollationError("PLINK pairwise-r2 duplicate directions disagree")
        observed[key] = r2
    missing = expected - set(observed)
    if missing:
        raise CollationError(
            "PLINK pairwise-r2 audit omitted required within-window lead pair(s): "
            + ",".join("/".join(pair) for pair in sorted(missing)[:5])
        )
    violations = {pair: value for pair, value in observed.items() if value >= strict_upper_bound}
    if violations:
        pair, value = sorted(violations.items())[0]
        raise CollationError(
            f"PLINK clump lead independence failed strict r2<0.1: {pair[0]}/{pair[1]} r2={value}"
        )
    maximum = max(observed.values()) if observed else None
    return {
        "expected_pairs": len(expected), "observed_pairs": len(observed),
        "maximum_r2": maximum,
        "status": "PASS_STRICT_R2_LT_0.1_WITHIN_1_MB",
    }


def _write_stage_text(path: Path, content: str) -> None:
    if path.exists() or path.is_symlink():
        raise CollationError(f"private stage artifact already exists: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())


def _default_runner(argv: Sequence[str], cwd: Path) -> subprocess.CompletedProcess[str]:
    """Run one PLINK command in an owned PGID with live aggregate RSS sampling."""
    started = time.monotonic()
    peak = 0
    child_seen = False
    process: subprocess.Popen[bytes] | None = None
    with tempfile.TemporaryFile(mode="w+b") as stdout_file, tempfile.TemporaryFile(mode="w+b") as stderr_file:
        try:
            process = subprocess.Popen(
                list(argv), cwd=cwd, stdout=stdout_file, stderr=stderr_file,
                start_new_session=True,
            )
            while True:
                sample = process_group_and_parent_rss_bytes(process.pid, os.getpid())
                peak = max(peak, sample["aggregate_bytes"])
                child_seen = child_seen or sample["process_group_bytes"] > 0
                if process.poll() is not None:
                    break
                time.sleep(0.05)
            lingering = process_group_rss_bytes(process.pid)
            if lingering > 0:
                confirmed = terminate_owned_process_group(process)
                returncode = 70
                stderr_file.write(
                    (
                        "\nPLINK owned process group outlived its leader; "
                        f"termination_confirmed={str(confirmed).upper()}\n"
                    ).encode("utf-8")
                )
            else:
                returncode = int(process.wait(timeout=1.0))
        except BaseException:
            if process is not None:
                terminate_owned_process_group(process)
            raise
        stdout_file.seek(0)
        stderr_file.seek(0)
        stdout = stdout_file.read().decode("utf-8", errors="replace")
        stderr = stderr_file.read().decode("utf-8", errors="replace")
    if not child_seen:
        peak = max(
            peak,
            normalize_ru_maxrss(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
            + normalize_ru_maxrss(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),
        )
    completed = subprocess.CompletedProcess(list(argv), returncode, stdout=stdout, stderr=stderr)
    completed.peak_aggregate_rss_bytes = peak
    completed.rss_measurement_method = (
        "LIVE_PARENT_PLUS_OWNED_FRESH_SESSION_PGID_PS_SAMPLED"
        if child_seen else "NORMALIZED_SELF_PLUS_CHILD_RUSAGE_CONSERVATIVE_FAST_EXIT_FALLBACK"
    )
    completed.runtime_seconds = time.monotonic() - started
    return completed


def normalize_ru_maxrss(value: int | float, *, platform: str | None = None) -> int:
    """Normalize getrusage maxrss to bytes (Darwin reports bytes; Linux KiB)."""
    observed = int(value)
    if observed < 0:
        raise CollationError("negative ru_maxrss is invalid")
    target = sys.platform if platform is None else platform
    return observed if target.startswith("darwin") else observed * 1024


def parse_ps_rss_snapshot(
    text: str, process_group_id: int, parent_pid: int,
) -> dict[str, int]:
    if process_group_id <= 0 or parent_pid <= 0:
        raise CollationError("invalid process-group/parent identity for RSS sampling")
    group_kib = parent_kib = 0
    seen_pids: set[int] = set()
    for line in text.splitlines():
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            pid, pgid, rss_kib = (int(value) for value in fields)
        except ValueError:
            continue
        if pid <= 0 or rss_kib < 0 or pid in seen_pids:
            continue
        seen_pids.add(pid)
        if pgid == process_group_id:
            group_kib += rss_kib
        if pid == parent_pid and pgid != process_group_id:
            parent_kib = rss_kib
    return {
        "process_group_bytes": group_kib * 1024,
        "parent_bytes": parent_kib * 1024,
        "aggregate_bytes": (group_kib + parent_kib) * 1024,
    }


def process_group_and_parent_rss_bytes(
    process_group_id: int, parent_pid: int,
) -> dict[str, int]:
    result = subprocess.run(
        ["ps", "-axo", "pid=,pgid=,rss="], check=False,
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise CollationError(f"ps failed while sampling PLINK RSS: {result.stderr.strip()}")
    return parse_ps_rss_snapshot(result.stdout, process_group_id, parent_pid)


def process_group_rss_bytes(process_group_id: int) -> int:
    return process_group_and_parent_rss_bytes(process_group_id, os.getpid())["process_group_bytes"]


def terminate_owned_process_group(process: subprocess.Popen[Any]) -> bool:
    """Best-effort bounded termination; true only once the exact PGID disappears."""
    process_group_id = process.pid

    def group_exists() -> bool:
        process.poll()
        try:
            os.killpg(process_group_id, 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True

    try:
        os.killpg(process_group_id, signal.SIGTERM)
    except ProcessLookupError:
        return True
    except PermissionError:
        return False
    deadline = time.monotonic() + 0.5
    while time.monotonic() < deadline:
        if not group_exists():
            return True
        time.sleep(0.05)
    try:
        os.killpg(process_group_id, signal.SIGKILL)
    except ProcessLookupError:
        return True
    except PermissionError:
        return False
    deadline = time.monotonic() + 1.0
    while time.monotonic() < deadline:
        if not group_exists():
            return True
        time.sleep(0.05)
    return False


def _stage_identity(
    path: Path, package_relative: Path, *, allow_empty: bool = False,
) -> dict[str, object]:
    digest = hashlib.sha256()
    before = os.lstat(path)
    if not stat.S_ISREG(before.st_mode) or (before.st_size == 0 and not allow_empty):
        raise CollationError(f"stage artifact is not a permitted regular file: {path}")
    with path.open("rb") as handle:
        opened = os.fstat(handle.fileno())
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
        after = os.fstat(handle.fileno())
    current = os.lstat(path)
    if len({_stat_tuple(item) for item in (before, opened, after, current)}) != 1:
        raise CollationError(f"PLINK artifact changed while hashing: {path}")
    return {"path": str(package_relative), "bytes": current.st_size, "sha256": digest.hexdigest()}


def _run_plink(
    argv: Sequence[str], cwd: Path, phase: str,
    runner: Callable[[Sequence[str], Path], subprocess.CompletedProcess[str]],
) -> dict[str, Any]:
    started = time.monotonic()
    result = runner(argv, cwd)
    elapsed = float(getattr(result, "runtime_seconds", time.monotonic() - started))
    transcript = cwd / f"{phase}.subprocess.log"
    _write_stage_text(
        transcript,
        "ARGV_JSON\t" + json.dumps(list(argv), separators=(",", ":")) + "\n"
        + "RETURN_CODE\t" + str(result.returncode) + "\n"
        + "STDOUT\n" + (result.stdout or "") + "\nSTDERR\n" + (result.stderr or ""),
    )
    if result.returncode != 0:
        raise CollationError(
            f"pinned PLINK {phase} failed with status {result.returncode}; transcript={transcript}"
        )
    return {
        "argv": list(argv), "runtime_seconds": elapsed, "returncode": result.returncode,
        "transcript": transcript,
        "peak_aggregate_rss_bytes": int(getattr(result, "peak_aggregate_rss_bytes", 0)),
        "rss_measurement_method": str(
            getattr(result, "rss_measurement_method", "CUSTOM_RUNNER_DID_NOT_REPORT_RSS")
        ),
    }


def _write_plink_input_files(
    connection: sqlite3.Connection, pair_id: str, association: Path, extract: Path,
) -> int:
    """Stream the eligible exact-reference pair family into PLINK inputs."""
    if any(path.exists() or path.is_symlink() for path in (association, extract)):
        raise CollationError("private PLINK input already exists")
    count = 0
    with association.open("xb") as association_handle, extract.open("xb") as extract_handle:
        association_handle.write(b"SNP\tP\n")
        cursor = connection.execute(
            "SELECT snp,p FROM eligible WHERE pair_id=? AND reference_status='MATCHED' "
            "ORDER BY CAST(p AS REAL),source_row_index", (pair_id,),
        )
        for count, (snp, p_value) in enumerate(cursor, start=1):
            association_handle.write(f"{snp}\t{p_value}\n".encode("ascii"))
            extract_handle.write(f"{snp}\n".encode("ascii"))
        association_handle.flush()
        extract_handle.flush()
        os.fsync(association_handle.fileno())
        os.fsync(extract_handle.fileno())
    _fsync_directory(association.parent)
    return count


def run_plink_clump_for_pair(
    connection: sqlite3.Connection, pair_id: str, plink_path: Path,
    reference_prefix: Path, directory: Path,
    *, root: Path = ROOT,
    runner: Callable[[Sequence[str], Path], subprocess.CompletedProcess[str]] = _default_runner,
) -> dict[str, Any]:
    """Run real pair-scoped PLINK clumping and an independent pairwise-r2 audit."""
    if pair_id not in PAIR_ORDER:
        raise CollationError(f"unknown clumping pair: {pair_id}")
    directory.mkdir(parents=False, exist_ok=False)
    association = directory / "eligible.assoc.tsv"
    extract = directory / "eligible.snps.txt"
    matched_rows = _write_plink_input_files(connection, pair_id, association, extract)
    if matched_rows == 0:
        return {
            "pair_id": pair_id, "matched_rows": 0, "leads": [], "clumps": [],
            "commands": [], "pairwise": {
                "expected_pairs": 0, "observed_pairs": 0, "maximum_r2": None,
                "status": "PASS_NO_REFERENCE_MATCHED_ELIGIBLE_SIGNAL",
            },
        }
    out = directory / "clump"
    clump_command = [
        str(plink_path), "--bfile", str(reference_prefix), "--extract", str(extract),
        "--clump", str(association), "--clump-snp-field", "SNP", "--clump-field", "P",
        "--clump-p1", "1", "--clump-p2", "1", "--clump-r2", "0.1",
        "--clump-kb", "1000", "--threads", "1", "--out", str(out),
    ]
    commands = [_run_plink(clump_command, directory, "clump", runner)]
    clumped_path = out.with_suffix(".clumped")
    if not clumped_path.is_file() or clumped_path.is_symlink():
        raise CollationError(f"PLINK did not publish a regular .clumped file for {pair_id}")
    clumped_identity = _stage_identity(
        clumped_path, clumped_path.relative_to(root), allow_empty=False,
    )
    with stable_plain_text(root, clumped_path.relative_to(root), clumped_identity) as handle:
        clumps = parse_plink_clumped_text(handle.read())
    for clump in clumps:
        lead = clump["lead"]
        lead_rows = connection.execute(
            "SELECT chr,bp,p FROM eligible WHERE pair_id=? AND snp=? "
            "AND reference_status='MATCHED'", (pair_id, lead),
        ).fetchall()
        if len(lead_rows) != 1:
            raise CollationError(f"PLINK returned a lead outside eligible matched {pair_id} rows")
        chromosome, position, source_p = lead_rows[0]
        if (clump["chr"], clump["bp"]) != (int(chromosome), int(position)):
            raise CollationError("PLINK lead coordinate differs from exact sealed BIM match")
        if not math.isclose(
            clump["plink_p"], float(source_p), rel_tol=1e-5, abs_tol=1e-300,
        ):
            raise CollationError("PLINK clump lead P differs from actual PLACO P")
        for member in clump["members"]:
            member_rows = connection.execute(
                "SELECT p FROM eligible WHERE pair_id=? AND snp=? "
                "AND reference_status='MATCHED'", (pair_id, member),
            ).fetchall()
            if len(member_rows) != 1:
                raise CollationError("PLINK clump member escaped eligible matched signal family")
            if float(member_rows[0][0]) < float(source_p):
                raise CollationError("PLINK clump lead is not ordered by actual PLACO P")
    connection.execute(
        "UPDATE eligible SET clump_status='PLINK_CLUMP_UNACCOUNTED',"
        "failure_reason='PLINK_SUCCESS_BUT_ELIGIBLE_MATCHED_VARIANT_NOT_IN_CLUMP_OUTPUT' "
        "WHERE pair_id=? AND reference_status='MATCHED'", (pair_id,),
    )
    for clump in clumps:
        lead = clump["lead"]
        secondaries = [member for member in clump["members"] if member != lead]
        if secondaries:
            connection.executemany(
                "UPDATE eligible SET clump_status=?,failure_reason='NONE' WHERE pair_id=? AND snp=?",
                [(f"CLUMPED_SECONDARY_TO:{lead}", pair_id, member) for member in secondaries],
            )
        connection.execute(
            "UPDATE eligible SET clump_status='INDEPENDENT_LEAD_PENDING_R2_AUDIT',"
            "failure_reason='NONE' WHERE pair_id=? AND snp=?", (pair_id, lead),
        )
    missing = [
        str(row[0]) for row in connection.execute(
            "SELECT snp FROM eligible WHERE pair_id=? AND reference_status='MATCHED' "
            "AND clump_status='PLINK_CLUMP_UNACCOUNTED' ORDER BY source_row_index LIMIT 5",
            (pair_id,),
        )
    ]
    if missing:
        raise CollationError(
            "successful PLINK clumping omitted eligible exact-reference variants: "
            + ",".join(missing)
        )
    leads_path = directory / "leads.snps.txt"
    _write_stage_text(leads_path, "".join(f"{clump['lead']}\n" for clump in clumps))
    lead_coordinates = {
        clump["lead"]: tuple(map(int, connection.execute(
            "SELECT chr,bp FROM eligible WHERE pair_id=? AND snp=?",
            (pair_id, clump["lead"]),
        ).fetchone()))
        for clump in clumps
    }
    expected_pair_count = sum(
        left_chr == right_chr and abs(left_bp - right_bp) <= 1_000_000
        for offset, (left_chr, left_bp) in enumerate(lead_coordinates.values())
        for right_chr, right_bp in list(lead_coordinates.values())[offset + 1:]
    )
    if len(clumps) >= 2 and expected_pair_count:
        r2_out = directory / "lead_pairwise"
        r2_command = [
            str(plink_path), "--bfile", str(reference_prefix), "--extract", str(leads_path),
            "--r2", "--ld-snp-list", str(leads_path), "--ld-window", "99999999",
            "--ld-window-kb", "1000", "--ld-window-r2", "0",
            "--threads", "1", "--out", str(r2_out),
        ]
        commands.append(_run_plink(r2_command, directory, "pairwise_r2", runner))
        ld_path = r2_out.with_suffix(".ld")
        if not ld_path.is_file() or ld_path.is_symlink():
            raise CollationError(f"PLINK did not publish pairwise .ld output for {pair_id}")
        ld_identity = _stage_identity(ld_path, ld_path.relative_to(root), allow_empty=False)
        with stable_plain_text(root, ld_path.relative_to(root), ld_identity) as handle:
            pairwise = parse_pairwise_r2_text(handle.read(), lead_coordinates)
    else:
        pairwise = {
            "expected_pairs": 0, "observed_pairs": 0, "maximum_r2": None,
            "status": "PASS_NO_WITHIN_WINDOW_LEAD_PAIR",
        }
    connection.executemany(
        "UPDATE eligible SET clump_status='INDEPENDENT_LEAD_VERIFIED_R2_LT_0.1',"
        "failure_reason='NONE' WHERE pair_id=? AND snp=?",
        [(pair_id, clump["lead"]) for clump in clumps],
    )
    connection.commit()
    return {
        "pair_id": pair_id, "matched_rows": matched_rows, "missing_from_clump": [],
        "leads": [clump["lead"] for clump in clumps], "clumps": clumps,
        "commands": commands, "pairwise": pairwise,
    }


def _read_small_tsv(
    root: Path, relative: Path | str, *, required_fields: Iterable[str] = (),
) -> tuple[list[str], list[dict[str, str]]]:
    identity = stable_identity(root, relative)
    if int(identity["bytes"]) > 256 * 1024**2:
        raise CollationError(f"small structured TSV exceeds 256-MiB bound: {relative}")
    path = root / str(relative)
    try:
        with stable_plain_text(
            root, relative,
            {key: identity[key] for key in ("path", "bytes", "sha256")},
        ) as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = list(reader.fieldnames or [])
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise CollationError(f"unreadable TSV {relative}: {error}") from error
    if (
        not fields or len(fields) != len(set(fields))
        or not set(required_fields).issubset(fields)
        or any(None in row or any(value is None for value in row.values()) for row in rows)
    ):
        raise CollationError(f"TSV schema/row shape drifted: {relative}")
    return fields, rows


def validate_conjfdr_blocked(root: Path) -> dict[str, Any]:
    readiness = Path("results/track_b/pleiotropy/readiness_gate.tsv")
    _, rows = _read_small_tsv(
        root, readiness,
        required_fields={
            "component_id", "pair_scope", "implementation_gate", "software_gate",
            "blocked_result_semantics", "claim_status",
        },
    )
    matches = [row for row in rows if row["component_id"] == "CONJFDR"]
    if len(matches) != 1:
        raise CollationError("frozen readiness lacks exactly one conjFDR component")
    row = matches[0]
    if (
        row["pair_scope"] != "A;B;CONTROL"
        or row["implementation_gate"]
        != "BLOCKED_BY_IMPLEMENTATION_TRACK_B_CONJFDR_RUNNER_AND_VALIDATOR"
        or not row["software_gate"].startswith("BLOCKED_BY_SOFTWARE")
        or "NOT_TESTED_METHOD_BLOCKED" not in row["blocked_result_semantics"]
        or "NO_ONLY_LABEL_IF_OTHER_METHOD_BLOCKED" not in row["blocked_result_semantics"]
        or row["claim_status"] != "NO_SCIENTIFIC_RESULT"
    ):
        raise CollationError("frozen readiness no longer proves blocked conjFDR semantics")
    root_relative = Path("results/track_b/pleiotropy/results/conjfdr")
    directory = root / root_relative
    if directory.exists() or directory.is_symlink():
        if directory.is_symlink() or not directory.is_dir():
            raise CollationError("conjFDR result namespace is a symlink/special file")
        artifacts = [path for path in directory.rglob("*") if path.is_file() or path.is_symlink()]
        if artifacts:
            raise CollationError(
                "conjFDR artifacts exist without a supported signal-level cross-method collation; "
                "refusing to call the method blocked"
            )
    return {
        "family_state": CONJFDR_BLOCKED,
        "pair_states": {pair: CONJFDR_BLOCKED for pair in PAIR_ORDER},
        "evidence_rows": 0,
        "frozen_readiness": portable_identity(root, readiness),
        "semantics": "NOT_NULL_NOT_NEGATIVE_AND_NO_PLACO_ONLY_OR_BOTH_LABEL",
    }


def load_lava_diagnostic_index(root: Path) -> dict[tuple[str, str], dict[str, str]]:
    """Load only small, already validated 142 accounting rows for annotations."""
    _, bivariate = _read_small_tsv(
        root, "results/track_b/lava/local_bivariate_results.tsv",
        required_fields={"LOC", "pair_id", "locus_status", "analysis_status", "terminal_scientific_state"},
    )
    _, failures = _read_small_tsv(
        root, "results/track_b/lava/local_failures.tsv",
        required_fields={"LOC", "pair_id", "stage", "analysis_status", "terminal_scientific_state"},
    )
    index: dict[tuple[str, str], dict[str, str]] = {}
    for row in bivariate:
        pair, locus = row["pair_id"], row["LOC"]
        if pair not in PAIR_ORDER or not locus.isdigit():
            raise CollationError("LAVA bivariate diagnostic has invalid pair/locus identity")
        key = pair, locus
        if key in index:
            raise CollationError("LAVA bivariate diagnostic duplicates a pair/locus row")
        if row["terminal_scientific_state"] != "TERMINAL_FAILED_QC_DIAGNOSTIC_ONLY":
            raise CollationError("LAVA bivariate row lost terminal failed-QC semantics")
        index[key] = {
            "locus_status": row["locus_status"],
            "bivariate_status": row["analysis_status"],
            "conditional_failures": "NONE_REPORTED",
        }
    conditional: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in failures:
        if row["terminal_scientific_state"] != "TERMINAL_FAILED_QC_DIAGNOSTIC_ONLY":
            raise CollationError("LAVA failure row lost terminal failed-QC semantics")
        if row["stage"] not in {"CONDITIONAL_MODEL", "CONDITIONAL_CHECKPOINT"}:
            continue
        if row["pair_id"] in PAIR_ORDER and row["LOC"].isdigit():
            conditional[(row["pair_id"], row["LOC"])].add(row["analysis_status"])
    for key, statuses in conditional.items():
        if key not in index:
            raise CollationError("conditional LAVA failure lacks its complete pair/locus grid row")
        index[key]["conditional_failures"] = ";".join(sorted(statuses))
    if len(index) != 2495 * len(PAIR_ORDER):
        raise CollationError(
            f"LAVA diagnostic pair/locus grid has {len(index)} rows, expected {2495 * len(PAIR_ORDER)}"
        )
    return index


def _stage_write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]],
) -> tuple[int, dict[str, object]]:
    if path.exists() or path.is_symlink():
        raise CollationError(f"private stage output already exists: {path}")
    count = 0
    digest = hashlib.sha256()
    with path.open("xb") as handle:
        header = _tsv_header(fields)
        handle.write(header)
        digest.update(header)
        for count, row in enumerate(rows, start=1):
            encoded = _encode_tsv_row(fields, row)
            handle.write(encoded)
            digest.update(encoded)
        handle.flush()
        os.fsync(handle.fileno())
    return count, {"path": path.name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def _coverage_rows(
    connection: sqlite3.Connection, *, pair_filter: str | None = None,
) -> Iterator[dict[str, str]]:
    fields = [
        "pair_id", "family_role", "evidence_id", "snp", "chr", "bp", "a1", "a2",
        "p", "q", "eligibility_basis", "reference_status", "reference_chr", "reference_bp",
        "reference_a1", "reference_a2", "clump_status", "failure_reason",
    ]
    sql = "SELECT " + ",".join(fields) + " FROM eligible"
    parameters: tuple[object, ...] = ()
    if pair_filter is not None:
        sql += " WHERE pair_id=?"
        parameters = (pair_filter,)
    sql += " ORDER BY CASE pair_id WHEN 'A' THEN 1 WHEN 'B' THEN 2 ELSE 3 END,source_row_index"
    for values in connection.execute(sql, parameters):
        row = dict(zip(fields, values, strict=True))
        yield {
            "analysis_id": "track-b-v1.0-pleiotropy",
            "pair_id": row["pair_id"], "family_role": row["family_role"],
            "evidence_id": row["evidence_id"], "SNP": row["snp"],
            "CHR": str(row["chr"]), "BP": str(row["bp"]), "A1": row["a1"], "A2": row["a2"],
            "P_PLACO_PLUS": row["p"], "PLACO_BH_Q": row["q"],
            "eligibility_basis": row["eligibility_basis"],
            "reference_status": row["reference_status"] or "UNCLASSIFIED",
            "reference_CHR": row["reference_chr"] or "NA",
            "reference_BP": row["reference_bp"] or "NA",
            "reference_A1": row["reference_a1"] or "NA",
            "reference_A2": row["reference_a2"] or "NA",
            "clump_status": row["clump_status"] or "NOT_CLUMPED_REFERENCE_COVERAGE_FAILURE",
            "failure_reason": row["failure_reason"] or "UNCLASSIFIED_FAILURE",
        }


def build_shared_locus_rows(
    connection: sqlite3.Connection,
    clump_results: Mapping[str, Mapping[str, Any]],
    block_index: Mapping[int, tuple[list[int], list[Block]]],
    chromosome_lengths: Mapping[str, int],
    lava_index: Mapping[tuple[str, str], Mapping[str, str]],
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    lead_records: list[dict[str, Any]] = []
    for pair_id in PAIR_ORDER:
        result = clump_results[pair_id]
        by_lead = {record["lead"]: record for record in result["clumps"]}
        for lead in result["leads"]:
            fields = [
                "evidence_id", "snp", "chr", "bp", "a1", "a2", "z1", "z2", "p1", "p2",
                "p", "q", "eligibility_basis", "evidence_labels", "reference_status",
            ]
            values = connection.execute(
                "SELECT " + ",".join(fields) + " FROM eligible WHERE pair_id=? AND snp=?",
                (pair_id, lead),
            ).fetchall()
            if len(values) != 1:
                raise CollationError("independent PLINK lead lacks one exact eligible source row")
            source = dict(zip(fields, values[0], strict=True))
            if source["reference_status"] != "MATCHED":
                raise CollationError("independent PLINK lead is not exact-matched to the sealed BIM")
            block = map_to_block(int(source["chr"]), int(source["bp"]), block_index)
            clump = by_lead[lead]
            member_rows: list[tuple[str, str]] = []
            for member in clump["members"]:
                found = connection.execute(
                    "SELECT evidence_id,p FROM eligible WHERE pair_id=? AND snp=?", (pair_id, member),
                ).fetchall()
                if len(found) != 1:
                    raise CollationError("PLINK clump member lacks exact eligible evidence identity")
                member_rows.append((str(found[0][0]), str(found[0][1])))
            lead_records.append({
                "pair_id": pair_id, "source": source, "block": block,
                "clump_members": [value[0] for value in member_rows],
                "pairwise": result["pairwise"],
            })
    block_leads: dict[tuple[str, str], list[str]] = defaultdict(list)
    for record in lead_records:
        block = record["block"]
        block_id = f"LOC{block.locus}" if block is not None else "UNMAPPED"
        block_leads[(record["pair_id"], block_id)].append(record["source"]["evidence_id"])
    primary: list[dict[str, str]] = []
    control: list[dict[str, str]] = []
    pair_rank = {pair: index for index, pair in enumerate(PAIR_ORDER)}
    for record in sorted(
        lead_records,
        key=lambda item: (
            pair_rank[item["pair_id"]], int(item["source"]["chr"]),
            int(item["source"]["bp"]), item["source"]["snp"],
        ),
    ):
        pair_id = record["pair_id"]
        trait1, trait2, role = PAIR_IDENTITIES[pair_id]
        source, block = record["source"], record["block"]
        chromosome, position = int(source["chr"]), int(source["bp"])
        chromosome_stop = int(chromosome_lengths[str(chromosome)])
        block_id = f"LOC{block.locus}" if block is not None else "UNMAPPED"
        if block is None:
            lava_accounting = "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;OFFICIAL_BLOCK=UNMAPPED"
            conditional = "BLOCKED_BY_OFFICIAL_LAVA_PARTITION_COVERAGE"
            block_start = block_stop = "NA"
        else:
            diagnostic = lava_index.get((pair_id, str(block.locus)))
            if diagnostic is None:
                raise CollationError("independent lead block lacks complete LAVA diagnostic accounting")
            lava_accounting = (
                "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;"
                f"LOCUS_STATUS={diagnostic['locus_status']};"
                f"BIVARIATE_STATUS={diagnostic['bivariate_status']};SCIENTIFIC_SUPPORT=FORBIDDEN"
            )
            conditional = (
                "TERMINAL_FAILED_QC_ACCOUNTING_ONLY;FAILURE_STATUSES="
                + diagnostic["conditional_failures"]
            )
            block_start, block_stop = str(block.start), str(block.stop)
        maximum_r2 = record["pairwise"]["maximum_r2"]
        row = {
            "pair_id": pair_id, "lead_variant": source["snp"],
            "chr": str(chromosome), "position": str(position),
            "PLACO_P": source["p"], "FDR": source["q"],
            "trait1_P": source["p1"], "trait2_P": source["p2"],
            "locus_start": str(max(1, position - 1_000_000)),
            "locus_end": str(min(chromosome_stop, position + 1_000_000)),
            "independent_signal": "TRUE",
            "annotations": (
                "STATISTICAL_CROSS_TRAIT_PLEIOTROPY_ONLY;"
                "NO_SHARED_CAUSAL_VARIANT_NO_MEDIATION_NO_CAUSAL_DIRECTION"
            ),
            "analysis_id": "track-b-v1.0-pleiotropy", "family_role": role,
            "trait1": trait1, "trait2": trait2, "lead_A1": source["a1"],
            "lead_A2": source["a2"], "lead_Z1": source["z1"], "lead_Z2": source["z2"],
            "evidence_labels": source["evidence_labels"],
            "primary_headline": (
                "NOT_APPLICABLE_CONTROL" if pair_id == "CONTROL"
                else str(float(source["p"]) <= 2.5e-8).upper()
            ),
            "eligibility_basis": source["eligibility_basis"],
            "clump_id": f"PLACO:{pair_id}:{source['snp']}",
            "clump_evidence_ids": ";".join(record["clump_members"]),
            "lava_block_id": block_id, "lava_block_start": block_start,
            "lava_block_end": block_stop,
            "pair_block_evidence_ids": ";".join(block_leads[(pair_id, block_id)]),
            "ld_reference_status": "EXACT_RSID_CHR_BP_UNORDERED_NONPALINDROMIC_ALLELES_MATCHED",
            "plink_clump_status": "INDEPENDENT_LEAD_REAL_PLINK_1_9",
            "pairwise_r2_status": record["pairwise"]["status"],
            "max_pairwise_r2_within_1mb": (
                "NA" if maximum_r2 is None else format(float(maximum_r2), ".17g")
            ),
            "lava_terminal_state": "TERMINAL_FAILED_QC_DIAGNOSTIC_ONLY",
            "lava_failed_qc_accounting": lava_accounting,
            "conditional_status": conditional,
            "method_availability": f"PLACO:COMPLETE_WITH_HITS;CONJFDR:{CONJFDR_BLOCKED}",
            "claim_limit": (
                "CONTROL_RECOVERY_ONLY_NO_NOVELTY_NO_MECHANISM"
                if pair_id == "CONTROL" else
                "STATISTICAL_PLEIOTROPY_ONLY_NO_SHARED_CAUSAL_VARIANT_NO_MECHANISM_OR_MEDIATION"
            ),
        }
        (control if pair_id == "CONTROL" else primary).append(row)
    return primary, control


def comparison_rows(shared: Iterable[Mapping[str, str]]) -> Iterator[dict[str, str]]:
    for row in shared:
        yield {
            "analysis_id": "track-b-v1.0-pleiotropy", "pair_id": row["pair_id"],
            "family_role": row["family_role"], "union_locus_id": row["clump_id"],
            "CHR": row["chr"], "START": row["locus_start"], "STOP": row["locus_end"],
            "placo_lead_snp": row["lead_variant"], "placo_P": row["PLACO_P"],
            "placo_BH_Q": row["FDR"], "conjfdr_lead_snp": "NA", "conjfdr": "NA",
            "evidence_labels": row["evidence_labels"],
            "method_availability": row["method_availability"],
            "comparison_label": "NOT_EMITTED_OTHER_METHOD_NOT_TESTED",
            "claim_limit": row["claim_limit"],
        }


def package_source_map(
    science_state: str, *, primary_has_loci: bool = False, control_has_loci: bool = False,
) -> dict[Path, Path]:
    mapping = {
        PLACO_07_PATH: Path("science/07_placo_plus_variants.tsv"),
        LD_COVERAGE_PATH: Path("qc/placo_ld_coverage.tsv"),
        CONTROL_COVERAGE_PATH: Path("control/ld_reference_coverage.tsv"),
        RAM_COMPONENT_PATH: Path("resource/RAM_BENCHMARK.tsv"),
        RAM_NAMESPACE_PATH: Path("resource/RAM_BENCHMARK.namespace.json"),
        RAM_REPORT_PATH: Path("resource/RAM_AWARE_EXECUTION_REPORT.md"),
        RESULT_PROVENANCE_PATH: Path("provenance/results.provenance.json"),
    }
    if science_state == "ZERO_ELIGIBLE_SIGNAL_FAMILY":
        mapping[ZERO_FAMILY_PATH] = Path("provenance/zero_family.provenance.json")
    elif science_state == "BLOCKED_BY_LD_REFERENCE_COVERAGE":
        mapping[BLOCKED_COVERAGE_PATH] = Path("provenance/ld_coverage_blocked_family.provenance.json")
    elif science_state == "COMPLETE_WITH_INDEPENDENT_LOCI":
        if primary_has_loci:
            mapping.update({
                SHARED_LOCI_08_PATH: Path("science/08_shared_loci.tsv"),
                COMPARISON_09_PATH: Path("science/09_pleiotropy_comparison.tsv"),
            })
        if control_has_loci:
            mapping.update({
                CONTROL_08_PATH: Path("control/08_shared_loci.tsv"),
                CONTROL_09_PATH: Path("control/09_pleiotropy_comparison.tsv"),
            })
        if not primary_has_loci and not control_has_loci:
            raise CollationError("complete locus state requires at least one primary/control lead")
    else:
        raise CollationError(f"unknown post-PLACO science state: {science_state}")
    return mapping


def collation_fingerprint(
    root: Path, contract: Mapping[str, Any], pairs: Sequence[PairInput],
) -> str:
    return _digest_json({
        "schema_version": SCHEMA,
        "contract": portable_identity(root, COLLATION_CONTRACT),
        "pair_order": list(PAIR_ORDER),
        "pairs": {
            pair.pair_id: {
                "ledger": dict(pair.ledger_identity),
                "provenance": dict(pair.provenance_identity),
                "terminal_state": pair.terminal_state,
                "rows": pair.expected_rows,
            }
            for pair in pairs
        },
        "ld_reference": contract["ld_reference"],
        "official_lava_partition": contract["official_lava_partition"],
        "collator": portable_identity(root, SCRIPT),
        "scientific_rules": {
            "eligible": "P_LE_5E_8_OR_BH_Q_LE_0_05",
            "headline_A_B": "P_LE_2_5E_8_LABEL_ONLY",
            "clumping": "REAL_PLINK_1_9_1MB_R2_0_1",
            "verification": "ALL_WITHIN_1MB_LEAD_PAIRS_R2_STRICTLY_LT_0_1",
            "control": "SEPARATE_NON_NOVELTY",
            "conjfdr": CONJFDR_BLOCKED,
        },
    })


def _package_manifest(package: Path) -> dict[str, Any]:
    manifest_path = package / "package.manifest.json"
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise CollationError("complete collation package lacks a real manifest")
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CollationError(f"collation package manifest is unreadable: {error}") from error
    if not isinstance(payload, dict):
        raise CollationError("collation package manifest is not an object")
    return payload


def verify_package(root: Path, fingerprint: str) -> dict[str, Any]:
    if not SHA256_RE.fullmatch(fingerprint):
        raise CollationError("invalid collation package fingerprint")
    package_relative = PACKAGE_ROOT / fingerprint
    package = _safe_path(root, package_relative, "collation package")
    observed = os.lstat(package)
    if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
        raise CollationError("collation package is not a real directory")
    manifest = read_json(root, package_relative / "package.manifest.json")
    if (
        manifest.get("schema_version") != "sleep-atlas-track-b-post-placo-package.2"
        or manifest.get("fingerprint") != fingerprint
        or not isinstance(manifest.get("files"), list)
    ):
        raise CollationError("collation package manifest identity drifted")
    expected_paths: set[str] = set()
    for record in manifest["files"]:
        if (
            not isinstance(record, Mapping) or not isinstance(record.get("path"), str)
            or record["path"] in expected_paths or Path(record["path"]).is_absolute()
            or ".." in Path(record["path"]).parts
        ):
            raise CollationError("collation package manifest has an unsafe/duplicate path")
        expected_paths.add(record["path"])
        path = package / record["path"]
        if path.is_symlink() or not path.is_file():
            raise CollationError(f"collation package artifact is absent/special: {record['path']}")
        identity = _stage_identity(path, Path(record["path"]), allow_empty=record.get("bytes") == 0)
        if any(identity.get(key) != record.get(key) for key in ("path", "bytes", "sha256")):
            raise CollationError(f"collation package artifact drifted: {record['path']}")
    observed_paths: set[str] = set()
    for path in package.rglob("*"):
        status = os.lstat(path)
        if stat.S_ISDIR(status.st_mode):
            if stat.S_ISLNK(status.st_mode):
                raise CollationError("collation package contains a directory symlink")
            continue
        if not stat.S_ISREG(status.st_mode):
            raise CollationError(f"collation package contains a special file: {path}")
        relative = str(path.relative_to(package))
        if relative != "package.manifest.json":
            observed_paths.add(relative)
    if observed_paths != expected_paths:
        raise CollationError("collation package file family differs from its manifest")
    return manifest


def _rename_directory_no_replace(source: Path, destination: Path) -> None:
    """Kernel no-replace directory commit on macOS/Linux."""
    import ctypes
    import ctypes.util

    try:
        libc = ctypes.CDLL(ctypes.util.find_library("c") or None, use_errno=True)
    except OSError as error:
        raise CollationError(f"could not load libc for package publication: {error}") from error
    source_raw, destination_raw = os.fsencode(source), os.fsencode(destination)
    if sys.platform == "darwin":
        try:
            rename = libc.renamex_np
        except AttributeError as error:
            raise CollationError("renamex_np unavailable for atomic no-replace package") from error
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
        raise CollationError("kernel atomic no-replace directory publication is unavailable")
    if result:
        number = ctypes.get_errno()
        if number in {errno.EEXIST, errno.ENOTEMPTY}:
            raise CollationError("collation package destination appeared concurrently")
        raise CollationError(f"collation package rename failed: {os.strerror(number)}")
    _fsync_directory(destination.parent)


def _peak_rss_bytes(sampled_parent_plus_child_peak: int = 0) -> int:
    parent_only_peak = normalize_ru_maxrss(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return max(parent_only_peak, int(sampled_parent_plus_child_peak))


def _ram_rows(total_rows: int, peak_bytes: int, elapsed: float, output_hash: str) -> list[dict[str, str]]:
    return [{
        "analysis": "TRACK_B_PLEIOTROPY_POST_PLACO_COLLATION_V2",
        "pair": "A;B;CONTROL", "locus": "GENOME_WIDE_COMPLETE_VALID_VARIANT_FAMILY",
        "chromosome": "ALL", "n_snps": str(total_rows),
        "peak_ram_gb": format(peak_bytes / 1024**3, ".12g"),
        "runtime_sec": format(elapsed, ".12g"), "exit_status": "0",
        "output_hash": output_hash,
    }]


def _ram_report(
    *, total_rows: int, peak_bytes: int, elapsed: float, feasible: bool,
    pair_results: Mapping[str, Mapping[str, Any]], measurement_backend: str,
) -> str:
    measured = peak_bytes / 1024**3
    plink_runtime = sum(
        float(command["runtime_seconds"])
        for result in pair_results.values() for command in result["commands"]
    )
    command_peak_lines = [
        (
            f"- {pair_id} PLINK {index}: live parent+PGID peak "
            f"{int(command['peak_aggregate_rss_bytes']) / 1024**3:.6f} GiB; "
            f"runtime {float(command['runtime_seconds']):.6f} seconds; "
            f"backend {command['rss_measurement_method']}"
        )
        for pair_id, result in pair_results.items()
        for index, command in enumerate(result["commands"], start=1)
    ] or ["- No PLINK child was required because no eligible exact-reference signal existed."]
    return "\n".join([
        "# Track B post-PLACO RAM-aware execution report",
        "",
        "This report covers only the versioned post-PLACO collation component. It does not append to or overwrite the legacy top-level RAM namespace.",
        "",
        "## Safe decomposition",
        "",
        "The three complete PLACO ledgers were streamed sequentially. PLINK clumping and pairwise-r2 verification ran one frozen pair at a time. No family was subset for memory or runtime.",
        "",
        "## Non-decomposable scientific stages",
        "",
        "The upstream PLACO global nuisance variance/correlation fit and complete within-pair BH denominator remain indivisible all-valid-variant operations. This stage only consumes their sealed outputs. The sealed BIM was scanned as one complete reference identity, and each pair's full independent-lead family was jointly audited within 1 Mb.",
        "",
        "## Measurements and 8-GB feasibility",
        "",
        f"- Complete 07 rows: {total_rows}",
        f"- Measured peak process/child RSS: {measured:.6f} GiB",
        f"- Authoritative RSS measurement backend: {measurement_backend}",
        f"- End-to-end component runtime: {elapsed:.6f} seconds",
        f"- PLINK subprocess runtime: {plink_runtime:.6f} seconds",
        f"- Peak plus 1-GiB reserve fits 8 GiB: {'YES' if feasible else 'NO'}",
        "",
        "### Per-command measured peaks",
        "",
        *command_peak_lines,
        "",
        "## Larger-server stages",
        "",
        "The independently blocked MATLAB conjFDR family still requires its policy minimum (16 GiB RAM and 20 GiB free storage). Any future cross-method r2>=0.6 reconciliation must be separately sealed; this component does not infer it. Upstream PLACO global nuisance fitting is not re-estimated here.",
        "",
        "The adjacent component benchmark uses the exact federatable columns `analysis, pair, locus, chromosome, n_snps, peak_ram_gb, runtime_sec, exit_status, output_hash`.",
        "",
    ])


def _normalize_command_records(
    attempt: Path, records: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    normalized: dict[str, Any] = {}
    for pair_id, result in records.items():
        commands = []
        for command in result["commands"]:
            transcript = Path(command["transcript"])
            commands.append({
                "argv": command["argv"], "runtime_seconds": command["runtime_seconds"],
                "returncode": command["returncode"],
                "peak_aggregate_rss_bytes": command["peak_aggregate_rss_bytes"],
                "rss_measurement_method": command["rss_measurement_method"],
                "transcript": _stage_identity(
                    transcript, transcript.relative_to(attempt), allow_empty=False,
                ),
            })
        normalized[pair_id] = {
            "matched_rows": result["matched_rows"],
            "missing_from_clump": result.get("missing_from_clump", []),
            "lead_count": len(result["leads"]), "pairwise": result["pairwise"],
            "commands": commands,
        }
    return normalized


def _stage_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    _write_stage_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _freeze_runtime_file(
    root: Path, expected: Mapping[str, Any], destination: Path,
) -> dict[str, Any]:
    """Hard-link one sealed runtime input so named-path replacement cannot reach PLINK."""
    expected_path = expected.get("path")
    if not isinstance(expected_path, str):
        raise CollationError("runtime identity lacks a repository-relative path")
    validated_identity = stable_identity(root, expected_path)
    if any(validated_identity.get(key) != expected.get(key) for key in ("path", "bytes", "sha256")):
        raise CollationError(f"runtime input differs from pre-result contract: {expected_path}")
    if destination.exists() or destination.is_symlink():
        raise CollationError(f"runtime snapshot destination already exists: {destination}")
    try:
        os.link(root / expected_path, destination, follow_symlinks=False)
    except OSError as error:
        raise CollationError(f"could not freeze runtime input {expected_path}: {error}") from error
    linked = os.lstat(destination)
    if (
        not stat.S_ISREG(linked.st_mode)
        or (linked.st_dev, linked.st_ino)
        != (validated_identity["device"], validated_identity["inode"])
    ):
        raise CollationError("runtime snapshot is not the exact sealed source inode")
    _fsync_directory(destination.parent)
    # Creating a hard link legitimately changes inode ctime, so freeze the
    # cross-command baseline only after the link exists.
    source_identity = stable_identity(root, expected_path)
    snapshot_identity = stable_identity(root, destination.relative_to(root))
    if (source_identity["device"], source_identity["inode"]) != (
        snapshot_identity["device"], snapshot_identity["inode"],
    ):
        raise CollationError("runtime snapshot stopped naming the sealed source inode")
    return {"source": source_identity, "snapshot": snapshot_identity}


def freeze_runtime_family(
    root: Path, contract: Mapping[str, Any], attempt: Path,
) -> dict[str, Any]:
    reference_directory = attempt / "runtime_reference"
    software_directory = attempt / "runtime_software"
    reference_directory.mkdir()
    software_directory.mkdir()
    reference_records: dict[str, Any] = {}
    for suffix in ("bed", "bim", "fam"):
        reference_records[suffix] = _freeze_runtime_file(
            root, contract["ld_reference"][suffix],
            reference_directory / f"g1000_eur.{suffix}",
        )
    plink_expected = contract["ld_reference"]["plink"]
    plink_record = _freeze_runtime_file(root, plink_expected, software_directory / "plink")
    plink_path = software_directory / "plink"
    if not os.access(plink_path, os.X_OK):
        raise CollationError("frozen pinned PLINK inode is not executable")
    result = subprocess.run(
        [str(plink_path), "--version"], check=False, text=True, capture_output=True,
    )
    version = (result.stdout + result.stderr).strip()
    if result.returncode or not version.startswith(str(plink_expected.get("version", ""))):
        raise CollationError(f"frozen pinned PLINK version drifted: {version!r}")
    return {
        "reference": reference_records, "plink": plink_record,
        "reference_prefix": reference_directory / "g1000_eur",
        "plink_path": plink_path, "plink_version": version,
    }


def verify_and_release_runtime_family(root: Path, frozen: Mapping[str, Any]) -> dict[str, Any]:
    """Rehash the exact inodes PLINK used, then remove only our private hard links."""
    audited: dict[str, Any] = {"reference": {}}
    records = [
        *(frozen["reference"][suffix] for suffix in ("bed", "bim", "fam")),
        frozen["plink"],
    ]
    for record in records:
        before_source, before_snapshot = record["source"], record["snapshot"]
        source = stable_identity(root, before_source["path"])
        snapshot = stable_identity(root, before_snapshot["path"])
        invariant_fields = (
            "device", "inode", "bytes", "mtime_ns", "ctime_ns", "sha256",
        )
        if any(source.get(key) != before_source.get(key) for key in invariant_fields):
            raise CollationError("sealed runtime source changed across PLINK execution")
        if any(snapshot.get(key) != before_snapshot.get(key) for key in invariant_fields):
            raise CollationError("frozen runtime inode changed across PLINK execution")
        if (source["device"], source["inode"]) != (snapshot["device"], snapshot["inode"]):
            raise CollationError("runtime source path stopped naming PLINK's frozen inode")
        portable = {key: source[key] for key in ("path", "bytes", "sha256")}
        if str(source["path"]).endswith((".bed", ".bim", ".fam")):
            audited["reference"][str(source["path"])[-3:]] = portable
        else:
            audited["plink"] = portable
    audited["snapshot_mode"] = (
        "HARDLINKED_EXACT_INODES_BEFORE_BIM_SCAN_AND_PLINK;FULL_HASH_AND_"
        "INODE_MTIME_CTIME_RECHECK_AFTER_ALL_COMMANDS"
    )
    snapshot_paths = [Path(record["snapshot"]["path"]) for record in records]
    for relative in snapshot_paths:
        path = root / relative
        current = os.lstat(path)
        if (current.st_dev, current.st_ino) != (
            next(item for item in records if item["snapshot"]["path"] == str(relative))["snapshot"]["device"],
            next(item for item in records if item["snapshot"]["path"] == str(relative))["snapshot"]["inode"],
        ):
            raise CollationError("runtime snapshot inode changed before private-link release")
        path.unlink()
        _fsync_directory(path.parent)
    for directory in (Path(frozen["reference_prefix"]).parent, Path(frozen["plink_path"]).parent):
        directory.rmdir()
        _fsync_directory(directory.parent)
    return audited


def _inventory_files(base: Path, *, exclude: Iterable[Path] = ()) -> list[dict[str, object]]:
    excluded = {str(path) for path in exclude}
    records: list[dict[str, object]] = []
    for path in sorted(base.rglob("*"), key=lambda item: str(item.relative_to(base))):
        observed = os.lstat(path)
        if stat.S_ISDIR(observed.st_mode):
            if stat.S_ISLNK(observed.st_mode):
                raise CollationError("private package contains a directory symlink")
            continue
        relative = path.relative_to(base)
        if str(relative) in excluded:
            continue
        if not stat.S_ISREG(observed.st_mode):
            raise CollationError(f"private package contains a special file: {relative}")
        records.append(_stage_identity(path, relative, allow_empty=observed.st_size == 0))
    return records


def _inventory_subtree(package: Path, subtree: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    base = package / subtree
    for path in sorted(base.rglob("*"), key=lambda item: str(item.relative_to(package))):
        observed = os.lstat(path)
        if stat.S_ISDIR(observed.st_mode):
            if stat.S_ISLNK(observed.st_mode):
                raise CollationError("private package subtree contains a directory symlink")
            continue
        if not stat.S_ISREG(observed.st_mode):
            raise CollationError("private package subtree contains a special file")
        records.append(_stage_identity(
            path, path.relative_to(package), allow_empty=observed.st_size == 0,
        ))
    return records


def _output_records(
    attempt: Path, mapping: Mapping[Path, Path], *, omit: Iterable[Path] = (),
) -> list[dict[str, Any]]:
    omitted = set(omit)
    records: list[dict[str, Any]] = []
    for canonical, packaged in mapping.items():
        if canonical in omitted:
            continue
        identity = _stage_identity(
            attempt / packaged, packaged,
            allow_empty=(attempt / packaged).stat().st_size == 0,
        )
        records.append({"canonical_path": str(canonical), **identity})
    return records


def _pair_input_payload(pairs: Sequence[PairInput]) -> list[dict[str, Any]]:
    return [
        {
            "pair_id": pair.pair_id, "terminal_state": pair.terminal_state,
            "expected_rows": pair.expected_rows,
            "ledger": dict(pair.ledger_identity),
            "provenance": dict(pair.provenance_identity),
        }
        for pair in pairs
    ]


def verify_pair_input_identities(root: Path, pairs: Sequence[PairInput]) -> None:
    for pair in pairs:
        if (
            not identity_matches(root, pair.ledger_identity)
            or not identity_matches(root, pair.provenance_identity)
        ):
            raise CollationError(
                f"authoritative PLACO {pair.pair_id} input changed during post-processing"
            )


def _existing_owned_outputs(root: Path) -> list[Path]:
    return [
        path for path in sorted(CANONICAL_SCIENCE_PATHS, key=str)
        if _artifact_present(root / path)
    ]


def _build_package_from_verified_inputs(
    root: Path,
    contract: Mapping[str, Any],
    pairs: Sequence[PairInput],
    fingerprint: str,
    archive_state: Mapping[str, Any],
    conjfdr: Mapping[str, Any],
    lava_index: Mapping[tuple[str, str], Mapping[str, str]],
    *,
    runner: Callable[[Sequence[str], Path], subprocess.CompletedProcess[str]] = _default_runner,
) -> dict[str, Any]:
    """Build one immutable package after every upstream family has been verified."""
    if [pair.pair_id for pair in pairs] != list(PAIR_ORDER):
        raise CollationError("package inputs are not the exact A/B/CONTROL family")
    if archive_state_invariant(root, archive_state) != contract["archive_reference_state_guard"]:
        raise CollationError("live LAVA archive/reference identity differs from the sealed invariant")
    if conjfdr.get("family_state") != CONJFDR_BLOCKED:
        raise CollationError("unsupported conjFDR state reached blocked-method collation")
    destination = root / PACKAGE_ROOT / fingerprint
    if _artifact_present(destination):
        manifest = verify_package(root, fingerprint)
        return {"fingerprint": fingerprint, "manifest": manifest, "reused": True}
    offenders = _existing_owned_outputs(root)
    if offenders:
        raise CollationError(
            "canonical result artifact exists without its immutable package: "
            + ", ".join(str(path) for path in offenders[:8])
        )
    competing = [path for path in sorted(COMPETING_PATHS, key=str) if _artifact_present(root / path)]
    if competing:
        raise CollationError(
            "fine-mapping/competing evidence exists before post-PLACO publication: "
            + ", ".join(str(path) for path in competing)
        )

    package_parent = _ensure_directories(root, PACKAGE_ROOT)
    staging_parent = _ensure_directories(root, PACKAGE_STAGING_ROOT)
    attempt = Path(tempfile.mkdtemp(prefix=f"{fingerprint}.", dir=staging_parent))
    _fsync_directory(staging_parent)
    for relative in (
        Path("science"), Path("qc"), Path("control"), Path("resource"),
        Path("provenance"), Path("plink"),
    ):
        (attempt / relative).mkdir()
    _fsync_directory(attempt)
    started = time.monotonic()
    database = attempt / "eligible.sqlite"
    connection = create_index_database(database)
    completed = False
    try:
        family_counts, projected_07 = index_complete_placo_family(root, pairs, connection)
        eligible_rows = int(connection.execute("SELECT COUNT(*) FROM eligible").fetchone()[0])
        storage = disk_preflight(root, projected_07, eligible_rows)
        complete_07 = attempt / "science/07_placo_plus_variants.tsv"
        complete_07_identity = write_complete_07(
            root, pairs, complete_07, projected_07,
        )
        complete_07_identity["path"] = "science/07_placo_plus_variants.tsv"

        frozen = freeze_runtime_family(root, contract, attempt)
        snapshot_bim_relative = Path(frozen["reference"]["bim"]["snapshot"]["path"])
        reference_coverage = match_eligible_to_bim(
            root, connection, snapshot_bim_relative,
            {
                key: frozen["reference"]["bim"]["snapshot"][key]
                for key in ("path", "bytes", "sha256")
            },
        )
        if reference_coverage["bim_rows_streamed"] != contract["ld_reference"]["variant_count"]:
            raise CollationError("sealed BIM row count differs from the pinned full reference")
        plink_results: dict[str, dict[str, Any]] = {}
        for pair_id in PAIR_ORDER:
            plink_results[pair_id] = run_plink_clump_for_pair(
                connection, pair_id, Path(frozen["plink_path"]),
                Path(frozen["reference_prefix"]), attempt / "plink" / pair_id,
                root=root, runner=runner,
            )
        runtime_audit = verify_and_release_runtime_family(root, frozen)

        _, block_index = load_official_blocks(root, contract["official_lava_partition"])
        primary_rows, control_rows = build_shared_locus_rows(
            connection, plink_results, block_index,
            contract["ld_reference"]["chromosome_lengths"], lava_index,
        )
        total_leads = len(primary_rows) + len(control_rows)
        if eligible_rows == 0:
            science_state = "ZERO_ELIGIBLE_SIGNAL_FAMILY"
        elif total_leads == 0:
            science_state = "BLOCKED_BY_LD_REFERENCE_COVERAGE"
        else:
            science_state = "COMPLETE_WITH_INDEPENDENT_LOCI"
        primary_has_loci, control_has_loci = bool(primary_rows), bool(control_rows)
        mapping = package_source_map(
            science_state, primary_has_loci=primary_has_loci,
            control_has_loci=control_has_loci,
        )

        _stage_write_tsv(
            attempt / "qc/placo_ld_coverage.tsv", LD_COVERAGE_FIELDS,
            _coverage_rows(connection),
        )
        _stage_write_tsv(
            attempt / "control/ld_reference_coverage.tsv", LD_COVERAGE_FIELDS,
            _coverage_rows(connection, pair_filter="CONTROL"),
        )
        if primary_has_loci:
            _stage_write_tsv(
                attempt / "science/08_shared_loci.tsv", SHARED_LOCI_08_FIELDS, primary_rows,
            )
            _stage_write_tsv(
                attempt / "science/09_pleiotropy_comparison.tsv", COMPARISON_09_FIELDS,
                comparison_rows(primary_rows),
            )
        if control_has_loci:
            _stage_write_tsv(
                attempt / "control/08_shared_loci.tsv", SHARED_LOCI_08_FIELDS, control_rows,
            )
            _stage_write_tsv(
                attempt / "control/09_pleiotropy_comparison.tsv", COMPARISON_09_FIELDS,
                comparison_rows(control_rows),
            )

        coverage_failures = eligible_rows - int(reference_coverage.get("MATCHED", 0))
        terminal_receipt: dict[str, Any] | None = None
        if science_state == "ZERO_ELIGIBLE_SIGNAL_FAMILY":
            terminal_receipt = {
                "schema_version": ZERO_SCHEMA, "analysis_id": "track-b-v1.0-pleiotropy",
                "fingerprint": fingerprint, "state": science_state,
                "complete_family_rows": sum(pair.expected_rows for pair in pairs),
                "eligible_signal_rows": 0, "independent_loci": 0,
                "pair_counts": family_counts,
                "interpretation": (
                    "SEALED_COMPLETE_PLACO_FAMILY_HAS_NO_POLICY_ELIGIBLE_SIGNAL;"
                    "NOT_A_CONJFDR_RESULT_AND_NOT_A_CROSS_METHOD_NULL"
                ),
            }
            _stage_write_json(attempt / "provenance/zero_family.provenance.json", terminal_receipt)
        elif science_state == "BLOCKED_BY_LD_REFERENCE_COVERAGE":
            terminal_receipt = {
                "schema_version": "sleep-atlas-track-b-post-placo-ld-coverage-block.2",
                "analysis_id": "track-b-v1.0-pleiotropy", "fingerprint": fingerprint,
                "state": science_state, "eligible_signal_rows": eligible_rows,
                "exact_reference_matched_rows": 0, "coverage_failure_rows": coverage_failures,
                "interpretation": (
                    "ELIGIBLE_PLACO_SIGNALS_RETAINED_IN_EXPLICIT_COVERAGE_LEDGER;"
                    "NO_CANONICAL_LOCUS_OR_NULL_CLAIM"
                ),
            }
            _stage_write_json(
                attempt / "provenance/ld_coverage_blocked_family.provenance.json",
                terminal_receipt,
            )

        elapsed = time.monotonic() - started
        commands_flat = [
            command for result in plink_results.values() for command in result["commands"]
        ]
        if runner is _default_runner and any(
            command["peak_aggregate_rss_bytes"] <= 0
            or command["rss_measurement_method"] not in {
                "LIVE_PARENT_PLUS_OWNED_FRESH_SESSION_PGID_PS_SAMPLED",
                "NORMALIZED_SELF_PLUS_CHILD_RUSAGE_CONSERVATIVE_FAST_EXIT_FALLBACK",
            }
            for command in commands_flat
        ):
            raise CollationError("production PLINK command lacks authoritative RSS measurement")
        sampled_peak = max(
            (int(command["peak_aggregate_rss_bytes"]) for command in commands_flat),
            default=0,
        )
        parent_only_peak = normalize_ru_maxrss(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        )
        peak = max(parent_only_peak, sampled_peak)
        measurement_backend = (
            "MAX_OF_NORMALIZED_PARENT_RUSAGE_AND_LIVE_PARENT_PLUS_OWNED_FRESH_"
            "SESSION_PGID_PS_SAMPLES;FAST_EXIT_USES_CONSERVATIVE_SUM_OF_NORMALIZED_MAXIMA"
            if commands_flat else "NORMALIZED_PARENT_RUSAGE_NO_PLINK_CHILD_REQUIRED"
        )
        feasible = peak + 1024**3 <= 8 * 1024**3
        if not feasible:
            raise CollationError(
                "BLOCKED_BY_COMPUTE measured post-PLACO peak plus 1-GiB reserve exceeds 8 GiB"
            )
        _stage_write_tsv(
            attempt / "resource/RAM_BENCHMARK.tsv", RAM_FIELDS,
            _ram_rows(
                sum(pair.expected_rows for pair in pairs), peak, elapsed,
                str(complete_07_identity["sha256"]),
            ),
        )
        ram_identity = _stage_identity(
            attempt / "resource/RAM_BENCHMARK.tsv", Path("resource/RAM_BENCHMARK.tsv"),
        )
        _stage_write_json(
            attempt / "resource/RAM_BENCHMARK.namespace.json",
            {
                "schema_version": "track-b-ram-benchmark-component-namespace.1",
                "component": "post_placo_collation_v2",
                "component_benchmark": ram_identity,
                "legacy_top_level_mutated": False,
                "federation_key": ["analysis", "pair", "locus", "chromosome"],
                "collision_policy": "FAIL_ON_DUPLICATE_FEDERATION_KEY_OR_DIFFERING_ROW_BYTES",
                "exact_schema": RAM_FIELDS,
            },
        )
        _write_stage_text(
            attempt / "resource/RAM_AWARE_EXECUTION_REPORT.md",
            _ram_report(
                total_rows=sum(pair.expected_rows for pair in pairs), peak_bytes=peak,
                elapsed=elapsed, feasible=feasible, pair_results=plink_results,
                measurement_backend=measurement_backend,
            ),
        )

        normalized_commands = _normalize_command_records(attempt, plink_results)
        verify_pair_input_identities(root, pairs)
        plink_artifacts = _inventory_subtree(attempt, Path("plink"))
        outputs_without_provenance = _output_records(
            attempt, mapping, omit={RESULT_PROVENANCE_PATH},
        )
        reference_status_counts = {
            str(status): int(count)
            for status, count in connection.execute(
                "SELECT reference_status,COUNT(*) FROM eligible GROUP BY reference_status"
            )
        }
        inputs = _pair_input_payload(pairs)
        provenance = {
            "schema_version": SCHEMA, "analysis_id": "track-b-v1.0-pleiotropy",
            "fingerprint": fingerprint, "science_state": science_state,
            "collation_contract": portable_identity(root, COLLATION_CONTRACT),
            "policy_sha256": contract["upstream"][str(POLICY)]["sha256"],
            "contract_lock_sha256": contract["upstream"][str(PLEIOTROPY_CONTRACT_LOCK)]["sha256"],
            "input_gate_lock_sha256": contract["upstream"][str(INPUT_GATE_LOCK)]["sha256"],
            "software_sha256_or_commit": contract["upstream"][str(SCRIPT)]["sha256"],
            "reference_sha256": contract["ld_reference"]["bim"]["sha256"],
            "input_sha256": _digest_json(inputs),
            "input_rows": sum(pair.expected_rows for pair in pairs),
            "output_sha256": complete_07_identity["sha256"],
            "output_rows": sum(pair.expected_rows for pair in pairs),
            "exact_schema": ",".join(PLACO_07_FIELDS),
            "complete_family_counts": family_counts, "qc_status": "PASS",
            "inputs": inputs,
            "complete_143_family_verification": {
                "verifier": contract["upstream"][str(SEQUENTIAL_SCRIPT)],
                "outcomes": [
                    {"pair": pair, "outcome": "VERIFIED_COMPLETE_CLEAN"}
                    for pair in PAIR_ORDER
                ],
                "performed_after_pre_result_contract": True,
                "all_valid_variants_and_global_bh_denominators_deeply_validated": True,
            },
            "complete_family_order": list(PAIR_ORDER),
            "complete_family_role": {
                "primary": ["A", "B"], "control": ["CONTROL"],
            },
            "eligible_signal_rows": eligible_rows,
            "independent_primary_loci": len(primary_rows),
            "independent_control_loci": len(control_rows),
            "reference_coverage": {
                **reference_coverage, "status_counts": reference_status_counts,
                "failure_rows": coverage_failures,
            },
            "thresholds": contract["thresholds"],
            "ld_reference": contract["ld_reference"],
            "runtime_reference_audit": runtime_audit,
            "plink_execution": normalized_commands,
            "plink_artifacts": plink_artifacts,
            "official_lava_partition": contract["official_lava_partition"],
            "lava_failed_qc_accounting": contract["lava_failed_qc_accounting"],
            "lava_use_limit": "FAILED_QC_DIAGNOSTIC_ACCOUNTING_ONLY_NOT_SCIENTIFIC_SUPPORT",
            "archive_reference_state": dict(archive_state),
            "archive_reference_state_guard": contract["archive_reference_state_guard"],
            "conjfdr": dict(conjfdr),
            "method_comparison_rule": (
                "NO_PLACO_ONLY_CONJFDR_ONLY_OR_PLACO_AND_CONJFDR_LABEL_WHILE_CONJFDR_BLOCKED"
            ),
            "primary_has_loci": primary_has_loci,
            "control_has_loci": control_has_loci,
            "control_separation": (
                "CONTROL_OUTPUTS_UNDER_RESULTS_TRACK_B_CONTROL_INSOMNIA_FRAILTY;"
                "EXCLUDED_FROM_NOVELTY_AND_PRIMARY_HEADLINE_FAMILY"
            ),
            "terminal_receipt": terminal_receipt,
            "storage_preflight": storage,
            "resource_metrics": {
                "peak_ram_bytes": peak, "runtime_seconds": elapsed,
                "eight_gib_sequential_feasible_with_one_gib_reserve": feasible,
                "benchmark_schema": RAM_FIELDS,
                "rss_measurement_backend": measurement_backend,
                "maximum_live_parent_plus_plink_pgid_rss_bytes": sampled_peak,
                "normalized_parent_only_peak_rss_bytes": parent_only_peak,
            },
            "outputs": outputs_without_provenance,
            "publication": {
                "mode": "IMMUTABLE_PACKAGE_THEN_HARDLINK_NO_REPLACE",
                "result_provenance_linked_last": True,
                "partial_exact_prefix_restartable": True,
                "canonical_overwrite": "FORBIDDEN",
                "competing_finemapping_artifacts_absent_at_commit_required": True,
            },
            "claim_limit": (
                "STATISTICAL_CROSS_TRAIT_PLEIOTROPY_ONLY;NO_SHARED_CAUSAL_VARIANT;"
                "NO_BIOLOGICAL_MECHANISM_MEDIATION_OR_CAUSAL_DIRECTION"
            ),
        }
        _stage_write_json(attempt / "provenance/results.provenance.json", provenance)

        connection.close()
        connection = None
        for candidate in (database, database.with_name(database.name + "-journal")):
            if candidate.exists():
                observed = os.lstat(candidate)
                if not stat.S_ISREG(observed.st_mode):
                    raise CollationError("private SQLite work artifact became non-regular")
                candidate.unlink()
        _fsync_directory(attempt)

        files = _inventory_files(attempt)
        manifest_payload = {
            "schema_version": "sleep-atlas-track-b-post-placo-package.2",
            "fingerprint": fingerprint, "science_state": science_state,
            "primary_has_loci": primary_has_loci,
            "control_has_loci": control_has_loci,
            "canonical_mapping": {
                str(canonical): str(packaged) for canonical, packaged in mapping.items()
            },
            "files": files,
            "publication_order": [str(path) for path in publication_order(mapping)],
        }
        _stage_write_json(attempt / "package.manifest.json", manifest_payload)
        _fsync_directory(attempt)
        _rename_directory_no_replace(attempt, destination)
        completed = True
        manifest = verify_package(root, fingerprint)
        return {"fingerprint": fingerprint, "manifest": manifest, "reused": False}
    finally:
        if connection is not None:
            connection.close()
        # Failed attempts are intentionally retained as non-authoritative forensic state.
        if completed and attempt.exists():
            raise CollationError("package commit reported success but private stage still exists")


def publication_order(mapping: Mapping[Path, Path]) -> list[Path]:
    order = [
        PLACO_07_PATH, LD_COVERAGE_PATH,
        SHARED_LOCI_08_PATH, COMPARISON_09_PATH,
        CONTROL_COVERAGE_PATH, CONTROL_08_PATH, CONTROL_09_PATH,
        RAM_COMPONENT_PATH, RAM_NAMESPACE_PATH, RAM_REPORT_PATH,
        ZERO_FAMILY_PATH, BLOCKED_COVERAGE_PATH, RESULT_PROVENANCE_PATH,
    ]
    selected = [path for path in order if path in mapping]
    if set(selected) != set(mapping) or selected[-1:] != [RESULT_PROVENANCE_PATH]:
        raise CollationError("canonical publication mapping/order is incomplete")
    return selected


def _manifest_mapping(manifest: Mapping[str, Any]) -> dict[Path, Path]:
    raw = manifest.get("canonical_mapping")
    if not isinstance(raw, Mapping):
        raise CollationError("package manifest lacks its canonical mapping")
    mapping: dict[Path, Path] = {}
    for canonical_raw, packaged_raw in raw.items():
        if not isinstance(canonical_raw, str) or not isinstance(packaged_raw, str):
            raise CollationError("package manifest mapping has non-string paths")
        canonical, packaged = Path(canonical_raw), Path(packaged_raw)
        if (
            canonical not in CANONICAL_SCIENCE_PATHS
            or canonical in mapping or packaged.is_absolute() or ".." in packaged.parts
        ):
            raise CollationError("package manifest mapping has an unsafe/duplicate target")
        mapping[canonical] = packaged
    expected = package_source_map(
        str(manifest.get("science_state", "")),
        primary_has_loci=manifest.get("primary_has_loci") is True,
        control_has_loci=manifest.get("control_has_loci") is True,
    )
    if mapping != expected:
        raise CollationError("package manifest canonical mapping differs from its science state")
    if manifest.get("publication_order") != [str(path) for path in publication_order(mapping)]:
        raise CollationError("package manifest publication order drifted")
    return mapping


def _same_regular_inode(left: Path, right: Path) -> bool:
    try:
        first, second = os.lstat(left), os.lstat(right)
    except OSError:
        return False
    return (
        stat.S_ISREG(first.st_mode) and stat.S_ISREG(second.st_mode)
        and (first.st_dev, first.st_ino) == (second.st_dev, second.st_ino)
    )


def publish_package(
    root: Path, fingerprint: str, *, crash_after_links: int | None = None,
) -> dict[str, Any]:
    """Link an immutable package to canonical paths, provenance last, without replacement."""
    manifest = verify_package(root, fingerprint)
    mapping = _manifest_mapping(manifest)
    package = root / PACKAGE_ROOT / fingerprint
    order = publication_order(mapping)
    if crash_after_links is not None and crash_after_links < 1:
        raise CollationError("crash simulation link count must be positive")
    forbidden = [
        path for path in sorted(CANONICAL_SCIENCE_PATHS - set(mapping), key=str)
        if _artifact_present(root / path)
    ]
    if forbidden:
        raise CollationError(
            "science-state-incompatible canonical artifact already exists: "
            + ", ".join(str(path) for path in forbidden)
        )
    competing = [path for path in sorted(COMPETING_PATHS, key=str) if _artifact_present(root / path)]
    if competing:
        raise CollationError("competing fine-mapping evidence blocks result publication")

    gap_seen = False
    for canonical in order:
        source = package / mapping[canonical]
        destination = root / canonical
        present = _artifact_present(destination)
        if present:
            if gap_seen:
                raise CollationError("canonical publication is not an exact restartable prefix")
            if not _same_regular_inode(source, destination):
                raise CollationError(
                    f"canonical no-replace target is not the immutable package inode: {canonical}"
                )
        else:
            gap_seen = True

    linked = 0
    for canonical in order:
        source = package / mapping[canonical]
        destination = root / canonical
        if _artifact_present(destination):
            continue
        if canonical == RESULT_PROVENANCE_PATH:
            late_competing = [
                path for path in sorted(COMPETING_PATHS, key=str)
                if _artifact_present(root / path)
            ]
            if late_competing:
                raise CollationError(
                    "competing fine-mapping evidence appeared before provenance commit"
                )
        parent = _ensure_directories(root, canonical.parent)
        source_before = os.lstat(source)
        if not stat.S_ISREG(source_before.st_mode):
            raise CollationError("package publication source became non-regular")
        try:
            os.link(source, destination, follow_symlinks=False)
        except FileExistsError as error:
            raise CollationError(f"canonical target appeared concurrently: {canonical}") from error
        destination_status = os.lstat(destination)
        if (
            not stat.S_ISREG(destination_status.st_mode)
            or (source_before.st_dev, source_before.st_ino)
            != (destination_status.st_dev, destination_status.st_ino)
        ):
            raise CollationError("canonical hard-link publication changed source identity")
        _fsync_directory(parent)
        linked += 1
        if crash_after_links == linked:
            raise CollationError("SIMULATED_CRASH_AFTER_CANONICAL_LINK")
    for canonical in order:
        if not _same_regular_inode(package / mapping[canonical], root / canonical):
            raise CollationError("canonical package publication failed final inode audit")
    return {
        "fingerprint": fingerprint, "science_state": manifest["science_state"],
        "linked_this_call": linked, "canonical_files": len(mapping),
        "result_provenance_committed_last": True,
    }


def build_and_publish_results(
    root: Path = ROOT,
    *,
    runner: Callable[[Sequence[str], Path], subprocess.CompletedProcess[str]] = _default_runner,
    crash_after_links: int | None = None,
) -> dict[str, Any]:
    """Production route: sealed contract first, then validated ledgers, package, publish."""
    contract = verify_pre_result_contract(root)
    modules = _load_production_modules(root)
    archive_state = verify_archive_state(modules["archive_guard"], root)
    if archive_state_invariant(root, archive_state) != contract["archive_reference_state_guard"]:
        raise CollationError("current archive/reference family differs from the pre-result contract")
    pairs = resolve_production_pair_inputs(root, modules["sequential"])
    fingerprint = collation_fingerprint(root, contract, pairs)
    conjfdr = validate_conjfdr_blocked(root)
    lava_index = load_lava_diagnostic_index(root)
    package = _build_package_from_verified_inputs(
        root, contract, pairs, fingerprint, archive_state, conjfdr, lava_index,
        runner=runner,
    )
    if verify_pre_result_contract(root) != contract:
        raise CollationError("pre-result contract changed during immutable package construction")
    final_archive = verify_archive_state(modules["archive_guard"], root)
    if archive_state_invariant(root, final_archive) != contract["archive_reference_state_guard"]:
        raise CollationError("archive/reference identity changed before canonical publication")
    verify_pair_input_identities(root, pairs)
    publication = publish_package(root, fingerprint, crash_after_links=crash_after_links)
    verified = verify_results(root)
    return {"package": package, "publication": publication, "verified": verified}


def _manifest_file_record(manifest: Mapping[str, Any], relative: Path) -> Mapping[str, Any]:
    matches = [
        record for record in manifest.get("files", [])
        if isinstance(record, Mapping) and record.get("path") == str(relative)
    ]
    if len(matches) != 1:
        raise CollationError(f"package manifest lacks one file record for {relative}")
    return matches[0]


@contextmanager
def _stable_package_text(
    root: Path, fingerprint: str, manifest: Mapping[str, Any], relative: Path,
) -> Iterator[TextIO]:
    record = _manifest_file_record(manifest, relative)
    package_relative = PACKAGE_ROOT / fingerprint / relative
    expected = {
        "path": str(package_relative), "bytes": record["bytes"], "sha256": record["sha256"],
    }
    with stable_plain_text(
        root, package_relative, expected, allow_empty=int(record["bytes"]) == 0,
    ) as handle:
        yield handle


def _pairs_from_result_provenance(root: Path, provenance: Mapping[str, Any]) -> list[PairInput]:
    raw_inputs = provenance.get("inputs")
    if not isinstance(raw_inputs, list) or len(raw_inputs) != len(PAIR_ORDER):
        raise CollationError("result provenance lacks the exact three-pair input family")
    pairs: list[PairInput] = []
    for wanted, record in zip(PAIR_ORDER, raw_inputs, strict=True):
        if not isinstance(record, Mapping) or record.get("pair_id") != wanted:
            raise CollationError("result provenance input order is not A/B/CONTROL")
        ledger, pair_provenance = record.get("ledger"), record.get("provenance")
        if not isinstance(ledger, Mapping) or not isinstance(pair_provenance, Mapping):
            raise CollationError("result provenance input record lacks sealed identities")
        if not identity_matches(root, ledger) or not identity_matches(root, pair_provenance):
            raise CollationError(f"published result input identity drifted for {wanted}")
        try:
            expected_rows = int(record.get("expected_rows", -1))
        except (TypeError, ValueError) as error:
            raise CollationError("result provenance has an invalid pair denominator") from error
        pair = PairInput(
            pair_id=wanted, ledger=root / str(ledger["path"]),
            provenance=root / str(pair_provenance["path"]),
            terminal_state=str(record.get("terminal_state", "")),
            expected_rows=expected_rows, ledger_identity=dict(ledger),
            provenance_identity=dict(pair_provenance),
        )
        _validate_pair_provenance(pair, read_json(root, str(pair_provenance["path"])))
        pairs.append(pair)
    if provenance.get("input_sha256") != _digest_json(_pair_input_payload(pairs)):
        raise CollationError("result provenance complete-input family digest drifted")
    return pairs


def _verify_complete_07_stream(
    root: Path, pairs: Sequence[PairInput], expected: Mapping[str, Any],
) -> dict[str, object]:
    canonical_expected = {
        "path": str(PLACO_07_PATH), "bytes": expected["bytes"], "sha256": expected["sha256"],
    }
    compared = 0
    with stable_plain_text(root, PLACO_07_PATH, canonical_expected) as output_handle:
        output_reader = csv.DictReader(output_handle, delimiter="\t")
        if list(output_reader.fieldnames or []) != PLACO_07_FIELDS:
            raise CollationError("canonical complete 07 schema drifted")
        for pair in pairs:
            role = PAIR_IDENTITIES[pair.pair_id][2]
            with stable_gzip_text(
                root, pair.ledger.relative_to(root), pair.ledger_identity,
            ) as ledger_handle:
                input_reader = csv.DictReader(ledger_handle, delimiter="\t")
                if list(input_reader.fieldnames or []) != PLACO_LEDGER_FIELDS:
                    raise CollationError("authoritative PLACO ledger schema drifted during 07 audit")
                pair_rows = 0
                for pair_rows, source in enumerate(input_reader, start=1):
                    observed = next(output_reader, None)
                    if observed is None:
                        raise CollationError("canonical 07 truncates the complete PLACO family")
                    expected_row = {
                        "analysis_id": source["analysis_id"], "pair_id": pair.pair_id,
                        "family_role": role,
                        **{field: source[field] for field in PLACO_LEDGER_FIELDS[2:]},
                    }
                    if observed != expected_row:
                        raise CollationError(
                            f"canonical 07 differs from authoritative {pair.pair_id} row {pair_rows}"
                        )
                    compared += 1
                if pair_rows != pair.expected_rows:
                    raise CollationError("authoritative pair row count changed during 07 comparison")
        if next(output_reader, None) is not None:
            raise CollationError("canonical 07 has rows outside the A/B/CONTROL complete family")
    if compared != sum(pair.expected_rows for pair in pairs):
        raise CollationError("canonical 07 complete-family denominator did not reconcile")
    return portable_identity(root, PLACO_07_PATH)


def _verify_plink_input_files(
    root: Path, fingerprint: str, manifest: Mapping[str, Any],
    connection: sqlite3.Connection, pair_id: str,
) -> int:
    association_relative = Path(f"plink/{pair_id}/eligible.assoc.tsv")
    extract_relative = Path(f"plink/{pair_id}/eligible.snps.txt")
    cursor = connection.execute(
        "SELECT snp,p FROM eligible WHERE pair_id=? AND reference_status='MATCHED' "
        "ORDER BY CAST(p AS REAL),source_row_index", (pair_id,),
    )
    count = 0
    with (
        _stable_package_text(root, fingerprint, manifest, association_relative) as assoc_handle,
        _stable_package_text(root, fingerprint, manifest, extract_relative) as extract_handle,
    ):
        reader = csv.DictReader(assoc_handle, delimiter="\t")
        if list(reader.fieldnames or []) != ["SNP", "P"]:
            raise CollationError("packaged PLINK association schema drifted")
        for count, (snp, p_value) in enumerate(cursor, start=1):
            observed = next(reader, None)
            extracted = extract_handle.readline()
            if observed != {"SNP": str(snp), "P": str(p_value)}:
                raise CollationError("packaged PLINK association input differs from eligible union")
            if extracted != f"{snp}\n":
                raise CollationError("packaged PLINK extract input differs from eligible union")
        if next(reader, None) is not None or extract_handle.readline() != "":
            raise CollationError("packaged PLINK inputs contain rows outside eligible matched union")
    return count


def _validate_packaged_clumps(
    root: Path, fingerprint: str, manifest: Mapping[str, Any],
    connection: sqlite3.Connection,
) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    for pair_id in PAIR_ORDER:
        matched = _verify_plink_input_files(root, fingerprint, manifest, connection, pair_id)
        if matched == 0:
            results[pair_id] = {
                "pair_id": pair_id, "matched_rows": 0, "leads": [], "clumps": [],
                "commands": [], "pairwise": {
                    "expected_pairs": 0, "observed_pairs": 0, "maximum_r2": None,
                    "status": "PASS_NO_REFERENCE_MATCHED_ELIGIBLE_SIGNAL",
                },
            }
            continue
        clumped_relative = Path(f"plink/{pair_id}/clump.clumped")
        with _stable_package_text(root, fingerprint, manifest, clumped_relative) as handle:
            clumps = parse_plink_clumped_text(handle.read())
        assigned: set[str] = set()
        lead_coordinates: dict[str, tuple[int, int]] = {}
        for clump in clumps:
            member_p: list[float] = []
            for member in clump["members"]:
                found = connection.execute(
                    "SELECT chr,bp,p,reference_status FROM eligible WHERE pair_id=? AND snp=?",
                    (pair_id, member),
                ).fetchall()
                if len(found) != 1 or found[0][3] != "MATCHED":
                    raise CollationError("packaged PLINK clump escaped exact-reference eligible union")
                member_p.append(float(found[0][2]))
            lead_row = connection.execute(
                "SELECT chr,bp,p FROM eligible WHERE pair_id=? AND snp=?",
                (pair_id, clump["lead"]),
            ).fetchone()
            if lead_row is None:
                raise CollationError("packaged PLINK lead has no eligible source")
            if (clump["chr"], clump["bp"]) != (int(lead_row[0]), int(lead_row[1])):
                raise CollationError("packaged PLINK lead coordinate differs from sealed BIM")
            if not math.isclose(clump["plink_p"], float(lead_row[2]), rel_tol=1e-5, abs_tol=1e-300):
                raise CollationError("packaged PLINK lead was not ordered using actual PLACO P")
            if any(value < float(lead_row[2]) for value in member_p):
                raise CollationError("packaged clump contains a lower-P member than its lead")
            assigned.update(clump["members"])
            lead_coordinates[clump["lead"]] = (int(lead_row[0]), int(lead_row[1]))
        matched_ids = {
            str(row[0]) for row in connection.execute(
                "SELECT snp FROM eligible WHERE pair_id=? AND reference_status='MATCHED'", (pair_id,),
            )
        }
        if assigned != matched_ids:
            raise CollationError("packaged PLINK clumps do not account for every matched eligible variant")
        expected_pairs = sum(
            left[0] == right[0] and abs(left[1] - right[1]) <= 1_000_000
            for offset, left in enumerate(lead_coordinates.values())
            for right in list(lead_coordinates.values())[offset + 1:]
        )
        if expected_pairs:
            ld_relative = Path(f"plink/{pair_id}/lead_pairwise.ld")
            with _stable_package_text(root, fingerprint, manifest, ld_relative) as handle:
                pairwise = parse_pairwise_r2_text(handle.read(), lead_coordinates)
        else:
            pairwise = {
                "expected_pairs": 0, "observed_pairs": 0, "maximum_r2": None,
                "status": "PASS_NO_WITHIN_WINDOW_LEAD_PAIR",
            }
        connection.execute(
            "UPDATE eligible SET clump_status='PLINK_CLUMP_UNACCOUNTED',"
            "failure_reason='PLINK_SUCCESS_BUT_ELIGIBLE_MATCHED_VARIANT_NOT_IN_CLUMP_OUTPUT' "
            "WHERE pair_id=? AND reference_status='MATCHED'", (pair_id,),
        )
        for clump in clumps:
            lead = clump["lead"]
            connection.execute(
                "UPDATE eligible SET clump_status='INDEPENDENT_LEAD_VERIFIED_R2_LT_0.1',"
                "failure_reason='NONE' WHERE pair_id=? AND snp=?", (pair_id, lead),
            )
            connection.executemany(
                "UPDATE eligible SET clump_status=?,failure_reason='NONE' "
                "WHERE pair_id=? AND snp=?",
                [
                    (f"CLUMPED_SECONDARY_TO:{lead}", pair_id, member)
                    for member in clump["members"] if member != lead
                ],
            )
        connection.commit()
        results[pair_id] = {
            "pair_id": pair_id, "matched_rows": matched,
            "leads": [record["lead"] for record in clumps], "clumps": clumps,
            "commands": [], "pairwise": pairwise, "missing_from_clump": [],
        }
    return results


def _verify_tsv_against_rows(
    root: Path, relative: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]],
    expected_identity: Mapping[str, Any],
) -> int:
    expected = {
        "path": str(relative), "bytes": expected_identity["bytes"],
        "sha256": expected_identity["sha256"],
    }
    count = 0
    with stable_plain_text(root, relative, expected) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if list(reader.fieldnames or []) != list(fields):
            raise CollationError(f"canonical TSV schema drifted: {relative}")
        for count, wanted in enumerate(rows, start=1):
            observed = next(reader, None)
            exact = {field: str(wanted[field]) for field in fields}
            if observed != exact:
                raise CollationError(f"canonical TSV differs at {relative} row {count}")
        if next(reader, None) is not None:
            raise CollationError(f"canonical TSV contains extra rows: {relative}")
    return count


def _read_exact_science_tsv(
    root: Path, relative: Path, fields: Sequence[str], expected_identity: Mapping[str, Any],
) -> list[dict[str, str]]:
    expected = {
        "path": str(relative), "bytes": expected_identity["bytes"],
        "sha256": expected_identity["sha256"],
    }
    rows: list[dict[str, str]] = []
    with stable_plain_text(root, relative, expected) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if list(reader.fieldnames or []) != list(fields):
            raise CollationError(f"canonical science TSV schema drifted: {relative}")
        for row_number, row in enumerate(reader, start=1):
            if None in row or any(value is None for value in row.values()):
                raise CollationError(f"canonical science TSV is ragged at {relative}:{row_number}")
            rows.append(dict(row))
    return rows


def _verify_command_provenance(
    provenance: Mapping[str, Any], results: Mapping[str, Mapping[str, Any]],
    manifest: Mapping[str, Any],
) -> None:
    recorded = provenance.get("plink_execution")
    if not isinstance(recorded, Mapping) or set(recorded) != set(PAIR_ORDER):
        raise CollationError("result provenance lacks all three PLINK execution records")
    recorded_artifacts = provenance.get("plink_artifacts")
    expected_artifacts = [
        record for record in manifest.get("files", [])
        if isinstance(record, Mapping) and str(record.get("path", "")).startswith("plink/")
    ]
    if recorded_artifacts != expected_artifacts:
        raise CollationError("result provenance does not bind every packaged PLINK input/log/output")
    for pair_id in PAIR_ORDER:
        pair = recorded[pair_id]
        if not isinstance(pair, Mapping):
            raise CollationError("PLINK execution record is not an object")
        observed = results[pair_id]
        if (
            pair.get("matched_rows") != observed["matched_rows"]
            or pair.get("lead_count") != len(observed["leads"])
            or pair.get("missing_from_clump", [])
            or pair.get("pairwise") != observed["pairwise"]
        ):
            raise CollationError("PLINK execution provenance differs from reparsed artifacts")
        commands = pair.get("commands")
        expected_commands = int(observed["matched_rows"] > 0) + int(
            observed["pairwise"]["expected_pairs"] > 0
        )
        if not isinstance(commands, list) or len(commands) != expected_commands:
            raise CollationError("PLINK command count differs from required clump/r2 execution")
        for offset, command in enumerate(commands):
            if not isinstance(command, Mapping) or command.get("returncode") != 0:
                raise CollationError("PLINK command provenance lacks a successful exit")
            if (
                not isinstance(command.get("peak_aggregate_rss_bytes"), int)
                or command["peak_aggregate_rss_bytes"] <= 0
                or command.get("rss_measurement_method") not in {
                    "LIVE_PARENT_PLUS_OWNED_FRESH_SESSION_PGID_PS_SAMPLED",
                    "NORMALIZED_SELF_PLUS_CHILD_RUSAGE_CONSERVATIVE_FAST_EXIT_FALLBACK",
                }
            ):
                raise CollationError("PLINK command lacks authoritative aggregate RSS evidence")
            argv = command.get("argv")
            transcript = command.get("transcript")
            if not isinstance(argv, list) or not all(isinstance(value, str) for value in argv):
                raise CollationError("PLINK command argv is invalid")
            if not argv or Path(argv[0]).name != "plink":
                raise CollationError("PLINK command did not use the frozen executable link")
            if not isinstance(transcript, Mapping):
                raise CollationError("PLINK command transcript identity is invalid")
            manifest_record = _manifest_file_record(manifest, Path(str(transcript.get("path", ""))))
            if any(
                manifest_record.get(key) != transcript.get(key)
                for key in ("path", "bytes", "sha256")
            ):
                raise CollationError("PLINK transcript differs from package manifest")
            options = {value: argv[index + 1] for index, value in enumerate(argv[:-1]) if value.startswith("--")}
            if not str(options.get("--bfile", "")).endswith("/runtime_reference/g1000_eur"):
                raise CollationError("PLINK command did not use the frozen reference prefix")
            if offset == 0:
                required = {
                    "--clump-p1": "1", "--clump-p2": "1", "--clump-r2": "0.1",
                    "--clump-kb": "1000", "--threads": "1",
                    "--clump-snp-field": "SNP", "--clump-field": "P",
                }
                if any(options.get(key) != value for key, value in required.items()):
                    raise CollationError("recorded PLINK clump command differs from frozen settings")
                if "--clump" not in argv or "--bfile" not in argv or "--extract" not in argv:
                    raise CollationError("recorded PLINK clump command lacks required real inputs")
                if (
                    not str(options.get("--clump", "")).endswith(
                        f"/plink/{pair_id}/eligible.assoc.tsv"
                    )
                    or not str(options.get("--extract", "")).endswith(
                        f"/plink/{pair_id}/eligible.snps.txt"
                    )
                    or not str(options.get("--out", "")).endswith(f"/plink/{pair_id}/clump")
                ):
                    raise CollationError("recorded PLINK clump paths escaped the pair stage")
                _manifest_file_record(manifest, Path(f"plink/{pair_id}/clump.log"))
            else:
                required = {
                    "--ld-window": "99999999", "--ld-window-kb": "1000",
                    "--ld-window-r2": "0", "--threads": "1",
                }
                if (
                    any(options.get(key) != value for key, value in required.items())
                    or "--r2" not in argv or "--ld-snp-list" not in argv
                    or "--bfile" not in argv or "--extract" not in argv
                ):
                    raise CollationError("recorded PLINK r2 command differs from frozen settings")
                if (
                    not str(options.get("--extract", "")).endswith(
                        f"/plink/{pair_id}/leads.snps.txt"
                    )
                    or not str(options.get("--ld-snp-list", "")).endswith(
                        f"/plink/{pair_id}/leads.snps.txt"
                    )
                    or not str(options.get("--out", "")).endswith(
                        f"/plink/{pair_id}/lead_pairwise"
                    )
                ):
                    raise CollationError("recorded PLINK r2 paths escaped the pair stage")
                _manifest_file_record(manifest, Path(f"plink/{pair_id}/lead_pairwise.log"))


def _verify_result_output_records(
    root: Path, fingerprint: str, manifest: Mapping[str, Any],
    provenance: Mapping[str, Any], mapping: Mapping[Path, Path],
) -> None:
    raw = provenance.get("outputs")
    if not isinstance(raw, list):
        raise CollationError("result provenance output family is missing")
    by_canonical: dict[str, Mapping[str, Any]] = {}
    for record in raw:
        if not isinstance(record, Mapping) or not isinstance(record.get("canonical_path"), str):
            raise CollationError("result provenance has an invalid output record")
        if record["canonical_path"] in by_canonical:
            raise CollationError("result provenance duplicates a canonical output")
        by_canonical[record["canonical_path"]] = record
    expected_canonical = {str(path) for path in mapping if path != RESULT_PROVENANCE_PATH}
    if set(by_canonical) != expected_canonical:
        raise CollationError("result provenance output family is incomplete")
    for canonical, packaged in mapping.items():
        if canonical == RESULT_PROVENANCE_PATH:
            continue
        result_record = by_canonical[str(canonical)]
        manifest_record = _manifest_file_record(manifest, packaged)
        if any(
            result_record.get(key) != manifest_record.get(key)
            for key in ("path", "bytes", "sha256")
        ):
            raise CollationError("result provenance output differs from immutable package")
        if not _same_regular_inode(
            root / canonical, root / PACKAGE_ROOT / fingerprint / packaged,
        ):
            raise CollationError("canonical output is not the package's exact immutable inode")


def verify_published_link_family(
    root: Path, fingerprint: str, manifest: Mapping[str, Any],
) -> dict[Path, Path]:
    """Verify this stage's links without rejecting later authorized consumers."""
    mapping = _manifest_mapping(manifest)
    for canonical, packaged in mapping.items():
        if not _same_regular_inode(
            root / canonical, root / PACKAGE_ROOT / fingerprint / packaged,
        ):
            raise CollationError(f"canonical result is not package-linked: {canonical}")
    forbidden = [
        path for path in CANONICAL_SCIENCE_PATHS - set(mapping)
        if _artifact_present(root / path)
    ]
    if forbidden:
        raise CollationError("canonical result family contains science-state-incompatible artifacts")
    return mapping


def _verify_resource_artifacts(
    root: Path, provenance: Mapping[str, Any], mapping: Mapping[Path, Path],
    manifest: Mapping[str, Any],
) -> None:
    ram_record = _manifest_file_record(manifest, mapping[RAM_COMPONENT_PATH])
    rows = _read_exact_science_tsv(root, RAM_COMPONENT_PATH, RAM_FIELDS, ram_record)
    if len(rows) != 1:
        raise CollationError("component RAM benchmark must contain exactly one aggregate row")
    row = rows[0]
    try:
        n_snps, peak, runtime, exit_status = (
            int(row["n_snps"]), float(row["peak_ram_gb"]),
            float(row["runtime_sec"]), int(row["exit_status"]),
        )
    except ValueError as error:
        raise CollationError("component RAM benchmark has invalid numeric fields") from error
    if (
        row["analysis"] != "TRACK_B_PLEIOTROPY_POST_PLACO_COLLATION_V2"
        or row["pair"] != "A;B;CONTROL" or row["chromosome"] != "ALL"
        or n_snps != provenance.get("input_rows") or peak < 0 or runtime < 0
        or exit_status != 0 or row["output_hash"] != provenance.get("output_sha256")
    ):
        raise CollationError("component RAM benchmark differs from complete-family provenance")
    metrics = provenance.get("resource_metrics", {})
    try:
        metric_peak = int(metrics.get("peak_ram_bytes", -1))
        sampled_peak = int(metrics.get("maximum_live_parent_plus_plink_pgid_rss_bytes", -1))
        parent_peak = int(metrics.get("normalized_parent_only_peak_rss_bytes", -1))
    except (AttributeError, TypeError, ValueError) as error:
        raise CollationError("resource provenance has invalid RSS metrics") from error
    backend = str(metrics.get("rss_measurement_backend", "")) if isinstance(metrics, Mapping) else ""
    if (
        not isinstance(metrics, Mapping)
        or not backend.startswith(
            "MAX_OF_NORMALIZED_PARENT_RUSAGE_AND_LIVE_PARENT_PLUS_OWNED_FRESH_SESSION_PGID_PS_SAMPLES"
        ) and backend != "NORMALIZED_PARENT_RUSAGE_NO_PLINK_CHILD_REQUIRED"
        or metrics.get("eight_gib_sequential_feasible_with_one_gib_reserve") is not True
        or metric_peak < sampled_peak or metric_peak < parent_peak
        or not abs(peak * 1024**3 - metric_peak) <= max(1.0, metric_peak * 1e-10)
    ):
        raise CollationError("resource provenance lacks authoritative parent-plus-PGID peak RSS")
    namespace = read_json(root, RAM_NAMESPACE_PATH)
    if (
        namespace.get("schema_version") != "track-b-ram-benchmark-component-namespace.1"
        or namespace.get("legacy_top_level_mutated") is not False
        or namespace.get("exact_schema") != RAM_FIELDS
        or namespace.get("component_benchmark", {}).get("sha256") != ram_record["sha256"]
    ):
        raise CollationError("RAM benchmark namespace is not collision-safe/federatable")
    report_record = _manifest_file_record(manifest, mapping[RAM_REPORT_PATH])
    expected = {
        "path": str(RAM_REPORT_PATH), "bytes": report_record["bytes"],
        "sha256": report_record["sha256"],
    }
    with stable_plain_text(root, RAM_REPORT_PATH, expected) as handle:
        report = handle.read()
    for phrase in (
        "Safe decomposition", "Non-decomposable scientific stages", "Measured peak",
        "8-GB feasibility", "End-to-end component runtime", "Larger-server stages",
    ):
        if phrase not in report:
            raise CollationError(f"RAM-aware report lacks required section/content: {phrase}")


def verify_results(root: Path = ROOT) -> dict[str, Any]:
    """Side-effect-free deep verifier for downstream script 145.

    Repository artifacts are only read.  A private system-temporary SQLite
    index is used so complete ledgers/07 are never materialized in RAM.
    """
    contract = verify_pre_result_contract(root)
    if not _artifact_present(root / RESULT_PROVENANCE_PATH):
        raise CollationError("post-PLACO result provenance is absent")
    provenance = read_json(root, RESULT_PROVENANCE_PATH)
    fingerprint = provenance.get("fingerprint")
    if not isinstance(fingerprint, str) or not SHA256_RE.fullmatch(fingerprint):
        raise CollationError("post-PLACO result provenance has an invalid fingerprint")
    manifest = verify_package(root, fingerprint)
    mapping = verify_published_link_family(root, fingerprint, manifest)
    if RESULT_PROVENANCE_PATH not in mapping or not _same_regular_inode(
        root / RESULT_PROVENANCE_PATH,
        root / PACKAGE_ROOT / fingerprint / mapping[RESULT_PROVENANCE_PATH],
    ):
        raise CollationError("result provenance is not the immutable package commit inode")

    required_scalar = {
        "schema_version": SCHEMA,
        "analysis_id": "track-b-v1.0-pleiotropy",
        "science_state": manifest["science_state"],
        "policy_sha256": contract["upstream"][str(POLICY)]["sha256"],
        "contract_lock_sha256": contract["upstream"][str(PLEIOTROPY_CONTRACT_LOCK)]["sha256"],
        "input_gate_lock_sha256": contract["upstream"][str(INPUT_GATE_LOCK)]["sha256"],
        "software_sha256_or_commit": contract["upstream"][str(SCRIPT)]["sha256"],
        "reference_sha256": contract["ld_reference"]["bim"]["sha256"],
        "exact_schema": ",".join(PLACO_07_FIELDS),
        "qc_status": "PASS",
    }
    if any(provenance.get(key) != value for key, value in required_scalar.items()):
        raise CollationError("result provenance differs from frozen contract/science state")
    if provenance.get("fingerprint") != fingerprint:
        raise CollationError("result provenance fingerprint drifted")
    if provenance.get("thresholds") != contract["thresholds"]:
        raise CollationError("result provenance thresholds drifted")
    upstream_family = provenance.get("complete_143_family_verification", {})
    if (
        not isinstance(upstream_family, Mapping)
        or upstream_family.get("verifier") != contract["upstream"][str(SEQUENTIAL_SCRIPT)]
        or upstream_family.get("outcomes")
        != [{"pair": pair, "outcome": "VERIFIED_COMPLETE_CLEAN"} for pair in PAIR_ORDER]
        or upstream_family.get("performed_after_pre_result_contract") is not True
        or upstream_family.get(
            "all_valid_variants_and_global_bh_denominators_deeply_validated"
        ) is not True
    ):
        raise CollationError("result provenance lost the complete 143 family verification")
    if (
        provenance.get("collation_contract") != portable_identity(root, COLLATION_CONTRACT)
        or provenance.get("publication", {}).get(
            "competing_finemapping_artifacts_absent_at_commit_required"
        ) is not True
        or provenance.get("publication", {}).get("result_provenance_linked_last") is not True
    ):
        raise CollationError("result provenance lost its pre-result/commit-order attestation")
    if provenance.get("ld_reference") != contract["ld_reference"]:
        raise CollationError("result provenance LD reference identity drifted")
    runtime_audit = provenance.get("runtime_reference_audit", {})
    if (
        not isinstance(runtime_audit, Mapping)
        or runtime_audit.get("reference")
        != {suffix: contract["ld_reference"][suffix] for suffix in ("bed", "bim", "fam")}
        or runtime_audit.get("plink")
        != {
            key: contract["ld_reference"]["plink"][key]
            for key in ("path", "bytes", "sha256")
        }
        or runtime_audit.get("snapshot_mode")
        != (
            "HARDLINKED_EXACT_INODES_BEFORE_BIM_SCAN_AND_PLINK;FULL_HASH_AND_"
            "INODE_MTIME_CTIME_RECHECK_AFTER_ALL_COMMANDS"
        )
    ):
        raise CollationError("result provenance lost its exact runtime reference inode audit")
    if provenance.get("official_lava_partition") != contract["official_lava_partition"]:
        raise CollationError("result provenance official block identity drifted")
    if provenance.get("lava_failed_qc_accounting") != contract["lava_failed_qc_accounting"]:
        raise CollationError("result provenance LAVA failed-QC accounting drifted")

    modules = _load_production_modules(root)
    current_archive = verify_archive_state(modules["archive_guard"], root)
    current_archive_invariant = archive_state_invariant(root, current_archive)
    if (
        current_archive_invariant != contract["archive_reference_state_guard"]
        or provenance.get("archive_reference_state_guard") != current_archive_invariant
    ):
        raise CollationError("current archive state is not one of the exact sealed identities")
    recorded_archive = provenance.get("archive_reference_state")
    if (
        not isinstance(recorded_archive, Mapping)
        or recorded_archive.get("state") not in ACCEPTED_ARCHIVE_STATES
        or archive_state_invariant(root, recorded_archive) != current_archive_invariant
    ):
        raise CollationError("result-time archive state is invalid")
    current_conjfdr = validate_conjfdr_blocked(root)
    if provenance.get("conjfdr") != current_conjfdr:
        raise CollationError("conjFDR blocked-method attestation drifted")

    pairs = _pairs_from_result_provenance(root, provenance)
    expected_total = sum(pair.expected_rows for pair in pairs)
    if provenance.get("input_rows") != expected_total or provenance.get("output_rows") != expected_total:
        raise CollationError("result provenance complete-family row denominator drifted")
    package_07 = _manifest_file_record(manifest, mapping[PLACO_07_PATH])
    with tempfile.TemporaryDirectory(prefix="track_b_147_verify_") as temporary:
        database = Path(temporary) / "eligible.sqlite"
        connection = create_index_database(database)
        try:
            family_counts, projected_07 = index_complete_placo_family(root, pairs, connection)
            if family_counts != provenance.get("complete_family_counts"):
                raise CollationError("complete PLACO family counts differ from result provenance")
            if projected_07 != package_07["bytes"]:
                raise CollationError("complete 07 exact projected bytes drifted")
            complete_07_identity = _verify_complete_07_stream(root, pairs, package_07)
            if (
                provenance.get("output_sha256") != complete_07_identity["sha256"]
                or provenance.get("output_rows") != expected_total
            ):
                raise CollationError("complete 07 identity differs from result provenance")

            reference_coverage = match_eligible_to_bim(
                root, connection, contract["ld_reference"]["bim"]["path"],
                contract["ld_reference"]["bim"],
            )
            if reference_coverage["bim_rows_streamed"] != contract["ld_reference"]["variant_count"]:
                raise CollationError("deep verifier did not stream the complete sealed BIM")
            clump_results = _validate_packaged_clumps(
                root, fingerprint, manifest, connection,
            )
            _verify_command_provenance(provenance, clump_results, manifest)
            coverage_record = _manifest_file_record(manifest, mapping[LD_COVERAGE_PATH])
            coverage_rows = _verify_tsv_against_rows(
                root, LD_COVERAGE_PATH, LD_COVERAGE_FIELDS, _coverage_rows(connection),
                coverage_record,
            )
            control_coverage_record = _manifest_file_record(
                manifest, mapping[CONTROL_COVERAGE_PATH],
            )
            control_coverage_rows = _verify_tsv_against_rows(
                root, CONTROL_COVERAGE_PATH, LD_COVERAGE_FIELDS,
                _coverage_rows(connection, pair_filter="CONTROL"), control_coverage_record,
            )
            eligible_rows = int(connection.execute("SELECT COUNT(*) FROM eligible").fetchone()[0])
            expected_control_eligible = int(connection.execute(
                "SELECT COUNT(*) FROM eligible WHERE pair_id='CONTROL'"
            ).fetchone()[0])
            if coverage_rows != eligible_rows or control_coverage_rows != expected_control_eligible:
                raise CollationError("LD coverage ledger denominator did not reconcile")
            status_counts = {
                str(status): int(count) for status, count in connection.execute(
                    "SELECT reference_status,COUNT(*) FROM eligible GROUP BY reference_status"
                )
            }
            provenance_coverage = provenance.get("reference_coverage", {})
            if (
                not isinstance(provenance_coverage, Mapping)
                or provenance_coverage.get("eligible_rows") != eligible_rows
                or provenance_coverage.get("bim_rows_streamed")
                != reference_coverage["bim_rows_streamed"]
                or provenance_coverage.get("status_counts") != status_counts
                or provenance_coverage.get("failure_rows")
                != eligible_rows - int(status_counts.get("MATCHED", 0))
            ):
                raise CollationError("LD reference coverage provenance did not reconcile")

            _, block_index = load_official_blocks(root, contract["official_lava_partition"])
            lava_index = load_lava_diagnostic_index(root)
            primary_rows, control_rows = build_shared_locus_rows(
                connection, clump_results, block_index,
                contract["ld_reference"]["chromosome_lengths"], lava_index,
            )
        finally:
            connection.close()

    total_leads = len(primary_rows) + len(control_rows)
    expected_state = (
        "ZERO_ELIGIBLE_SIGNAL_FAMILY" if eligible_rows == 0 else
        "BLOCKED_BY_LD_REFERENCE_COVERAGE" if total_leads == 0 else
        "COMPLETE_WITH_INDEPENDENT_LOCI"
    )
    if provenance.get("science_state") != expected_state:
        raise CollationError("post-PLACO terminal science state does not follow exact counts")
    if (
        provenance.get("eligible_signal_rows") != eligible_rows
        or provenance.get("independent_primary_loci") != len(primary_rows)
        or provenance.get("independent_control_loci") != len(control_rows)
        or provenance.get("primary_has_loci") is not bool(primary_rows)
        or provenance.get("control_has_loci") is not bool(control_rows)
    ):
        raise CollationError("result provenance locus/eligible counts drifted")

    if primary_rows:
        primary_record = _manifest_file_record(manifest, mapping[SHARED_LOCI_08_PATH])
        _verify_tsv_against_rows(
            root, SHARED_LOCI_08_PATH, SHARED_LOCI_08_FIELDS, primary_rows, primary_record,
        )
        comparison_record = _manifest_file_record(manifest, mapping[COMPARISON_09_PATH])
        _verify_tsv_against_rows(
            root, COMPARISON_09_PATH, COMPARISON_09_FIELDS,
            comparison_rows(primary_rows), comparison_record,
        )
    if control_rows:
        control_record = _manifest_file_record(manifest, mapping[CONTROL_08_PATH])
        _verify_tsv_against_rows(
            root, CONTROL_08_PATH, SHARED_LOCI_08_FIELDS, control_rows, control_record,
        )
        control_comparison_record = _manifest_file_record(manifest, mapping[CONTROL_09_PATH])
        _verify_tsv_against_rows(
            root, CONTROL_09_PATH, COMPARISON_09_FIELDS,
            comparison_rows(control_rows), control_comparison_record,
        )
    if any(row["pair_id"] == "CONTROL" for row in primary_rows):
        raise CollationError("CONTROL leaked into the primary A/B locus family")
    if any(row["pair_id"] != "CONTROL" for row in control_rows):
        raise CollationError("primary discovery pair leaked into CONTROL outputs")
    all_locus_ids = [row["clump_id"] for row in [*primary_rows, *control_rows]]
    if len(all_locus_ids) != len(set(all_locus_ids)):
        raise CollationError("independent PLACO leads were collapsed or duplicated across pairs")
    for row in [*primary_rows, *control_rows]:
        if (
            row["pairwise_r2_status"] not in {
                "PASS_STRICT_R2_LT_0.1_WITHIN_1_MB", "PASS_NO_WITHIN_WINDOW_LEAD_PAIR",
            }
            or f"CONJFDR:{CONJFDR_BLOCKED}" not in row["method_availability"]
            or "MEDIATION" not in row["claim_limit"]
        ):
            raise CollationError("shared-locus row lost LD/method/claim constraints")
    invalid_comparison_labels = {"PLACO_ONLY", "CONJFDR_ONLY", "PLACO_AND_CONJFDR"}
    if any(
        row["comparison_label"] in invalid_comparison_labels
        for row in comparison_rows([*primary_rows, *control_rows])
    ):
        raise CollationError("blocked conjFDR state emitted a forbidden cross-method label")

    terminal = provenance.get("terminal_receipt")
    zero_family = expected_state == "ZERO_ELIGIBLE_SIGNAL_FAMILY"
    if zero_family:
        receipt = read_json(root, ZERO_FAMILY_PATH)
        if receipt != terminal or receipt.get("state") != expected_state:
            raise CollationError("zero-family terminal receipt drifted")
    elif expected_state == "BLOCKED_BY_LD_REFERENCE_COVERAGE":
        receipt = read_json(root, BLOCKED_COVERAGE_PATH)
        if (
            receipt != terminal or receipt.get("state") != expected_state
            or receipt.get("coverage_failure_rows") != eligible_rows
        ):
            raise CollationError("LD-coverage-blocked terminal receipt drifted")
    elif terminal is not None:
        raise CollationError("complete locus family unexpectedly carries a zero/blocked receipt")

    _verify_result_output_records(root, fingerprint, manifest, provenance, mapping)
    _verify_resource_artifacts(root, provenance, mapping, manifest)
    verify_pair_input_identities(root, pairs)
    return {
        "schema_version": VERIFIER_SCHEMA,
        "state": expected_state, "fingerprint": fingerprint,
        "zero_family": {"is_zero": zero_family, "state": expected_state},
        "primary_leads": primary_rows,
        "control_leads": control_rows,
        "complete_family": {
            "order": list(PAIR_ORDER), "rows": expected_total,
            "canonical_07": complete_07_identity,
            "pairs": {
                pair.pair_id: {
                    "terminal_state": pair.terminal_state,
                    "rows": pair.expected_rows,
                    "ledger": dict(pair.ledger_identity),
                    "provenance": dict(pair.provenance_identity),
                    "counts": family_counts[pair.pair_id],
                }
                for pair in pairs
            },
        },
        "archive_state": current_archive["state"],
        "archive_reference": current_archive,
        "ld_reference_coverage": {
            "eligible_rows": eligible_rows, "status_counts": status_counts,
            "coverage_failure_rows": eligible_rows - int(status_counts.get("MATCHED", 0)),
        },
        "method_availability": {"PLACO": "COMPLETE", "CONJFDR": CONJFDR_BLOCKED},
    }


deep_verify_published_family = verify_results


def preflight(root: Path = ROOT) -> dict[str, Any]:
    """Validate only pre-result dependencies; never open a PLACO result ledger."""
    assert_no_competing_finemapping(root)
    contract_path = root / COLLATION_CONTRACT
    if _artifact_present(contract_path):
        contract = verify_pre_result_contract(root)
        state = "SEALED_PRE_RESULT_CONTRACT_VERIFIED"
    else:
        assert_no_competing_publication(root)
        contract = build_pre_result_contract(root)
        state = "READY_TO_SEAL_PRE_RESULT_CONTRACT"
    return {
        "schema_version": CONTRACT_SCHEMA, "state": state,
        "placo_result_ledgers_opened": False,
        "pair_order": contract["pair_order"],
        "ld_reference": contract["ld_reference"],
        "archive_reference_state_guard": contract["archive_reference_state_guard"],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--preflight", action="store_true")
    action.add_argument("--materialize-ld", action="store_true")
    action.add_argument("--seal-contract", action="store_true")
    action.add_argument("--run", action="store_true")
    action.add_argument("--verify", action="store_true")
    parser.add_argument(
        "--execute", action="store_true",
        help="required for every production mutation (LD materialization, seal, or run)",
    )
    arguments = parser.parse_args(argv)
    try:
        if arguments.preflight:
            payload = preflight(ROOT)
        elif arguments.verify:
            payload = verify_results(ROOT)
        else:
            if not arguments.execute:
                raise CollationError("explicit --execute is required for production mutation")
            if arguments.materialize_ld:
                payload = materialize_ld_reference(ROOT)
            elif arguments.seal_contract:
                payload = seal_pre_result_contract(ROOT)
            else:
                payload = build_and_publish_results(ROOT)
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    except (CollationError, OSError, ValueError, KeyError) as error:
        print(f"TRACK_B_POST_PLACO_COLLATION_ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
