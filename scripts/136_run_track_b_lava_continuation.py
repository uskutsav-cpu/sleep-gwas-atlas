#!/usr/bin/env python3
"""Resume Track B LAVA after V1 discovery under an additive V2 identity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import re
import resource
import stat
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "scripts/134_track_b_lava_continuation_contract.py"
RUNNER = ROOT / "scripts/135_run_track_b_lava_postdiscovery_v2.R"
CHECKPOINT_VALIDATOR = ROOT / "scripts/137_validate_track_b_lava_checkpoint_v2.R"
RESULT_VALIDATOR = ROOT / "scripts/138_validate_track_b_lava_results_v2.py"
R_SCRIPT = ROOT / ".r-env/bin/Rscript"
RUNTIME_VALIDATOR = ROOT / "scripts/133_validate_track_b_lava_runtime.R"
BENCHMARK_ROOT = ROOT / "results/track_b/lava_continuations/benchmarks"
SUPERSESSION_VALIDATOR = ROOT / "scripts/144_freeze_track_b_lava_supersession.py"

EXECUTION_PHASES = {"aggregate-discovery", "conditional", "finalize"}
BUNDLE_PHASES = EXECUTION_PHASES | {"terminal-qc"}
LOCUS_PHASES = {"conditional"}
MARKER_PREFIX = "TRACK_B_LAVA_CONTINUATION_WORKER\t"
SEMANTIC_PREFIX = "TRACK_B_LAVA_CONTINUATION_CHECKPOINT\t"
RECEIPT_SCHEMA = "track-b-lava-continuation-receipt.2"
READY_SCHEMA = "track-b-lava-continuation-ready.1"
FINGERPRINT_RE = re.compile(r"^[0-9a-f]{64}$")
BENCHMARK_FIELDS = [
    "source_discovery_fingerprint", "continuation_execution_fingerprint", "analysis",
    "phase", "locus_index", "locus", "chromosome", "pair", "qc",
    "peak_ram_gb", "runtime_sec", "exit_status", "output_sha256",
]
CONDITIONAL_CANDIDATE_FIELDS = [
    "locus_index", "locus", "chromosome", "start", "stop", "n_snps",
    "eligible_models", "eligible_pairs",
]
SEMANTIC_ARTIFACT_NAMES = {
    "aggregate-discovery": (
        "conditional_candidates.tsv",
        "result.rds",
    ),
    "conditional": ("result.rds",),
    "finalize": (
        "04_lava_local_results.tsv",
        "05_local_conditional_results.tsv",
        "lava_bivariate.tsv",
        "lava_conditional.tsv",
        "lava_locus_status.tsv",
        "lava_univariate.tsv",
        "result.rds",
        "staging_validation.txt",
    ),
    "terminal-qc": (
        "lava_bivariate.tsv",
        "lava_conditional.tsv",
        "lava_locus_status.tsv",
        "lava_univariate.tsv",
        "result.rds",
        "terminal_qc.json",
    ),
}
PUBLICATION_METADATA_NAMES = {
    "READY", "receipt.json", "semantic_validation.txt", "worker.log",
}
_SOURCE_LOCK_CACHE: dict[tuple[str, int, str, str], dict[str, Any]] = {}


class ExecutionError(RuntimeError):
    pass


class TerminalQCError(ExecutionError):
    pass


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise ExecutionError(f"could not load module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


CONTRACT = _load_module("track_b_lava_continuation_contract", CONTRACT_PATH)


def file_identity(path: Path) -> tuple[int, str]:
    return CONTRACT.file_identity(path)


def sha256(path: Path) -> str:
    return file_identity(path)[1]


def fsync_directory(path: Path) -> None:
    CONTRACT.fsync_directory(path)


def publish_no_replace(path: Path, value: bytes) -> None:
    try:
        CONTRACT.publish_no_replace(path, value)
    except CONTRACT.ContinuationError as error:
        raise ExecutionError(str(error)) from error


def atomic_text(path: Path, value: str) -> None:
    if path == ROOT or ROOT in path.parents:
        CONTRACT.require_real_repository_path(path.parent, "atomic publication parent")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path == ROOT or ROOT in path.parents:
        CONTRACT.require_real_repository_path(path.parent, "atomic publication parent")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("x", encoding="utf-8") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def stable_stat(path: Path) -> dict[str, int]:
    try:
        if path.is_symlink():
            raise ExecutionError(f"active artifact may not be a symlink: {path}")
        before = path.stat()
        if not path.is_file() or before.st_size <= 0:
            raise ExecutionError(f"active artifact is missing or empty: {path}")
        after = path.stat()
    except OSError as error:
        raise ExecutionError(f"could not stat active artifact {path}: {error}") from error
    fields = lambda value: {
        "device": int(value.st_dev), "inode": int(value.st_ino), "bytes": int(value.st_size),
        "mtime_ns": int(value.st_mtime_ns), "ctime_ns": int(value.st_ctime_ns),
    }
    if fields(before) != fields(after):
        raise ExecutionError(f"active artifact changed while being captured: {path}")
    return fields(after)


def capture_semantic_artifacts(
    attempt: Path, phase: str, *, publication_metadata_allowed: bool = False,
) -> list[dict[str, object]]:
    """Identify exactly the scientific files consumed by semantic validation.

    The result is deliberately stronger than a checksum list: each member binds
    the direct file's device, inode, size, mtime, ctime, and SHA-256.  Comparing
    snapshots on both sides of semantic validation closes the otherwise-open
    interval in which different bytes could be receipted than were validated.
    """

    names = SEMANTIC_ARTIFACT_NAMES.get(phase)
    if names is None:
        raise ExecutionError(f"unsupported semantic-artifact phase: {phase}")
    if attempt.is_symlink() or not attempt.is_dir():
        raise ExecutionError("semantic artifact directory is missing or symlinked")
    live: set[str] = set()
    for path in attempt.iterdir():
        if path.is_symlink() or not path.is_file():
            raise ExecutionError(f"semantic artifact inventory contains a non-file: {path}")
        live.add(path.name)
    required = set(names)
    allowed = required | PUBLICATION_METADATA_NAMES if publication_metadata_allowed else required
    if not required.issubset(live) or not live.issubset(allowed):
        raise ExecutionError(
            f"semantic artifact inventory drifted for {phase}: "
            f"required={sorted(required)} allowed={sorted(allowed)} observed={sorted(live)}"
        )

    records: list[dict[str, object]] = []
    for name in names:
        path = attempt / name
        before = stable_stat(path)
        size, digest = file_identity(path)
        after = stable_stat(path)
        if before != after or size != before["bytes"]:
            raise ExecutionError(f"semantic artifact changed while being identified: {path}")
        records.append({
            "path": name,
            "bytes": size,
            "sha256": digest,
            "stable_stat": before,
        })
    return records


def require_semantic_artifacts_unchanged(
    attempt: Path, phase: str, expected: list[dict[str, object]], *,
    publication_metadata_allowed: bool = False,
) -> list[dict[str, object]]:
    observed = capture_semantic_artifacts(
        attempt, phase, publication_metadata_allowed=publication_metadata_allowed,
    )
    if observed != expected:
        raise ExecutionError(
            f"scientific artifact family changed across semantic validation/publication: {phase}"
        )
    return observed


def bind_semantic_validation(
    before: list[dict[str, object]], after: list[dict[str, object]], phase: str,
) -> list[dict[str, object]]:
    if before != after:
        raise ExecutionError(
            f"scientific artifact family changed across semantic validation: {phase}"
        )
    return [
        {
            "path": item["path"],
            "bytes": item["bytes"],
            "sha256": item["sha256"],
            "pre_semantic_stat": item["stable_stat"],
            "post_semantic_stat": observed["stable_stat"],
        }
        for item, observed in zip(before, after, strict=True)
    ]


def semantic_snapshot_from_binding(
    binding: list[dict[str, object]], phase: str,
) -> list[dict[str, object]]:
    expected_names = list(SEMANTIC_ARTIFACT_NAMES.get(phase, ()))
    required = {
        "path", "bytes", "sha256", "pre_semantic_stat", "post_semantic_stat",
    }
    if (
        len(binding) != len(expected_names)
        or [item.get("path") if isinstance(item, dict) else None for item in binding]
        != expected_names
    ):
        raise ExecutionError(f"semantic artifact binding family drifted: {phase}")
    snapshot: list[dict[str, object]] = []
    for item in binding:
        if not isinstance(item, dict) or set(item) != required:
            raise ExecutionError(f"semantic artifact binding schema drifted: {phase}")
        if item["pre_semantic_stat"] != item["post_semantic_stat"]:
            raise ExecutionError(f"semantic artifact stat identity changed during validation: {phase}")
        snapshot.append({
            "path": item["path"],
            "bytes": item["bytes"],
            "sha256": item["sha256"],
            "stable_stat": item["post_semantic_stat"],
        })
    return snapshot


def _relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError as error:
        raise ExecutionError(f"path escapes repository: {path}") from error


def _safe_relative(value: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ExecutionError(f"unsafe continuation dependency path: {value}")
    path = ROOT / relative
    try:
        path.resolve(strict=False).relative_to(ROOT.resolve())
    except ValueError as error:
        raise ExecutionError(f"continuation dependency escapes repository: {value}") from error
    return path


def run_root(continuation_fingerprint: str) -> Path:
    try:
        return CONTRACT.continuation_run_root(continuation_fingerprint)
    except CONTRACT.ContinuationError as error:
        raise ExecutionError(str(error)) from error


def bundle_name(phase: str, index: int | None) -> str:
    if phase == "conditional":
        if index is None or not 1 <= index <= 2495:
            raise ExecutionError("conditional phase requires a locus index in 1..2495")
        return f"locus_{index:04d}"
    if phase in {"aggregate-discovery", "finalize", "terminal-qc"} and index is None:
        return "unit_all"
    raise ExecutionError(f"unsupported continuation phase/index: {phase}/{index}")


def bundle_path(phase: str, index: int | None, continuation_fingerprint: str) -> Path:
    if phase not in BUNDLE_PHASES:
        raise ExecutionError(f"unsupported bundle phase: {phase}")
    return run_root(continuation_fingerprint) / phase / bundle_name(phase, index)


def _attempt_directory(phase: str, index: int | None, continuation_fingerprint: str) -> Path:
    namespace = run_root(continuation_fingerprint)
    namespace.mkdir(parents=True, exist_ok=True)
    CONTRACT.require_real_repository_path(namespace, "continuation run namespace")
    attempts = namespace / ".attempts"
    attempts.mkdir(parents=True, exist_ok=True)
    CONTRACT.require_real_repository_path(attempts, "continuation attempt namespace")
    if not attempts.is_dir():
        raise ExecutionError("continuation attempt namespace is not a directory")
    attempt = attempts / f"{phase}_{bundle_name(phase, index)}.{os.getpid()}.{time.time_ns()}"
    attempt.mkdir()
    CONTRACT.require_real_repository_path(attempt, "continuation attempt directory")
    fsync_directory(attempts)
    return attempt


def _safe_attempt_target(bundle: Path, continuation_fingerprint: str) -> Path:
    if not bundle.is_symlink():
        raise ExecutionError(f"continuation checkpoint is not an atomic READY link: {bundle}")
    namespace = run_root(continuation_fingerprint)
    CONTRACT.require_real_repository_path(namespace, "continuation run namespace")
    CONTRACT.require_real_repository_path(bundle.parent, "continuation bundle parent")
    attempt_path = namespace / ".attempts"
    CONTRACT.require_real_repository_path(attempt_path, "continuation attempt namespace")
    try:
        attempts = attempt_path.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ExecutionError("continuation attempt namespace is unreadable") from error
    if not attempt_path.is_dir() or attempts != attempt_path:
        raise ExecutionError("continuation attempt namespace is not a real in-tree directory")
    try:
        target = bundle.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ExecutionError(f"continuation checkpoint READY link is broken: {bundle}") from error
    if target.parent != attempts:
        raise ExecutionError(f"continuation checkpoint READY link escapes its run namespace: {bundle}")
    try:
        target_stat = os.lstat(target)
    except OSError as error:
        raise ExecutionError(f"continuation checkpoint target is unreadable: {bundle}") from error
    if not stat.S_ISDIR(target_stat.st_mode):
        raise ExecutionError(f"continuation checkpoint target is not a real directory: {bundle}")
    return target


def _publish_ready_bundle(
    attempt: Path, phase: str, index: int | None, continuation_fingerprint: str,
) -> None:
    bundle = bundle_path(phase, index, continuation_fingerprint)
    namespace = run_root(continuation_fingerprint)
    attempts = namespace / ".attempts"
    CONTRACT.require_real_repository_path(attempts, "continuation attempt namespace")
    CONTRACT.require_real_repository_path(attempt, "continuation attempt directory")
    if attempt.parent != attempts or not attempt.is_dir():
        raise ExecutionError("refusing to publish an attempt outside the real attempt namespace")
    CONTRACT.require_real_repository_path(bundle.parent, "continuation bundle parent")
    bundle.parent.mkdir(parents=True, exist_ok=True)
    CONTRACT.require_real_repository_path(bundle.parent, "continuation bundle parent")
    try:
        os.symlink(os.path.relpath(attempt, start=bundle.parent), bundle, target_is_directory=True)
    except FileExistsError as error:
        raise ExecutionError(f"continuation READY link appeared during publication: {bundle}") from error
    fsync_directory(bundle.parent)


def quick_contract(source_fingerprint: str, continuation_fingerprint: str) -> None:
    try:
        observed = CONTRACT.continuation_fingerprint(source_fingerprint)
        if observed != continuation_fingerprint:
            raise ExecutionError("live continuation code/runtime differs from supplied fingerprint")
        CONTRACT.validate_lineage_record(source_fingerprint, continuation_fingerprint)
    except CONTRACT.ContinuationError as error:
        raise ExecutionError(str(error)) from error


def _source_lock(source_fingerprint: str) -> dict[str, Any]:
    path = CONTRACT.source_lock_path(source_fingerprint)
    try:
        size, digest = file_identity(path)
        legacy_digest = CONTRACT._digest(CONTRACT._validate_legacy_code())
        key = (source_fingerprint, size, digest, legacy_digest)
        cached = _SOURCE_LOCK_CACHE.get(key)
        if cached is not None:
            return cached
        payload = CONTRACT.validate_source_lock(
            source_fingerprint, validate_bundles=False, require_live_fingerprint=False,
        )
        if file_identity(path) != (size, digest):
            raise ExecutionError("source-family lock changed while entering the validation cache")
        _SOURCE_LOCK_CACHE.clear()
        _SOURCE_LOCK_CACHE[key] = payload
        return payload
    except CONTRACT.ContinuationError as error:
        raise ExecutionError(str(error)) from error


def _artifact_spec(path: Path, role: str, expected: tuple[int, str] | None = None) -> dict[str, object]:
    size, digest = expected or file_identity(path)
    return {"path": _relative(path), "bytes": int(size), "sha256": digest, "role": role}


def _source_discovery_specs(source_fingerprint: str) -> list[dict[str, object]]:
    lock = _source_lock(source_fingerprint)
    specs = []
    for entry in lock["checkpoints"]:
        record = entry["result"]
        specs.append(_artifact_spec(
            _safe_relative(str(record["path"])), "SOURCE_DISCOVERY_CHECKPOINT",
            (int(record["bytes"]), str(record["sha256"])),
        ))
    if len(specs) != 2495:
        raise ExecutionError("source lock does not expose the exact 2,495-result family")
    return specs


def _source_discovery_spec(
    source_fingerprint: str, locus_index: int,
) -> dict[str, object]:
    if not 1 <= locus_index <= 2495:
        raise ExecutionError("source discovery dependency index is outside 1..2495")
    lock = _source_lock(source_fingerprint)
    entries = lock.get("checkpoints")
    if not isinstance(entries, list) or len(entries) != 2495:
        raise ExecutionError("source lock does not expose the exact 2,495-result family")
    entry = entries[locus_index - 1]
    if not isinstance(entry, dict) or entry.get("locus_index") != locus_index:
        raise ExecutionError("source lock local checkpoint identity drifted")
    record = entry.get("result")
    if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
        raise ExecutionError("source lock local result identity is malformed")
    return _artifact_spec(
        _safe_relative(str(record["path"])), "SOURCE_DISCOVERY_CHECKPOINT",
        (int(record["bytes"]), str(record["sha256"])),
    )


def _scientific_input_specs(
    phase: str, index: int | None, source_fingerprint: str, *,
    include_conditional_chromosome: bool = False,
) -> list[dict[str, object]]:
    """Return the exact live scientific inputs read by this phase and validator."""

    legacy = CONTRACT.legacy_supervisor()
    try:
        baseline = legacy.validate_input_baseline(source_fingerprint)
    except Exception as error:
        raise ExecutionError(f"source active-input baseline is invalid: {error}") from error
    entries = baseline.get("entries")
    if not isinstance(entries, list) or not entries:
        raise ExecutionError("source active-input baseline has no entries")

    chromosome: int | None = None
    if phase == "conditional":
        if index is None or not 1 <= index <= 2495:
            raise ExecutionError("conditional scientific inputs require one locus index")
        if include_conditional_chromosome:
            chromosome = int(legacy.load_loci()[index - 1]["CHR"])
    elif include_conditional_chromosome:
        raise ExecutionError("only a conditional candidate may bind chromosome inputs")

    selected: list[dict[str, object]] = []
    chromosome_counts = {
        "REFERENCE_CHROMOSOME": 0,
        "CHROMOSOME_SHARD": 0,
        "CHROMOSOME_INPUT_INFO": 0,
    }
    for item in entries:
        if not isinstance(item, dict):
            raise ExecutionError("source active-input baseline entry is malformed")
        role = str(item.get("role", ""))
        include = role == "COMMON_INPUT"
        if chromosome is not None and role in chromosome_counts:
            try:
                same_chromosome = int(item.get("chromosome")) == chromosome
            except (TypeError, ValueError):
                same_chromosome = False
            if same_chromosome:
                include = True
                chromosome_counts[role] += 1
        if include:
            required = {"path", "bytes", "sha256", "role"}
            if not required.issubset(item):
                raise ExecutionError("source active-input baseline entry lacks identity fields")
            selected.append({key: item[key] for key in ("path", "bytes", "sha256", "role")})

    if not any(item["role"] == "COMMON_INPUT" for item in selected):
        raise ExecutionError("continuation phase lacks its common scientific inputs")
    if chromosome is not None and chromosome_counts != {
        "REFERENCE_CHROMOSOME": 2,
        "CHROMOSOME_SHARD": 8,
        "CHROMOSOME_INPUT_INFO": 1,
    }:
        raise ExecutionError(
            f"chromosome {chromosome} active-input family is incomplete: {chromosome_counts}"
        )
    return selected


def _execution_code_specs() -> list[dict[str, object]]:
    paths = [
        CONTRACT_PATH,
        Path(__file__).resolve(),
        RUNNER,
        CHECKPOINT_VALIDATOR,
        RESULT_VALIDATOR,
        SUPERSESSION_VALIDATOR,
        RUNTIME_VALIDATOR,
        R_SCRIPT,
    ]
    return [_artifact_spec(path, "CONTINUATION_EXECUTABLE") for path in paths]


def _result_spec(
    phase: str, index: int | None, source_fingerprint: str, continuation_fingerprint: str,
) -> dict[str, object]:
    receipt = validate_bundle(
        phase, index, source_fingerprint, continuation_fingerprint,
        revalidate_semantics=False, reconcile_active=False,
    )
    if receipt is None:
        raise ExecutionError(f"missing continuation prerequisite: {phase}/{index}")
    target = _safe_attempt_target(bundle_path(phase, index, continuation_fingerprint), continuation_fingerprint)
    record = next((item for item in receipt["artifacts"] if item["path"] == "result.rds"), None)
    if record is None:
        raise ExecutionError(f"continuation prerequisite lacks result.rds: {phase}/{index}")
    return _artifact_spec(
        target / "result.rds", "CONTINUATION_CHECKPOINT",
        (int(record["bytes"]), str(record["sha256"])),
    )


def _candidate_loci(path: Path) -> frozenset[int]:
    """Parse the sealed aggregate candidate attestation fail-closed."""

    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != CONDITIONAL_CANDIDATE_FIELDS:
                raise ExecutionError("conditional candidate attestation schema drifted")
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise ExecutionError(f"conditional candidate attestation is unreadable: {error}") from error

    official = CONTRACT.legacy_supervisor().load_loci()
    if len(official) != 2495:
        raise ExecutionError("official locus family drifted while reading conditional candidates")
    integer_pattern = re.compile(r"^(?:0|[1-9][0-9]*)$")
    observed: list[int] = []
    for row_number, row in enumerate(rows, start=2):
        if set(row) != set(CONDITIONAL_CANDIDATE_FIELDS) or any(
            row[field] is None or not integer_pattern.fullmatch(row[field])
            for field in CONDITIONAL_CANDIDATE_FIELDS[:-1]
        ):
            raise ExecutionError(f"malformed conditional candidate row {row_number}")
        values = {field: int(row[field]) for field in CONDITIONAL_CANDIDATE_FIELDS[:-1]}
        index = values["locus_index"]
        if not 1 <= index <= 2495 or values["eligible_models"] <= 0 or values["n_snps"] <= 0:
            raise ExecutionError(f"invalid conditional candidate row {row_number}")
        pairs = row["eligible_pairs"].split(";")
        if pairs not in (["A"], ["B"], ["A", "B"]):
            raise ExecutionError(f"invalid conditional candidate pair family at row {row_number}")
        locus = official[index - 1]
        expected = {
            "locus": int(locus["LOC"]), "chromosome": int(locus["CHR"]),
            "start": int(locus["START"]), "stop": int(locus["STOP"]),
        }
        if any(values[field] != expected[field] for field in expected):
            raise ExecutionError(f"conditional candidate coordinates drifted at row {row_number}")
        if observed and index <= observed[-1]:
            raise ExecutionError("conditional candidates are duplicated or out of locus order")
        observed.append(index)
    return frozenset(observed)


def _aggregate_dependency_specs(
    source_fingerprint: str, continuation_fingerprint: str,
) -> tuple[list[dict[str, object]], frozenset[int]]:
    """Bind the aggregate result, candidate artifact, and sealing receipt."""

    receipt = validate_bundle(
        "aggregate-discovery", None, source_fingerprint, continuation_fingerprint,
        revalidate_semantics=False, reconcile_active=False,
    )
    if receipt is None:
        raise ExecutionError("missing continuation prerequisite: aggregate-discovery/None")
    target = _safe_attempt_target(
        bundle_path("aggregate-discovery", None, continuation_fingerprint),
        continuation_fingerprint,
    )
    records = {str(item["path"]): item for item in receipt["artifacts"]}
    required = {"result.rds", "conditional_candidates.tsv"}
    if not required.issubset(records):
        raise ExecutionError("aggregate prerequisite lacks its result or candidate attestation")
    candidate_path = target / "conditional_candidates.tsv"
    candidate_record = records["conditional_candidates.tsv"]
    candidate_identity = (int(candidate_record["bytes"]), str(candidate_record["sha256"]))
    if file_identity(candidate_path) != candidate_identity:
        raise ExecutionError("conditional candidate attestation differs from aggregate receipt")
    result_record = records["result.rds"]
    specs = [
        _artifact_spec(
            target / "result.rds", "CONTINUATION_CHECKPOINT",
            (int(result_record["bytes"]), str(result_record["sha256"])),
        ),
        _artifact_spec(
            candidate_path, "CONDITIONAL_CANDIDATE_ATTESTATION", candidate_identity,
        ),
        _artifact_spec(target / "receipt.json", "CONTINUATION_CHECKPOINT_RECEIPT"),
    ]
    return specs, _candidate_loci(candidate_path)


def active_specs(
    phase: str, index: int | None, source_fingerprint: str, continuation_fingerprint: str,
) -> list[dict[str, object]]:
    if phase not in EXECUTION_PHASES and phase != "terminal-qc":
        raise ExecutionError(f"unsupported active-input phase: {phase}")
    lineage = CONTRACT.lineage_path(continuation_fingerprint)
    source_lock = CONTRACT.source_lock_path(source_fingerprint)
    aggregate_specs: list[dict[str, object]] = []
    candidate_loci: frozenset[int] = frozenset()
    if phase != "aggregate-discovery":
        aggregate_specs, candidate_loci = _aggregate_dependency_specs(
            source_fingerprint, continuation_fingerprint,
        )
    specs = [
        _artifact_spec(lineage, "CONTINUATION_LINEAGE"),
        _artifact_spec(source_lock, "SOURCE_FAMILY_LOCK"),
        _artifact_spec(
            CONTRACT.predecessor_supersession_lock_path(),
            "PREDECESSOR_SUPERSESSION_LOCK",
        ),
    ]
    specs.extend(_execution_code_specs())
    specs.extend(_scientific_input_specs(
        phase, index, source_fingerprint,
        include_conditional_chromosome=(
            phase == "conditional" and index is not None and index in candidate_loci
        ),
    ))
    if phase == "aggregate-discovery":
        specs.extend(_source_discovery_specs(source_fingerprint))
    elif phase == "conditional":
        if index is None:
            raise ExecutionError("conditional continuation requires one source checkpoint")
        specs.append(_source_discovery_spec(source_fingerprint, index))
        specs.extend(aggregate_specs)
    else:
        specs.extend(_source_discovery_specs(source_fingerprint))
        specs.extend(aggregate_specs)
        specs.extend(
            _result_spec("conditional", locus_index, source_fingerprint, continuation_fingerprint)
            for locus_index in range(1, 2496)
        )
    specs.sort(key=lambda item: str(item["path"]))
    if len({str(item["path"]) for item in specs}) != len(specs):
        raise ExecutionError(f"duplicate active-input dependency in continuation {phase}")
    return specs


def capture_active_inputs(specs: list[dict[str, object]]) -> list[dict[str, object]]:
    captured = []
    for spec in specs:
        path = _safe_relative(str(spec["path"]))
        before = stable_stat(path)
        if before["bytes"] != int(spec["bytes"]) or sha256(path) != spec["sha256"]:
            raise ExecutionError(f"continuation dependency differs from expected identity: {spec['path']}")
        captured.append(spec | {"pre_stat": before})
    return captured


def complete_active_inputs(captured: list[dict[str, object]]) -> list[dict[str, object]]:
    completed = []
    for item in captured:
        path = _safe_relative(str(item["path"]))
        after = stable_stat(path)
        if after != item["pre_stat"] or sha256(path) != item["sha256"]:
            raise ExecutionError(f"continuation dependency changed during execution: {item['path']}")
        completed.append(item | {"post_stat": after})
    return completed


def validate_phase_prerequisites(
    phase: str, source_fingerprint: str, continuation_fingerprint: str,
) -> None:
    if phase == "aggregate-discovery":
        _source_lock(source_fingerprint)
    elif phase == "conditional":
        if validate_bundle(
            "aggregate-discovery", None, source_fingerprint, continuation_fingerprint,
            revalidate_semantics=False, reconcile_active=False,
        ) is None:
            raise ExecutionError("conditional continuation lacks full-family discovery aggregation")
    elif phase in {"finalize", "terminal-qc"}:
        if validate_bundle(
            "aggregate-discovery", None, source_fingerprint, continuation_fingerprint,
        ) is None:
            raise ExecutionError("finalization lacks full-family discovery aggregation")
        for locus_index in range(1, 2496):
            if validate_bundle(
                "conditional", locus_index, source_fingerprint, continuation_fingerprint,
                revalidate_semantics=False, reconcile_active=False,
            ) is None:
                raise ExecutionError(f"finalization lacks conditional checkpoint {locus_index}")
    else:
        raise ExecutionError(f"unsupported continuation prerequisite phase: {phase}")


def parse_marker(log: str) -> dict[str, str]:
    lines = [line for line in log.splitlines() if line.startswith(MARKER_PREFIX)]
    if len(lines) != 1:
        raise ExecutionError("continuation worker did not emit exactly one metadata marker")
    marker: dict[str, str] = {}
    for field in lines[0].split("\t")[1:]:
        if "=" not in field:
            raise ExecutionError("malformed continuation worker marker")
        key, value = field.split("=", 1)
        marker[key] = value
    required = {
        "phase", "index", "locus", "chromosome", "n_snps", "pair", "qc",
        "construction_complete", "output", "source_fingerprint", "continuation_fingerprint",
    }
    if set(marker) != required:
        raise ExecutionError("continuation worker marker fields drifted")
    return marker


def validate_marker(
    marker: dict[str, str], phase: str, index: int | None, result_path: Path,
    source_fingerprint: str, continuation_fingerprint: str,
) -> None:
    expected_index = 0 if index is None else index
    if (
        marker["phase"] != phase or marker["index"] != str(expected_index)
        or marker["source_fingerprint"] != source_fingerprint
        or marker["continuation_fingerprint"] != continuation_fingerprint
    ):
        raise ExecutionError("continuation worker marker identity drifted")
    if phase == "conditional":
        locus = CONTRACT.legacy_supervisor().load_loci()[expected_index - 1]
        if marker["locus"] != str(locus["LOC"]) or marker["chromosome"] != str(locus["CHR"]):
            raise ExecutionError("conditional worker marker has wrong official locus")
    elif marker["locus"] != "ALL" or marker["chromosome"] != "ALL":
        raise ExecutionError("aggregate continuation marker is malformed")
    try:
        marker_output = Path(marker["output"]).resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ExecutionError("continuation marker names an invalid result") from error
    if marker_output != result_path.resolve(strict=True):
        raise ExecutionError("continuation marker result differs from supervised attempt")
    if marker["construction_complete"] not in {"TRUE", "FALSE"}:
        raise ExecutionError("continuation construction state is invalid")


def _semantic_validation(
    phase: str, index: int | None, source_fingerprint: str,
    continuation_fingerprint: str, result_path: Path,
) -> tuple[dict[str, str], str]:
    command = [
        str(R_SCRIPT), str(CHECKPOINT_VALIDATOR), "--phase", phase,
        "--source-fingerprint", source_fingerprint,
        "--continuation-fingerprint", continuation_fingerprint,
        "--source-checkpoint-root", str(SOURCE_CHECKPOINT_ROOT(source_fingerprint)),
        "--continuation-checkpoint-root", str(run_root(continuation_fingerprint)),
        "--rds", str(result_path),
    ]
    if index is not None:
        command.extend(["--locus-index", str(index)])
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    combined = "# command\n" + " ".join(command) + "\n# stdout\n" + result.stdout + "# stderr\n" + result.stderr
    markers = [line for line in result.stdout.splitlines() if line.startswith(SEMANTIC_PREFIX)]
    if result.returncode != 0 or len(markers) != 1:
        raise ExecutionError(f"continuation checkpoint failed semantic validation:\n{combined}")
    values: dict[str, str] = {}
    for field in markers[0].split("\t")[1:]:
        if "=" not in field:
            raise ExecutionError("malformed continuation semantic marker")
        key, value = field.split("=", 1)
        values[key] = value
    required = {
        "phase", "index", "locus", "chromosome", "n_snps", "pair", "qc",
        "construction_complete", "source_fingerprint", "continuation_fingerprint", "status",
    }
    if set(values) != required or values["status"] != "PASS":
        raise ExecutionError("continuation semantic marker scope drifted")
    if phase == "terminal-qc":
        terminal_command = [
            sys.executable, str(RESULT_VALIDATOR), "--verify-terminal-qc",
            "--staging-dir", str(result_path.parent),
            "--source-fingerprint", source_fingerprint,
            "--continuation-fingerprint", continuation_fingerprint,
        ]
        terminal = subprocess.run(
            terminal_command, cwd=ROOT, text=True, capture_output=True, check=False,
        )
        terminal_text = (
            "# terminal-qc command\n" + " ".join(terminal_command)
            + "\n# terminal-qc stdout\n" + terminal.stdout
            + "# terminal-qc stderr\n" + terminal.stderr
        )
        combined += terminal_text
        if (
            terminal.returncode != 0
            or "TRACK_B_LAVA_V2_TERMINAL_QC_VERIFIED" not in terminal.stdout
        ):
            raise ExecutionError(f"terminal QC attestation failed validation:\n{combined}")
    return values, combined


def _marker_agreement(worker: dict[str, str], semantic: dict[str, str]) -> None:
    keys = {
        "phase", "index", "locus", "chromosome", "pair", "qc", "construction_complete",
        "source_fingerprint", "continuation_fingerprint",
    }
    keys.add("n_snps")
    for key in keys:
        if worker[key] != semantic[key]:
            raise ExecutionError(f"worker and semantic markers differ: {key}")


def validate_bundle(
    phase: str, index: int | None, source_fingerprint: str,
    continuation_fingerprint: str, *, revalidate_semantics: bool = True,
    reconcile_active: bool = True,
) -> dict[str, Any] | None:
    bundle = bundle_path(phase, index, continuation_fingerprint)
    if not bundle.exists() and not bundle.is_symlink():
        return None
    target = _safe_attempt_target(bundle, continuation_fingerprint)
    required_files = [target / name for name in ("READY", "receipt.json", "worker.log", "result.rds", "semantic_validation.txt")]
    for path in required_files:
        if not path.is_file() or path.stat().st_size <= 0:
            raise ExecutionError(f"continuation checkpoint bundle is incomplete: {bundle}/{path.name}")
    try:
        ready = json.loads((target / "READY").read_text(encoding="utf-8"))
        receipt = json.loads((target / "receipt.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ExecutionError(f"continuation checkpoint metadata is unreadable: {bundle}") from error
    expected_ready = {
        "schema_version": READY_SCHEMA, "state": "READY",
        "receipt_sha256": sha256(target / "receipt.json"),
    }
    if ready != expected_ready:
        raise ExecutionError(f"continuation READY marker differs from receipt: {bundle}")
    analysis = {
        "aggregate-discovery": "LAVA_DISCOVERY_AGGREGATION_V2",
        "conditional": "LAVA_CONDITIONAL_V2", "finalize": "LAVA_FINALIZATION_V2",
        "terminal-qc": "LAVA_TERMINAL_QC_V2",
    }[phase]
    expected = {
        "schema_version": RECEIPT_SCHEMA, "phase": phase, "locus_index": index,
        "analysis": analysis,
        "requested_phase": "finalize" if phase == "terminal-qc" else phase,
        "source_discovery_fingerprint": source_fingerprint,
        "continuation_execution_fingerprint": continuation_fingerprint,
        "exit_status": 78 if phase == "terminal-qc" else 0,
        "bundle_path": _relative(bundle),
        "lineage_sha256": sha256(CONTRACT.lineage_path(continuation_fingerprint)),
        "source_family_lock_sha256": sha256(CONTRACT.source_lock_path(source_fingerprint)),
        "semantic_validator_sha256": sha256(CHECKPOINT_VALIDATOR),
        "terminal_qc_validator_sha256": sha256(RESULT_VALIDATOR) if phase == "terminal-qc" else None,
    }
    for key, value in expected.items():
        if receipt.get(key) != value:
            raise ExecutionError(f"continuation receipt identity drifted: {bundle} ({key})")
    required_receipt_fields = {
        "schema_version", "analysis", "phase", "requested_phase", "locus_index",
        "source_discovery_fingerprint", "continuation_execution_fingerprint", "marker",
        "command", "runtime_sec", "peak_rss_bytes", "exit_status", "bundle_path",
        "log_bytes", "log_sha256", "semantic_validator_sha256",
        "terminal_qc_validator_sha256", "lineage_sha256", "source_family_lock_sha256",
        "active_inputs", "semantic_artifact_binding", "artifacts",
    }
    if set(receipt) != required_receipt_fields:
        raise ExecutionError(f"continuation receipt schema drifted: {bundle}")
    marker = receipt.get("marker")
    if not isinstance(marker, dict):
        raise ExecutionError(f"continuation receipt lacks marker: {bundle}")
    validate_marker(
        marker, phase, index, target / "result.rds", source_fingerprint, continuation_fingerprint,
    )
    if (
        receipt.get("log_bytes") != (target / "worker.log").stat().st_size
        or receipt.get("log_sha256") != sha256(target / "worker.log")
    ):
        raise ExecutionError(f"continuation worker log differs from receipt: {bundle}")
    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise ExecutionError(f"continuation receipt lacks artifacts: {bundle}")
    recorded: set[str] = set()
    recorded_artifacts: dict[str, dict[str, object]] = {}
    for item in artifacts:
        if not isinstance(item, dict) or set(item) != {"path", "bytes", "sha256"}:
            raise ExecutionError(f"malformed continuation artifact receipt: {bundle}")
        name = str(item["path"])
        if name in recorded or "/" in name or name.startswith("."):
            raise ExecutionError(f"unsafe continuation artifact name: {bundle}/{name}")
        if file_identity(target / name) != (int(item["bytes"]), str(item["sha256"])):
            raise ExecutionError(f"continuation artifact differs from receipt: {bundle}/{name}")
        recorded.add(name)
        recorded_artifacts[name] = item
    live = {path.name for path in target.iterdir() if path.is_file()}
    if live != recorded | {"READY", "receipt.json", "worker.log"}:
        raise ExecutionError(f"continuation bundle contains unreceipted files: {bundle}")
    semantic_binding = receipt.get("semantic_artifact_binding")
    if not isinstance(semantic_binding, list):
        raise ExecutionError(f"continuation receipt lacks semantic-artifact binding: {bundle}")
    semantic_snapshot = semantic_snapshot_from_binding(semantic_binding, phase)
    for item in semantic_snapshot:
        name = str(item["path"])
        expected_artifact = {
            "path": name, "bytes": item["bytes"], "sha256": item["sha256"],
        }
        if (
            recorded_artifacts.get(name) != expected_artifact
            or stable_stat(target / name) != item["stable_stat"]
        ):
            raise ExecutionError(
                f"semantic-artifact publication identity drifted: {bundle}/{name}"
            )
    if revalidate_semantics:
        semantic, text = _semantic_validation(
            phase, index, source_fingerprint, continuation_fingerprint, target / "result.rds",
        )
        if (target / "semantic_validation.txt").read_text(encoding="utf-8") != text:
            raise ExecutionError(f"continuation semantic attestation drifted: {bundle}")
        _marker_agreement(marker, semantic)
        require_semantic_artifacts_unchanged(
            target, phase, semantic_snapshot, publication_metadata_allowed=True,
        )
    active = receipt.get("active_inputs")
    if not isinstance(active, list) or not active:
        raise ExecutionError(f"continuation receipt lacks active-input family: {bundle}")
    observed_specs = []
    for item in active:
        required = {"path", "bytes", "sha256", "role", "pre_stat", "post_stat"}
        if not isinstance(item, dict) or set(item) != required:
            raise ExecutionError(f"malformed continuation active-input receipt: {bundle}")
        path = _safe_relative(str(item["path"]))
        current = stable_stat(path)
        if current != item["pre_stat"] or current != item["post_stat"] or sha256(path) != item["sha256"]:
            raise ExecutionError(f"continuation active input differs from receipt: {item['path']}")
        observed_specs.append({key: item[key] for key in ("path", "bytes", "sha256", "role")})
    if reconcile_active:
        expected_specs = active_specs(
            phase, index, source_fingerprint, continuation_fingerprint,
        )
        if observed_specs != expected_specs:
            raise ExecutionError(f"continuation active-input family drifted: {bundle}")
    return receipt


def SOURCE_CHECKPOINT_ROOT(source_fingerprint: str) -> Path:
    return CONTRACT.SOURCE_CHECKPOINT_ROOT / source_fingerprint


def _tsv(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def benchmark_paths(continuation_fingerprint: str) -> tuple[Path, Path, Path]:
    directory = BENCHMARK_ROOT / continuation_fingerprint
    return directory / "RAM_BENCHMARK.tsv", directory / "RAM_BENCHMARK.namespace.json", directory / "RAM_BENCHMARK.provenance.json"


def ensure_benchmark_namespace(source_fingerprint: str, continuation_fingerprint: str) -> None:
    table, namespace, _ = benchmark_paths(continuation_fingerprint)
    expected = {
        "schema_version": "track-b-lava-continuation-ram-namespace.1",
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source_fingerprint,
        "continuation_execution_fingerprint": continuation_fingerprint,
    }
    if namespace.exists():
        if json.loads(namespace.read_text(encoding="utf-8")) != expected:
            raise ExecutionError("continuation RAM namespace identity drifted")
    else:
        publish_no_replace(namespace, json.dumps(expected, indent=2, sort_keys=True).encode() + b"\n")
    if not table.exists():
        publish_no_replace(table, _tsv(BENCHMARK_FIELDS, []).encode())


def update_benchmark(receipt: dict[str, Any]) -> None:
    source = str(receipt["source_discovery_fingerprint"])
    continuation = str(receipt["continuation_execution_fingerprint"])
    ensure_benchmark_namespace(source, continuation)
    table, _, provenance = benchmark_paths(continuation)
    with table.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    marker = receipt["marker"]
    result_sha = next(item["sha256"] for item in receipt["artifacts"] if item["path"] == "result.rds")
    row: dict[str, object] = {
        "source_discovery_fingerprint": source,
        "continuation_execution_fingerprint": continuation,
        "analysis": receipt["analysis"], "phase": receipt["phase"],
        "locus_index": "NA" if receipt["locus_index"] is None else receipt["locus_index"],
        "locus": marker["locus"], "chromosome": marker["chromosome"], "pair": marker["pair"],
        "qc": marker["qc"], "peak_ram_gb": f"{int(receipt['peak_rss_bytes']) / 1024**3:.6f}",
        "runtime_sec": f"{float(receipt['runtime_sec']):.3f}", "exit_status": receipt["exit_status"],
        "output_sha256": result_sha,
    }
    key = (str(row["phase"]), str(row["locus_index"]), str(row["output_sha256"]))
    rows = [item for item in rows if (item["phase"], item["locus_index"], item["output_sha256"]) != key]
    rows.append(row)
    rows.sort(key=lambda item: (item["phase"], int(item["locus_index"]) if item["locus_index"] != "NA" else 0))
    atomic_text(table, _tsv(BENCHMARK_FIELDS, rows))
    payload = {
        "schema_version": "track-b-lava-continuation-ram-provenance.1",
        "source_discovery_fingerprint": source,
        "continuation_execution_fingerprint": continuation,
        "namespace_sha256": sha256(benchmark_paths(continuation)[1]),
        "row_count": len(rows), "table_sha256": sha256(table),
    }
    atomic_text(provenance, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def rss_bytes(value: int) -> int:
    return int(value if sys.platform == "darwin" else value * 1024)


def internal_measure(args: argparse.Namespace) -> int:
    # This hidden entry point is callable outside the parent supervisor.  It
    # therefore performs its own live contract/lineage check before R can read
    # or write inside a continuation namespace.
    quick_contract(args.source_fingerprint, args.continuation_fingerprint)
    command = [
        str(R_SCRIPT), str(RUNNER), "--phase", args.phase,
        "--source-fingerprint", args.source_fingerprint,
        "--continuation-fingerprint", args.continuation_fingerprint,
        "--source-checkpoint-root", str(SOURCE_CHECKPOINT_ROOT(args.source_fingerprint)),
        "--continuation-checkpoint-root", str(run_root(args.continuation_fingerprint)),
        "--worker-output", args.worker_output,
    ]
    if args.locus_index is not None:
        command.extend(["--locus-index", str(args.locus_index)])
    started = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    elapsed = time.perf_counter() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    Path(args.internal_result).write_text(json.dumps({
        "command": command, "returncode": result.returncode, "runtime_sec": elapsed,
        "peak_rss_bytes": rss_bytes(int(usage.ru_maxrss)),
        "stdout": result.stdout, "stderr": result.stderr,
    }), encoding="utf-8")
    return 0


def _publish_receipt(
    attempt: Path, phase: str, index: int | None, source_fingerprint: str,
    continuation_fingerprint: str, marker: dict[str, str], metrics: dict[str, Any],
    active_inputs: list[dict[str, object]], semantic_binding: list[dict[str, object]],
    semantic_text: str, log: str,
) -> dict[str, Any]:
    semantic_artifacts = semantic_snapshot_from_binding(semantic_binding, phase)
    require_semantic_artifacts_unchanged(
        attempt, phase, semantic_artifacts,
    )
    publish_no_replace(attempt / "worker.log", log.encode())
    publish_no_replace(attempt / "semantic_validation.txt", semantic_text.encode())
    artifacts = []
    for path in sorted(
        (item for item in attempt.iterdir() if item.is_file() and item.name not in {"worker.log", "receipt.json", "READY"}),
        key=lambda item: item.name,
    ):
        size, digest = file_identity(path)
        artifacts.append({"path": path.name, "bytes": size, "sha256": digest})
    artifact_by_name = {str(item["path"]): item for item in artifacts}
    expected_semantic_receipt = {
        str(item["path"]): {
            "path": item["path"], "bytes": item["bytes"], "sha256": item["sha256"],
        }
        for item in semantic_artifacts
    }
    if any(
        artifact_by_name.get(name) != expected
        for name, expected in expected_semantic_receipt.items()
    ):
        raise ExecutionError(
            f"receipt artifacts differ from the semantically validated bytes: {phase}"
        )
    analysis = {
        "aggregate-discovery": "LAVA_DISCOVERY_AGGREGATION_V2",
        "conditional": "LAVA_CONDITIONAL_V2", "finalize": "LAVA_FINALIZATION_V2",
        "terminal-qc": "LAVA_TERMINAL_QC_V2",
    }[phase]
    receipt: dict[str, Any] = {
        "schema_version": RECEIPT_SCHEMA, "analysis": analysis, "phase": phase,
        "requested_phase": "finalize" if phase == "terminal-qc" else phase,
        "locus_index": index,
        "source_discovery_fingerprint": source_fingerprint,
        "continuation_execution_fingerprint": continuation_fingerprint,
        "marker": marker, "command": metrics["command"],
        "runtime_sec": float(metrics["runtime_sec"]),
        "peak_rss_bytes": int(metrics["peak_rss_bytes"]),
        "exit_status": 78 if phase == "terminal-qc" else 0,
        "bundle_path": _relative(bundle_path(phase, index, continuation_fingerprint)),
        "log_bytes": (attempt / "worker.log").stat().st_size,
        "log_sha256": sha256(attempt / "worker.log"),
        "semantic_validator_sha256": sha256(CHECKPOINT_VALIDATOR),
        "terminal_qc_validator_sha256": sha256(RESULT_VALIDATOR) if phase == "terminal-qc" else None,
        "lineage_sha256": sha256(CONTRACT.lineage_path(continuation_fingerprint)),
        "source_family_lock_sha256": sha256(CONTRACT.source_lock_path(source_fingerprint)),
        "active_inputs": active_inputs,
        "semantic_artifact_binding": semantic_binding,
        "artifacts": artifacts,
    }
    publish_no_replace(attempt / "receipt.json", json.dumps(receipt, indent=2, sort_keys=True).encode() + b"\n")
    require_semantic_artifacts_unchanged(
        attempt, phase, semantic_artifacts, publication_metadata_allowed=True,
    )
    ready = {"schema_version": READY_SCHEMA, "receipt_sha256": sha256(attempt / "receipt.json"), "state": "READY"}
    publish_no_replace(attempt / "READY", json.dumps(ready, sort_keys=True).encode() + b"\n")
    require_semantic_artifacts_unchanged(
        attempt, phase, semantic_artifacts, publication_metadata_allowed=True,
    )
    fsync_directory(attempt)
    _publish_ready_bundle(attempt, phase, index, continuation_fingerprint)
    published = validate_bundle(
        phase, index, source_fingerprint, continuation_fingerprint,
        revalidate_semantics=False, reconcile_active=False,
    )
    if published is None:
        raise ExecutionError("published continuation checkpoint disappeared")
    if published.get("active_inputs") != active_inputs:
        raise ExecutionError("published continuation active-input receipt differs from capture")
    if published.get("semantic_artifact_binding") != semantic_binding:
        raise ExecutionError("published semantic-artifact binding differs from validation")
    if (attempt / "semantic_validation.txt").read_text(encoding="utf-8") != semantic_text:
        raise ExecutionError("published continuation semantic attestation differs from validation")
    published_artifacts = {
        str(item["path"]): item for item in published.get("artifacts", [])
    }
    if any(
        published_artifacts.get(name) != expected
        for name, expected in expected_semantic_receipt.items()
    ):
        raise ExecutionError(
            f"published receipt differs from the semantically validated bytes: {phase}"
        )
    require_semantic_artifacts_unchanged(
        attempt, phase, semantic_artifacts, publication_metadata_allowed=True,
    )
    update_benchmark(published)
    return published


def run_measured(
    phase: str, index: int | None, source_fingerprint: str, continuation_fingerprint: str,
) -> dict[str, Any]:
    if phase not in EXECUTION_PHASES:
        raise ExecutionError("V2 supervisor cannot execute discovery or synthetic phases")
    validate_phase_prerequisites(phase, source_fingerprint, continuation_fingerprint)
    existing = validate_bundle(phase, index, source_fingerprint, continuation_fingerprint)
    if existing is not None:
        update_benchmark(existing)
        return existing
    if phase == "finalize":
        terminal = validate_bundle("terminal-qc", None, source_fingerprint, continuation_fingerprint)
        if terminal is not None:
            update_benchmark(terminal)
            raise TerminalQCError("Track B LAVA continuation reached immutable TERMINAL_FAILED_QC")
    specs = active_specs(phase, index, source_fingerprint, continuation_fingerprint)
    captured = capture_active_inputs(specs)
    attempt = _attempt_directory(phase, index, continuation_fingerprint)
    output = attempt / "result.rds"
    with tempfile.NamedTemporaryFile(prefix="track-b-lava-v2-measure-", suffix=".json", delete=False) as handle:
        metrics_path = Path(handle.name)
    command = [
        sys.executable, str(Path(__file__).resolve()), "--internal-measure", "--phase", phase,
        "--source-fingerprint", source_fingerprint,
        "--continuation-fingerprint", continuation_fingerprint,
        "--worker-output", str(output), "--internal-result", str(metrics_path),
    ]
    if index is not None:
        command.extend(["--locus-index", str(index)])
    try:
        helper = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
        if helper.returncode != 0 or not metrics_path.is_file():
            raise ExecutionError(f"continuation RSS helper failed: {helper.stdout}{helper.stderr}")
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
    finally:
        metrics_path.unlink(missing_ok=True)
    log = (
        "# command\n" + " ".join(str(value) for value in metrics["command"])
        + "\n# stdout\n" + str(metrics["stdout"]) + "\n# stderr\n" + str(metrics["stderr"])
    )
    quick_contract(source_fingerprint, continuation_fingerprint)
    active_inputs = complete_active_inputs(captured)
    returncode = int(metrics["returncode"])
    if returncode == 78 and phase == "finalize":
        marker = parse_marker(log)
        validate_marker(marker, "terminal-qc", None, output, source_fingerprint, continuation_fingerprint)
        if not (attempt / "terminal_qc.json").is_file():
            raise ExecutionError("terminal QC worker omitted its machine-readable attestation")
        semantic_before = capture_semantic_artifacts(attempt, "terminal-qc")
        semantic, semantic_text = _semantic_validation(
            "terminal-qc", None, source_fingerprint, continuation_fingerprint, output,
        )
        semantic_after = capture_semantic_artifacts(attempt, "terminal-qc")
        semantic_binding = bind_semantic_validation(
            semantic_before, semantic_after, "terminal-qc",
        )
        _marker_agreement(marker, semantic)
        _publish_receipt(
            attempt, "terminal-qc", None, source_fingerprint, continuation_fingerprint,
            marker, metrics, active_inputs, semantic_binding, semantic_text, log,
        )
        raise TerminalQCError("Track B LAVA continuation reached immutable TERMINAL_FAILED_QC")
    if returncode != 0:
        failure_dir = run_root(continuation_fingerprint) / "failed_attempts"
        failure_dir.mkdir(parents=True, exist_ok=True)
        failure = failure_dir / f"{phase}_{index if index is not None else 'all'}_{time.time_ns()}.log"
        publish_no_replace(failure, log.encode())
        raise ExecutionError(f"continuation {phase} worker failed with exit status {returncode}: {failure}")
    if not output.is_file() or output.stat().st_size <= 0:
        raise ExecutionError("successful continuation worker omitted result.rds")
    marker = parse_marker(log)
    validate_marker(marker, phase, index, output, source_fingerprint, continuation_fingerprint)
    semantic_before = capture_semantic_artifacts(attempt, phase)
    semantic, semantic_text = _semantic_validation(
        phase, index, source_fingerprint, continuation_fingerprint, output,
    )
    semantic_after = capture_semantic_artifacts(attempt, phase)
    semantic_binding = bind_semantic_validation(
        semantic_before, semantic_after, phase,
    )
    _marker_agreement(marker, semantic)
    return _publish_receipt(
        attempt, phase, index, source_fingerprint, continuation_fingerprint,
        marker, metrics, active_inputs, semantic_binding, semantic_text, log,
    )


def publish_results(source_fingerprint: str, continuation_fingerprint: str) -> None:
    command = [
        sys.executable, str(RESULT_VALIDATOR), "--publish-staged-results",
        "--source-fingerprint", source_fingerprint,
        "--continuation-fingerprint", continuation_fingerprint,
    ]
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=False)
    if result.returncode != 0:
        raise ExecutionError("validated continuation could not publish canonical results:\n" + result.stdout + result.stderr)


def run_all(source_fingerprint: str, continuation_fingerprint: str) -> None:
    run_measured("aggregate-discovery", None, source_fingerprint, continuation_fingerprint)
    for index in range(1, 2496):
        run_measured("conditional", index, source_fingerprint, continuation_fingerprint)
        if index % 100 == 0 or index == 2495:
            print(f"LAVA_CONTINUATION_CONDITIONAL_PROGRESS completed_through={index} total=2495")
    run_measured("finalize", None, source_fingerprint, continuation_fingerprint)
    publish_results(source_fingerprint, continuation_fingerprint)
    print(
        "TRACK_B_LAVA_CONTINUATION_COMPLETE "
        f"source={source_fingerprint} continuation={continuation_fingerprint} discovery_recomputed=FALSE"
    )


def preflight(source_fingerprint: str) -> str:
    try:
        _, continuation = CONTRACT.preflight(source_fingerprint)
    except CONTRACT.ContinuationError as error:
        raise ExecutionError(str(error)) from error
    ensure_benchmark_namespace(source_fingerprint, continuation)
    return continuation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--verify", action="store_true")
    mode.add_argument("--internal-measure", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--source-fingerprint", default=CONTRACT.PINNED_SOURCE_FINGERPRINT)
    parser.add_argument("--continuation-fingerprint")
    parser.add_argument("--phase", choices=sorted(EXECUTION_PHASES))
    parser.add_argument("--locus-index", type=int)
    parser.add_argument("--worker-output", help=argparse.SUPPRESS)
    parser.add_argument("--internal-result", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.internal_measure:
        required = [args.phase, args.source_fingerprint, args.continuation_fingerprint, args.worker_output, args.internal_result]
        if any(value is None for value in required):
            parser.error("internal measurement requires complete continuation identity and paths")
        raise SystemExit(internal_measure(args))
    continuation = preflight(args.source_fingerprint)
    if args.continuation_fingerprint and args.continuation_fingerprint != continuation:
        raise ExecutionError("supplied continuation fingerprint differs from preflight")
    if args.preflight:
        print(
            "TRACK_B_LAVA_CONTINUATION_PREFLIGHT_PASS "
            f"source={args.source_fingerprint} continuation={continuation}"
        )
    elif args.verify:
        aggregate = validate_bundle("aggregate-discovery", None, args.source_fingerprint, continuation)
        terminal = validate_bundle("terminal-qc", None, args.source_fingerprint, continuation)
        conditional_count = sum(
            validate_bundle(
                "conditional", index, args.source_fingerprint, continuation,
                revalidate_semantics=False, reconcile_active=False,
            ) is not None
            for index in range(1, 2496)
        )
        state = "TERMINAL_FAILED_QC" if terminal else "COMPLETE" if validate_bundle(
            "finalize", None, args.source_fingerprint, continuation,
        ) else "IN_PROGRESS"
        print(
            "TRACK_B_LAVA_CONTINUATION_STATUS "
            f"state={state} aggregate={int(aggregate is not None)} conditional={conditional_count}/2495"
        )
    else:
        run_all(args.source_fingerprint, continuation)


if __name__ == "__main__":
    try:
        main()
    except TerminalQCError as error:
        print(f"TRACK_B_LAVA_CONTINUATION_TERMINAL_FAILED_QC message={error}", file=sys.stderr)
        raise SystemExit(78) from error
    except (ExecutionError, CONTRACT.ContinuationError) as error:
        raise SystemExit(f"ERROR: {error}") from error
