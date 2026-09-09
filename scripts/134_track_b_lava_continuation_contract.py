#!/usr/bin/env python3
"""Seal and verify an additive Track B LAVA post-discovery continuation.

The completed V1 discovery family remains owned by its original execution
fingerprint.  Post-discovery code has an independent fingerprint and may only
consume that family through a no-replace, content-addressed source lock.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import stat
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
LEGACY_CONTRACT_PATH = ROOT / "scripts/119_track_b_lava_contract.py"
LEGACY_SUPERVISOR_PATH = ROOT / "scripts/131_run_track_b_lava_sequential.py"
V2_RUNNER = ROOT / "scripts/135_run_track_b_lava_postdiscovery_v2.R"
V2_SUPERVISOR = ROOT / "scripts/136_run_track_b_lava_continuation.py"
V2_CHECKPOINT_VALIDATOR = ROOT / "scripts/137_validate_track_b_lava_checkpoint_v2.R"
V2_RESULT_VALIDATOR = ROOT / "scripts/138_validate_track_b_lava_results_v2.py"
SUPERSESSION_VALIDATOR = ROOT / "scripts/144_freeze_track_b_lava_supersession.py"

SOURCE_CHECKPOINT_ROOT = ROOT / "results/track_b/checkpoints/lava"
CONTINUATION_CHECKPOINT_ROOT = ROOT / "results/track_b/checkpoints/lava_continuations"
SOURCE_LOCK_ROOT = CONTINUATION_CHECKPOINT_ROOT / "source_locks"
RUN_ROOT = CONTINUATION_CHECKPOINT_ROOT / "runs"

SOURCE_LOCK_SCHEMA = "track-b-lava-discovery-family.1"
LINEAGE_SCHEMA = "track-b-lava-continuation-lineage.2"
CONTINUATION_CONTRACT_SCHEMA = "track-b-lava-continuation-contract.1"
PINNED_SOURCE_FINGERPRINT = "caeb6b5a1188560f27609cad77801715e0bedfb998660a531c5491afad1177c7"
PREDECESSOR_CONTINUATION_FINGERPRINT = "09cd41888ba405ae2b54b8584c5bc15c7d0213e81e7c1cd8f661136bca34a152"
PREDECESSOR_SUPERSESSION_LOCK = (
    ROOT / "results/track_b/lava_continuations/supersessions/final_locks"
    / f"{PREDECESSOR_CONTINUATION_FINGERPRINT}.lock.json"
)
PREDECESSOR_READY_THROUGH = 22
PREDECESSOR_STOP_REASON = "PERFORMANCE_ONLY_REDUNDANT_SUPERVISOR_VALIDATION_OVERHEAD"
PREDECESSOR_PERFORMANCE_RATIONALE = (
    "OBSERVED_APPROXIMATELY_20_SECOND_CONDITIONAL_START_CADENCE_WITH_2_TO_3_SECOND_"
    "WORKERS;OPTIMIZE_DUPLICATE_CONTRACT_LINEAGE_AGGREGATE_SOURCE_LOOKUP_AND_POST_"
    "PUBLICATION_VALIDATION_WITHOUT_CHANGING_SCIENTIFIC_GATES"
)
PREDECESSOR_STOP_BOUNDARY = (
    "CLEAN_INTERRUPT_BETWEEN_LOCI_DURING_NEXT_AGGREGATE_ACTIVE_INPUT_HASH_VERIFICATION"
)
FINGERPRINT_RE = re.compile(r"^[0-9a-f]{64}$")

# These are the exact live V1 bytes from which caeb6b5a... was computed.  They
# deliberately remain outside the V2 fingerprint and must never be edited as
# part of the continuation.
LEGACY_CODE_SHA256 = {
    "scripts/119_track_b_lava_contract.py": "9f2f8e1708a040bcfb70cf31f08dee7ab4c3bcd914fd32a1347dde24cc5abd52",
    "scripts/120_run_track_b_lava.R": "ce40e9c9c32908cc7b396d74a17c12d49ec7b65eaa8b6eaf0ef78ece990237c3",
    "scripts/121_validate_track_b_lava.py": "0524c56177d64412ee5563d8e6702957315822f7933b9e71a0ac63860c06f029",
    "scripts/131_run_track_b_lava_sequential.py": "7c7001036d14fdc5842d153900f8cb0b114abb219038089a146dc1c37445c889",
    "scripts/132_validate_track_b_lava_checkpoint.R": "c516a546d4300d7093721d1141282bc7b61cdc5fb0d35f649ef31b69e9b128f1",
    "scripts/133_validate_track_b_lava_runtime.R": "05baf8d3abf8766afd680ecd8c44bc48dd94a965d31281dc25a8f7c79c748c0e",
    "scripts/129_prepare_track_b_lava_chromosome_inputs.py": "67740a3eb2ff452a3c9c3d5785228f6770d98dc6d9aadc057e5a913a782ac8c8",
    "scripts/30_setup_lava.sh": "a1070fb11a81c6cc0d73e14388a971fa851759c05858ca8f83542bb8975ea2bb",
    "scripts/lava_contract.py": "8b497a53d376a2ebfec9128b2e2b8f430f6dd1f7f0fd258ffd6b838ce590cfd7",
}


class ContinuationError(RuntimeError):
    """A fail-closed continuation-contract violation."""


def _load_module(name: str, path: Path):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise ContinuationError(f"could not load required module: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def legacy_contract():
    return _load_module("track_b_lava_legacy_contract", LEGACY_CONTRACT_PATH)


def legacy_supervisor():
    return _load_module("track_b_lava_legacy_supervisor", LEGACY_SUPERVISOR_PATH)


def supersession_validator():
    return _load_module("track_b_lava_supersession_validator", SUPERSESSION_VALIDATOR)


def _canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def file_identity(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size <= 0:
                raise ContinuationError(f"artifact is not a real non-empty file: {path}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        raise ContinuationError(f"could not hash artifact {path}: {error}") from error
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before) != identity(after) or identity(after) != identity(current):
        raise ContinuationError(f"artifact changed while hashing: {path}")
    return int(after.st_size), digest.hexdigest()


def sha256(path: Path) -> str:
    return file_identity(path)[1]


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def require_real_repository_path(path: Path, label: str) -> Path:
    """Reject symlinks in an existing repository-relative path prefix.

    Continuation namespaces are writable, so checking only ``resolve()`` against
    a resolved parent would allow a pre-positioned namespace symlink to redefine
    that parent outside the repository.  Missing suffixes are permitted so the
    caller can create a new fingerprint namespace, but every existing component
    below the repository root must be a real directory (or the final real file).
    """

    try:
        relative = path.relative_to(ROOT)
    except ValueError as error:
        raise ContinuationError(f"{label} escapes repository: {path}") from error
    current = ROOT
    for offset, part in enumerate(relative.parts):
        current /= part
        try:
            current_stat = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise ContinuationError(f"could not inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(current_stat.st_mode):
            raise ContinuationError(f"{label} contains a symbolic link: {current}")
        if offset < len(relative.parts) - 1 and not stat.S_ISDIR(current_stat.st_mode):
            raise ContinuationError(f"{label} has a non-directory ancestor: {current}")
    try:
        path.resolve(strict=False).relative_to(ROOT.resolve(strict=True))
    except (OSError, RuntimeError, ValueError) as error:
        raise ContinuationError(f"{label} escapes repository: {path}") from error
    return path


def publish_no_replace(path: Path, content: bytes) -> None:
    if path == ROOT or ROOT in path.parents:
        require_real_repository_path(path.parent, "publication parent")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path == ROOT or ROOT in path.parents:
        require_real_repository_path(path.parent, "publication parent")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise ContinuationError(f"refusing to replace existing artifact: {path}") from error
        fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def _relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError as error:
        raise ContinuationError(f"path escapes repository: {path}") from error


def _safe_relative(value: str, label: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ContinuationError(f"unsafe {label} path: {value}")
    path = ROOT / relative
    try:
        path.resolve(strict=False).relative_to(ROOT.resolve())
    except ValueError as error:
        raise ContinuationError(f"{label} path escapes repository: {value}") from error
    return path


def _identity_record(path: Path) -> dict[str, object]:
    size, digest = file_identity(path)
    return {"path": _relative(path), "bytes": size, "sha256": digest}


def _require_source_fingerprint(value: str) -> str:
    if value != PINNED_SOURCE_FINGERPRINT or not FINGERPRINT_RE.fullmatch(value):
        raise ContinuationError("unsupported source discovery fingerprint")
    return value


def predecessor_supersession_lock_path() -> Path:
    require_real_repository_path(
        PREDECESSOR_SUPERSESSION_LOCK, "predecessor supersession-lock path"
    )
    return PREDECESSOR_SUPERSESSION_LOCK


def validate_predecessor_supersession_lock() -> dict[str, object]:
    """Verify the no-replace predecessor seal before deriving a successor."""

    validator = supersession_validator()
    try:
        payload = validator.verify(PREDECESSOR_CONTINUATION_FINGERPRINT)
    except validator.SupersessionError as error:
        raise ContinuationError(f"predecessor supersession lock is invalid: {error}") from error
    snapshot = payload.get("verified_snapshot") if isinstance(payload, dict) else None
    family = snapshot.get("conditional_partial_family") if isinstance(snapshot, dict) else None
    aggregate = snapshot.get("aggregate") if isinstance(snapshot, dict) else None
    expected_binding = {
        "mode": "SUCCESSOR_CONTRACT_BINDS_THIS_LOCK",
        "successor_fingerprint": None,
    }
    if (
        payload.get("schema_version") != validator.DIRECT_LOCK_SCHEMA
        or payload.get("state") != "FINAL_LOCKED"
        or payload.get("analysis_id") != "track-b-v1.0-local"
        or payload.get("source_discovery_fingerprint") != PINNED_SOURCE_FINGERPRINT
        or payload.get("superseded_continuation_fingerprint")
        != PREDECESSOR_CONTINUATION_FINGERPRINT
        or payload.get("status") != "SUPERSEDED_PARTIAL_PERFORMANCE_OPTIMIZATION"
        or payload.get("stop_reason") != PREDECESSOR_STOP_REASON
        or payload.get("performance_rationale") != PREDECESSOR_PERFORMANCE_RATIONALE
        or payload.get("stop_boundary") != PREDECESSOR_STOP_BOUNDARY
        or payload.get("scientific_failure") is not False
        or payload.get("scientific_results_relabelled") is not False
        or payload.get("successor_binding") != expected_binding
        or not isinstance(family, dict)
        or family.get("complete") is not False
        or family.get("ready_count") != PREDECESSOR_READY_THROUGH
        or family.get("ready_indices") != list(range(1, PREDECESSOR_READY_THROUGH + 1))
        or family.get("no_other_ready_attempts") is not True
        or not isinstance(aggregate, dict)
        or aggregate.get("candidate_family") != "HEADER_ONLY_ZERO_CANDIDATES"
        or snapshot.get("failed_attempt_count") != 0
        or snapshot.get("finalize_present") is not False
        or snapshot.get("terminal_qc_present") is not False
    ):
        raise ContinuationError("predecessor supersession-lock scope drifted")
    path = predecessor_supersession_lock_path()
    if validator.lock_path(PREDECESSOR_CONTINUATION_FINGERPRINT) != path:
        raise ContinuationError("predecessor supersession-lock canonical path drifted")
    return _identity_record(path)


def source_lock_path(source_fingerprint: str) -> Path:
    path = SOURCE_LOCK_ROOT / f"{_require_source_fingerprint(source_fingerprint)}.json"
    require_real_repository_path(path, "source-lock namespace")
    return path


def continuation_run_root(continuation_fingerprint: str) -> Path:
    if not FINGERPRINT_RE.fullmatch(continuation_fingerprint):
        raise ContinuationError("invalid continuation fingerprint")
    path = RUN_ROOT / continuation_fingerprint
    require_real_repository_path(path, "continuation run namespace")
    return path


def lineage_path(continuation_fingerprint: str) -> Path:
    return continuation_run_root(continuation_fingerprint) / "lineage.json"


def _validate_legacy_code() -> dict[str, dict[str, object]]:
    observed: dict[str, dict[str, object]] = {}
    for relative, expected in LEGACY_CODE_SHA256.items():
        path = ROOT / relative
        size, digest = file_identity(path)
        if digest != expected:
            raise ContinuationError(f"legacy V1 code drifted: {relative}")
        observed[relative] = {"bytes": size, "sha256": digest}
    return observed


def _require_live_source_fingerprint(source_fingerprint: str) -> None:
    _validate_legacy_code()
    observed = legacy_contract().run_fingerprint()
    if observed != source_fingerprint:
        raise ContinuationError(
            "live V1 inputs/runtime no longer reproduce the source discovery fingerprint"
        )


def _source_entry(index: int, supervisor: Any, source_fingerprint: str, *, deep: bool) -> dict[str, object]:
    receipt = supervisor.validate_bundle(
        "discovery", index, source_fingerprint,
        _reconcile_active=True, _revalidate_semantic=deep,
    )
    if receipt is None:
        raise ContinuationError(f"missing source discovery checkpoint: {index}")
    bundle = supervisor.bundle_path("discovery", index, source_fingerprint)
    target = supervisor._safe_attempt_target(bundle)
    artifact_by_name = {str(item["path"]): item for item in receipt["artifacts"]}
    if "result.rds" not in artifact_by_name or "semantic_validation.txt" not in artifact_by_name:
        raise ContinuationError(f"source checkpoint lacks required artifacts: {index}")
    result_record = _identity_record(target / "result.rds")
    if (
        result_record["bytes"] != artifact_by_name["result.rds"]["bytes"]
        or result_record["sha256"] != artifact_by_name["result.rds"]["sha256"]
    ):
        raise ContinuationError(f"source result differs from receipt: {index}")
    return {
        "locus_index": index,
        "locus": str(receipt["marker"]["locus"]),
        "chromosome": str(receipt["marker"]["chromosome"]),
        "qc": str(receipt["marker"]["qc"]),
        "bundle_path": _relative(bundle),
        "attempt_path": _relative(target),
        "ready": _identity_record(target / "READY"),
        "receipt": _identity_record(target / "receipt.json"),
        "result": result_record,
        "log": _identity_record(target / "worker.log"),
        "semantic_validation": _identity_record(target / "semantic_validation.txt"),
    }


def build_source_lock(source_fingerprint: str, *, deep: bool = False) -> dict[str, object]:
    source_fingerprint = _require_source_fingerprint(source_fingerprint)
    _require_live_source_fingerprint(source_fingerprint)
    supervisor = legacy_supervisor()
    baseline = supervisor.validate_input_baseline(source_fingerprint)
    del baseline
    entries = [
        _source_entry(index, supervisor, source_fingerprint, deep=deep)
        for index in range(1, 2496)
    ]
    baseline_record = _identity_record(supervisor.input_baseline_path(source_fingerprint))
    payload: dict[str, object] = {
        "schema_version": SOURCE_LOCK_SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source_fingerprint,
        "compatibility_scope": {
            "source_phase": "discovery",
            "allowed_continuation_phases": ["aggregate-discovery", "conditional", "finalize"],
            "discovery_reexecution_allowed": False,
            "reason": "POST_DISCOVERY_EMPTY_CONDITIONAL_CANDIDATE_INTEGER_CONSTRUCTOR_COLLISION",
        },
        "legacy_code": _validate_legacy_code(),
        "active_input_baseline": baseline_record,
        "checkpoint_count": len(entries),
        "checkpoint_family_sha256": _digest(entries),
        "checkpoints": entries,
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    path = source_lock_path(source_fingerprint)
    if path.exists():
        observed = validate_source_lock(
            source_fingerprint, validate_bundles=True,
            revalidate_semantics=deep, require_live_fingerprint=True,
        )
        if observed != payload or path.read_bytes() != encoded:
            raise ContinuationError("existing source-family lock differs; overwrite is forbidden")
        return observed
    publish_no_replace(path, encoded)
    return validate_source_lock(
        source_fingerprint, validate_bundles=True,
        revalidate_semantics=deep, require_live_fingerprint=True,
    )


def _validate_record(record: object, label: str) -> tuple[Path, int, str]:
    if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
        raise ContinuationError(f"malformed {label} identity")
    path = _safe_relative(str(record["path"]), label)
    try:
        expected_size = int(record["bytes"])
    except (TypeError, ValueError) as error:
        raise ContinuationError(f"invalid {label} byte count") from error
    expected_digest = str(record["sha256"])
    if expected_size <= 0 or not FINGERPRINT_RE.fullmatch(expected_digest):
        raise ContinuationError(f"invalid {label} content identity")
    if file_identity(path) != (expected_size, expected_digest):
        raise ContinuationError(f"{label} differs from source-family lock")
    return path, expected_size, expected_digest


def _validate_record_shape(record: object, label: str) -> tuple[Path, int, str]:
    if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
        raise ContinuationError(f"malformed {label} identity")
    path = _safe_relative(str(record["path"]), label)
    try:
        expected_size = int(record["bytes"])
    except (TypeError, ValueError) as error:
        raise ContinuationError(f"invalid {label} byte count") from error
    expected_digest = str(record["sha256"])
    if expected_size <= 0 or not FINGERPRINT_RE.fullmatch(expected_digest):
        raise ContinuationError(f"invalid {label} content identity")
    return path, expected_size, expected_digest


def _require_canonical_source_entry(
    entry: dict[str, object], index: int, supervisor: Any,
    source_fingerprint: str, *, revalidate_semantics: bool,
) -> None:
    """Bind every lock field to the canonical READY bundle and its receipt."""

    expected = _source_entry(
        index, supervisor, source_fingerprint, deep=revalidate_semantics,
    )
    if entry != expected:
        raise ContinuationError(
            f"source checkpoint lock entry differs from canonical READY receipt: {index}"
        )


def validate_source_lock(
    source_fingerprint: str, *, validate_bundles: bool = True,
    revalidate_semantics: bool = False, require_live_fingerprint: bool = True,
) -> dict[str, object]:
    source_fingerprint = _require_source_fingerprint(source_fingerprint)
    if require_live_fingerprint:
        _require_live_source_fingerprint(source_fingerprint)
    path = source_lock_path(source_fingerprint)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContinuationError(f"source-family lock is unreadable: {error}") from error
    expected_top = {
        "schema_version", "analysis_id", "source_discovery_fingerprint", "compatibility_scope",
        "legacy_code", "active_input_baseline", "checkpoint_count",
        "checkpoint_family_sha256", "checkpoints",
    }
    if not isinstance(payload, dict) or set(payload) != expected_top:
        raise ContinuationError("source-family lock schema drifted")
    expected_scope = {
        "source_phase": "discovery",
        "allowed_continuation_phases": ["aggregate-discovery", "conditional", "finalize"],
        "discovery_reexecution_allowed": False,
        "reason": "POST_DISCOVERY_EMPTY_CONDITIONAL_CANDIDATE_INTEGER_CONSTRUCTOR_COLLISION",
    }
    if (
        payload["schema_version"] != SOURCE_LOCK_SCHEMA
        or payload["analysis_id"] != "track-b-v1.0-local"
        or payload["source_discovery_fingerprint"] != source_fingerprint
        or payload["compatibility_scope"] != expected_scope
        or payload["legacy_code"] != _validate_legacy_code()
    ):
        raise ContinuationError("source-family lock scope or V1 identity drifted")
    _validate_record(payload["active_input_baseline"], "active-input baseline")
    entries = payload["checkpoints"]
    if (
        not isinstance(entries, list) or len(entries) != 2495
        or payload["checkpoint_count"] != 2495
        or payload["checkpoint_family_sha256"] != _digest(entries)
    ):
        raise ContinuationError("source discovery checkpoint-family digest or count drifted")
    if [item.get("locus_index") if isinstance(item, dict) else None for item in entries] != list(range(1, 2496)):
        raise ContinuationError("source discovery indices are incomplete, duplicated, or unordered")
    loci = legacy_supervisor().load_loci()
    supervisor = legacy_supervisor() if validate_bundles else None
    required_entry = {
        "locus_index", "locus", "chromosome", "qc", "bundle_path", "attempt_path",
        "ready", "receipt", "result", "log", "semantic_validation",
    }
    for index, (entry, locus) in enumerate(zip(entries, loci, strict=True), start=1):
        if not isinstance(entry, dict) or set(entry) != required_entry:
            raise ContinuationError(f"source checkpoint lock entry malformed: {index}")
        if entry["locus"] != str(locus["LOC"]) or entry["chromosome"] != str(locus["CHR"]):
            raise ContinuationError(f"source checkpoint coordinate drifted: {index}")
        bundle = _safe_relative(str(entry["bundle_path"]), "source bundle")
        attempt = _safe_relative(str(entry["attempt_path"]), "source attempt")
        for field in ("ready", "receipt", "result", "log", "semantic_validation"):
            _validate_record_shape(entry[field], f"source checkpoint {index} {field}")
        if validate_bundles:
            _require_canonical_source_entry(
                entry, index, supervisor, source_fingerprint,
                revalidate_semantics=revalidate_semantics,
            )
    return payload


def continuation_contract_payload(
    source_fingerprint: str, *, validate_source_bundles: bool = False,
) -> dict[str, object]:
    source_fingerprint = _require_source_fingerprint(source_fingerprint)
    source_lock = validate_source_lock(
        source_fingerprint, validate_bundles=validate_source_bundles,
        require_live_fingerprint=validate_source_bundles,
    )
    source_lock_record = _identity_record(source_lock_path(source_fingerprint))
    predecessor_lock = validate_predecessor_supersession_lock()
    v2_files = [
        Path(__file__).resolve(), V2_RUNNER, V2_SUPERVISOR,
        V2_CHECKPOINT_VALIDATOR, V2_RESULT_VALIDATOR, SUPERSESSION_VALIDATOR,
    ]
    code: dict[str, dict[str, object]] = {}
    for path in v2_files:
        size, digest = file_identity(path)
        code[_relative(path)] = {"bytes": size, "sha256": digest}
    legacy = legacy_contract()
    return {
        "schema_version": CONTINUATION_CONTRACT_SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source_fingerprint,
        "source_family_lock": source_lock_record,
        "source_checkpoint_family_sha256": source_lock["checkpoint_family_sha256"],
        "active_input_baseline": source_lock["active_input_baseline"],
        "predecessor_supersession": {
            "relationship": "SUCCESSOR_CONTRACT_BINDS_FINAL_PREDECESSOR_LOCK",
            "superseded_continuation_fingerprint": PREDECESSOR_CONTINUATION_FINGERPRINT,
            "final_lock": predecessor_lock,
        },
        "scientific_inputs": {
            "policy_sha256": sha256(ROOT / "config/track_b_local_analysis_policy.json"),
            "input_lock_sha256": sha256(ROOT / "results/track_b/local_analysis_input.lock.json"),
            "reference_sha256": sha256(ROOT / "ref/lava/ukb_v1.1/reference.provenance.json"),
            "chromosome_input_lock_sha256": sha256(ROOT / "results/track_b/lava_chromosome_inputs.provenance.json"),
        },
        "runtime_contract": legacy.runtime_contract_payload(),
        "continuation_code": code,
        "allowed_phases": ["aggregate-discovery", "conditional", "finalize"],
        "scientific_equivalence": {
            "discovery_recomputed": False,
            "aggregate_rule": "ALL_2495_SOURCE_CHECKPOINTS_THEN_FULL_FAMILY_BH",
            "failed_family_member_rule": "FAILED_OR_UNTESTED_ELIGIBLE_MEMBERS_CONTRIBUTE_P_EQUALS_1",
            "patch_scope": (
                "RENAME_POLICY_INTEGER_ACCESSOR_USE_BASE_INTEGER_CONSTRUCTORS_"
                "SKIP_UNUSED_CHROMOSOME_IO_FOR_SEALED_NONCANDIDATE_LOCI_"
                "AND_REMOVE_REDUNDANT_SUPERVISOR_VALIDATION_WITHOUT_CHANGING_"
                "SCIENTIFIC_OR_ACTIVE_INPUT_GATES"
            ),
        },
    }


def continuation_fingerprint(
    source_fingerprint: str, *, validate_source_bundles: bool = False,
) -> str:
    return _digest(continuation_contract_payload(
        source_fingerprint, validate_source_bundles=validate_source_bundles,
    ))


def _lineage_payload(
    source_fingerprint: str, continuation_fingerprint_value: str,
) -> dict[str, object]:
    return {
        "schema_version": LINEAGE_SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source_fingerprint,
        "continuation_execution_fingerprint": continuation_fingerprint_value,
        "source_family_lock": _identity_record(source_lock_path(source_fingerprint)),
        "predecessor_supersession_lock": validate_predecessor_supersession_lock(),
        "contract_payload_sha256": continuation_fingerprint_value,
        "allowed_phases": ["aggregate-discovery", "conditional", "finalize"],
    }


def freeze_lineage(source_fingerprint: str, continuation_fingerprint_value: str) -> dict[str, object]:
    expected_fingerprint = continuation_fingerprint(source_fingerprint, validate_source_bundles=True)
    if continuation_fingerprint_value != expected_fingerprint:
        raise ContinuationError("supplied continuation fingerprint differs from live V2 contract")
    payload = _lineage_payload(source_fingerprint, expected_fingerprint)
    encoded = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    path = lineage_path(continuation_fingerprint_value)
    if path.exists():
        observed = validate_lineage_record(source_fingerprint, continuation_fingerprint_value)
        if observed != payload or path.read_bytes() != encoded:
            raise ContinuationError("existing continuation lineage differs; overwrite is forbidden")
        return observed
    publish_no_replace(path, encoded)
    return validate_lineage_record(source_fingerprint, continuation_fingerprint_value)


def validate_lineage_record(
    source_fingerprint: str, continuation_fingerprint_value: str,
) -> dict[str, object]:
    """Validate lineage bytes after the caller has checked the live contract."""

    source_fingerprint = _require_source_fingerprint(source_fingerprint)
    if (
        not FINGERPRINT_RE.fullmatch(continuation_fingerprint_value)
        or continuation_fingerprint_value == source_fingerprint
    ):
        raise ContinuationError("invalid continuation lineage fingerprint")
    path = lineage_path(continuation_fingerprint_value)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContinuationError(f"continuation lineage is unreadable: {error}") from error
    expected = _lineage_payload(source_fingerprint, continuation_fingerprint_value)
    if payload != expected:
        raise ContinuationError("continuation lineage differs from live contract identity")
    return payload


def validate_lineage(
    source_fingerprint: str, continuation_fingerprint_value: str, *, deep: bool = False,
) -> dict[str, object]:
    expected_fingerprint = continuation_fingerprint(
        source_fingerprint, validate_source_bundles=deep,
    )
    if continuation_fingerprint_value != expected_fingerprint:
        raise ContinuationError("continuation fingerprint differs from live V2 contract")
    return validate_lineage_record(source_fingerprint, continuation_fingerprint_value)


def preflight(source_fingerprint: str) -> tuple[str, str]:
    validate_source_lock(
        source_fingerprint, validate_bundles=True,
        revalidate_semantics=False, require_live_fingerprint=True,
    )
    continuation = continuation_fingerprint(source_fingerprint)
    freeze_lineage(source_fingerprint, continuation)
    return source_fingerprint, continuation


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--build-source-lock", action="store_true")
    mode.add_argument("--validate-source-lock", action="store_true")
    mode.add_argument("--continuation-fingerprint", action="store_true")
    mode.add_argument("--quick-continuation-fingerprint", action="store_true")
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--validate-lineage", action="store_true")
    parser.add_argument("--source-fingerprint", default=PINNED_SOURCE_FINGERPRINT)
    parser.add_argument("--expected-continuation-fingerprint")
    parser.add_argument("--deep", action="store_true")
    args = parser.parse_args()
    try:
        if args.build_source_lock:
            payload = build_source_lock(args.source_fingerprint, deep=args.deep)
            print(
                "TRACK_B_LAVA_SOURCE_FAMILY_LOCKED "
                f"source={payload['source_discovery_fingerprint']} checkpoints={payload['checkpoint_count']}"
            )
        elif args.validate_source_lock:
            payload = validate_source_lock(
                args.source_fingerprint, validate_bundles=True,
                revalidate_semantics=args.deep, require_live_fingerprint=True,
            )
            print(
                "TRACK_B_LAVA_SOURCE_FAMILY_VALIDATED "
                f"source={payload['source_discovery_fingerprint']} checkpoints={payload['checkpoint_count']}"
            )
        elif args.continuation_fingerprint or args.quick_continuation_fingerprint:
            print(continuation_fingerprint(
                args.source_fingerprint,
                validate_source_bundles=args.continuation_fingerprint and args.deep,
            ))
        elif args.preflight:
            source, continuation = preflight(args.source_fingerprint)
            print(
                "TRACK_B_LAVA_CONTINUATION_PREFLIGHT_PASS "
                f"source={source} continuation={continuation} discovery_recomputed=FALSE"
            )
        else:
            if not args.expected_continuation_fingerprint:
                parser.error("--validate-lineage requires --expected-continuation-fingerprint")
            validate_lineage(
                args.source_fingerprint, args.expected_continuation_fingerprint,
                deep=args.deep,
            )
            print("TRACK_B_LAVA_CONTINUATION_LINEAGE_VALIDATED")
    except ContinuationError as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
