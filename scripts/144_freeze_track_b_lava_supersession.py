#!/usr/bin/env python3
"""Freeze or verify the final evidence for a superseded LAVA continuation.

The supersession JSON is a pre-freeze narrative record.  It was edited while
the performance optimization was being reviewed, so its own ``immutable`` and
``NO_REPLACE`` fields are historical claims, not proof of immutability.  This
tool leaves those bytes untouched and publishes a separate no-replace lock
whose disclosure makes that limitation explicit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import stat
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SUPERSESSION_ROOT = ROOT / "results/track_b/lava_continuations/supersessions"
FINAL_LOCK_ROOT = SUPERSESSION_ROOT / "final_locks"
RUN_ROOT = ROOT / "results/track_b/checkpoints/lava_continuations/runs"
LOCK_SCHEMA = "track-b-lava-continuation-supersession-final-lock.1"
DIRECT_LOCK_SCHEMA = "track-b-lava-continuation-supersession-final-lock.2"
READY_SCHEMA = "track-b-lava-continuation-ready.1"
FINGERPRINT_RE = re.compile(r"^[0-9a-f]{64}$")
CANDIDATE_HEADER = (
    "locus_index\tlocus\tchromosome\tstart\tstop\tn_snps\t"
    "eligible_models\teligible_pairs\n"
).encode("utf-8")
RECORD_FIELDS = {
    "aggregate", "benchmark", "conditional_partial_family",
    "continuation_execution_fingerprint", "immutable", "lineage",
    "performance_evidence", "recorded_date", "registry_rule",
    "resume_old_namespace", "schema_version", "scientific_failure",
    "scientific_results_relabelled", "source_discovery_fingerprint", "status",
    "stop_reason", "successor_continuation_fingerprint",
}
V1_VALIDATOR_IDENTITY = {
    "path": "scripts/144_freeze_track_b_lava_supersession.py",
    "bytes": 23828,
    "sha256": "5ab203388415faaee27120e46542ca5632b22007e924dfe924f277a5015221a9",
}


class SupersessionError(RuntimeError):
    """A supersession record or frozen namespace failed closed."""


def _canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _digest(value: object) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _file_identity(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    try:
        if path.is_symlink():
            raise SupersessionError(f"frozen artifact may not be a symlink: {path}")
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size <= 0:
                raise SupersessionError(f"frozen artifact is missing or empty: {path}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        raise SupersessionError(f"could not identify frozen artifact {path}: {error}") from error
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before) != identity(after) or identity(after) != identity(current):
        raise SupersessionError(f"frozen artifact changed while hashing: {path}")
    return int(after.st_size), digest.hexdigest()


def _relative(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError as error:
        raise SupersessionError(f"path escapes repository: {path}") from error


def _safe_record_path(value: object, label: str) -> Path:
    if not isinstance(value, str):
        raise SupersessionError(f"{label} path is not a string")
    relative = Path(value)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise SupersessionError(f"unsafe {label} path: {value}")
    path = ROOT / relative
    try:
        path.resolve(strict=True).relative_to(ROOT.resolve(strict=True))
    except (OSError, RuntimeError, ValueError) as error:
        raise SupersessionError(f"{label} path is missing or escapes repository: {value}") from error
    return path


def _identity(path: Path) -> dict[str, object]:
    size, digest = _file_identity(path)
    return {"path": _relative(path), "bytes": size, "sha256": digest}


def _verify_record_identity(record: object, label: str) -> Path:
    if not isinstance(record, dict) or set(record) != {"path", "bytes", "sha256"}:
        raise SupersessionError(f"malformed {label} identity")
    if isinstance(record["bytes"], bool) or not isinstance(record["bytes"], int):
        raise SupersessionError(f"malformed {label} byte count")
    if record["bytes"] <= 0 or not isinstance(record["sha256"], str) or not FINGERPRINT_RE.fullmatch(record["sha256"]):
        raise SupersessionError(f"malformed {label} content identity")
    path = _safe_record_path(record["path"], label)
    if _identity(path) != record:
        raise SupersessionError(f"{label} differs from the supersession record")
    return path


def _load_json(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SupersessionError(f"{label} is unreadable: {error}") from error
    if not isinstance(value, dict):
        raise SupersessionError(f"{label} is not a JSON object")
    return value


def _direct_attempt(path: Path, attempts: Path, label: str) -> Path:
    if path.is_symlink():
        raise SupersessionError(f"{label} attempt may not be a symlink")
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise SupersessionError(f"{label} attempt is unreadable") from error
    if resolved.parent != attempts.resolve(strict=True) or not resolved.is_dir():
        raise SupersessionError(f"{label} attempt is outside the frozen namespace")
    return resolved


def _bundle_target(bundle: Path, attempts: Path, label: str) -> Path:
    if not bundle.is_symlink():
        raise SupersessionError(f"{label} is not an atomic READY link")
    try:
        target = bundle.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise SupersessionError(f"{label} READY link is broken") from error
    return _direct_attempt(target, attempts, label)


def _validate_bundle(
    target: Path, phase: str, index: int | None, source: str, continuation: str,
) -> dict[str, object]:
    required = {"READY", "receipt.json", "worker.log", "result.rds", "semantic_validation.txt"}
    live = {item.name for item in target.iterdir() if item.is_file() and not item.is_symlink()}
    if not required.issubset(live):
        raise SupersessionError(f"incomplete frozen bundle: {target}")
    receipt_path = target / "receipt.json"
    receipt = _load_json(receipt_path, f"{phase} receipt")
    ready = _load_json(target / "READY", f"{phase} READY marker")
    receipt_identity = _identity(receipt_path)
    expected_ready = {
        "schema_version": READY_SCHEMA,
        "receipt_sha256": receipt_identity["sha256"],
        "state": "READY",
    }
    if ready != expected_ready:
        raise SupersessionError(f"{phase} READY marker differs from its receipt")
    expected_index = index
    if (
        receipt.get("phase") != phase
        or receipt.get("locus_index") != expected_index
        or receipt.get("source_discovery_fingerprint") != source
        or receipt.get("continuation_execution_fingerprint") != continuation
        or receipt.get("exit_status") != 0
    ):
        raise SupersessionError(f"{phase} receipt identity drifted")
    marker = receipt.get("marker")
    if not isinstance(marker, dict) or marker.get("phase") != phase:
        raise SupersessionError(f"{phase} receipt marker drifted")
    if phase == "conditional" and marker.get("index") != str(index):
        raise SupersessionError(f"conditional marker index drifted: {index}")
    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise SupersessionError(f"{phase} receipt has no artifact family")
    artifact_names: set[str] = set()
    result_record: dict[str, object] | None = None
    for artifact in artifacts:
        if not isinstance(artifact, dict) or set(artifact) != {"path", "bytes", "sha256"}:
            raise SupersessionError(f"malformed {phase} artifact receipt")
        name = artifact["path"]
        if not isinstance(name, str) or Path(name).name != name or name.startswith(".") or name in artifact_names:
            raise SupersessionError(f"unsafe or duplicated {phase} artifact name")
        path = target / name
        observed = _identity(path)
        if observed["bytes"] != artifact["bytes"] or observed["sha256"] != artifact["sha256"]:
            raise SupersessionError(f"{phase} artifact differs from receipt: {name}")
        artifact_names.add(name)
        if name == "result.rds":
            result_record = observed
    if live != artifact_names | {"READY", "receipt.json", "worker.log"}:
        raise SupersessionError(f"{phase} bundle has unreceipted files")
    if result_record is None:
        raise SupersessionError(f"{phase} bundle lacks result.rds")
    if (
        receipt.get("log_bytes") != (target / "worker.log").stat().st_size
        or receipt.get("log_sha256") != _identity(target / "worker.log")["sha256"]
    ):
        raise SupersessionError(f"{phase} worker log differs from receipt")
    active_inputs = receipt.get("active_inputs")
    if not isinstance(active_inputs, list) or not active_inputs:
        raise SupersessionError(f"{phase} receipt lacks its historical active-input family")
    active_paths: set[str] = set()
    for item in active_inputs:
        expected_fields = {"path", "bytes", "sha256", "role", "pre_stat", "post_stat"}
        if not isinstance(item, dict) or set(item) != expected_fields:
            raise SupersessionError(f"malformed historical active-input record in {phase}")
        if (
            not isinstance(item["path"], str) or item["path"] in active_paths
            or isinstance(item["bytes"], bool) or not isinstance(item["bytes"], int)
            or item["bytes"] <= 0 or not isinstance(item["sha256"], str)
            or not FINGERPRINT_RE.fullmatch(item["sha256"])
            or item["pre_stat"] != item["post_stat"]
        ):
            raise SupersessionError(f"invalid historical active-input identity in {phase}")
        active_paths.add(item["path"])
    return {
        "ready": _identity(target / "READY"),
        "receipt": receipt_identity,
        "result": result_record,
    }


def _validate_supersession_record(fingerprint: str) -> tuple[Path, dict[str, Any]]:
    if not FINGERPRINT_RE.fullmatch(fingerprint):
        raise SupersessionError("invalid superseded continuation fingerprint")
    record_path = SUPERSESSION_ROOT / f"{fingerprint}.json"
    record = _load_json(record_path, "supersession record")
    if set(record) != RECORD_FIELDS:
        raise SupersessionError("supersession record schema drifted")
    source = record.get("source_discovery_fingerprint")
    successor = record.get("successor_continuation_fingerprint")
    partial = record.get("conditional_partial_family")
    if (
        record.get("schema_version") != "track-b-lava-continuation-supersession.1"
        or record.get("continuation_execution_fingerprint") != fingerprint
        or not isinstance(source, str) or not FINGERPRINT_RE.fullmatch(source)
        or not isinstance(successor, str) or not FINGERPRINT_RE.fullmatch(successor)
        or len({source, fingerprint, successor}) != 3
        or record.get("status") != "SUPERSEDED_PARTIAL_PERFORMANCE_OPTIMIZATION"
        or record.get("resume_old_namespace") is not False
        or record.get("scientific_failure") is not False
        or record.get("scientific_results_relabelled") is not False
        or record.get("immutable") is not True
        or record.get("registry_rule") != "NO_REPLACE_ONE_RECORD_PER_SUPERSEDED_CONTINUATION_FINGERPRINT"
        or not isinstance(partial, dict)
        or set(partial) != {"complete", "ready_count", "ready_indices", "receipt_result_family_sha256"}
        or partial.get("complete") is not False
        or partial.get("ready_count") != 10
        or partial.get("ready_indices") != list(range(1, 11))
        or not isinstance(partial.get("receipt_result_family_sha256"), str)
        or not FINGERPRINT_RE.fullmatch(partial["receipt_result_family_sha256"])
    ):
        raise SupersessionError("supersession identity or partial-family declaration drifted")
    return record_path, record


def build_lock_payload(fingerprint: str) -> dict[str, object]:
    record_path, record = _validate_supersession_record(fingerprint)
    source = str(record["source_discovery_fingerprint"])
    run = RUN_ROOT / fingerprint
    if run.is_symlink() or not run.is_dir():
        raise SupersessionError("superseded continuation namespace is missing or symlinked")
    attempts = run / ".attempts"
    if attempts.is_symlink() or not attempts.is_dir():
        raise SupersessionError("superseded attempt namespace is missing or symlinked")

    lineage_path = _verify_record_identity(record["lineage"], "lineage")
    benchmark_path = _verify_record_identity(record["benchmark"], "benchmark")
    if lineage_path != run / "lineage.json":
        raise SupersessionError("lineage record does not name the superseded run")
    lineage = _load_json(lineage_path, "lineage")
    if (
        lineage.get("source_discovery_fingerprint") != source
        or lineage.get("continuation_execution_fingerprint") != fingerprint
        or lineage.get("contract_payload_sha256") != fingerprint
    ):
        raise SupersessionError("superseded lineage identity drifted")

    aggregate = record.get("aggregate")
    if not isinstance(aggregate, dict) or set(aggregate) != {
        "candidate_attestation", "ready", "receipt", "result"
    } or aggregate.get("ready") is not True:
        raise SupersessionError("supersession aggregate declaration drifted")
    aggregate_receipt = _verify_record_identity(aggregate["receipt"], "aggregate receipt")
    aggregate_result = _verify_record_identity(aggregate["result"], "aggregate result")
    candidate_path = _verify_record_identity(
        aggregate["candidate_attestation"], "aggregate candidate attestation"
    )
    aggregate_target = _direct_attempt(aggregate_receipt.parent, attempts, "aggregate")
    if aggregate_result.parent != aggregate_target or candidate_path.parent != aggregate_target:
        raise SupersessionError("aggregate evidence does not share one frozen attempt")
    aggregate_bundle = run / "aggregate-discovery/unit_all"
    if _bundle_target(aggregate_bundle, attempts, "aggregate") != aggregate_target:
        raise SupersessionError("aggregate READY link differs from the recorded attempt")
    aggregate_verified = _validate_bundle(
        aggregate_target, "aggregate-discovery", None, source, fingerprint,
    )
    if aggregate_verified["receipt"] != aggregate["receipt"] or aggregate_verified["result"] != aggregate["result"]:
        raise SupersessionError("aggregate receipt/result differs from the narrative record")
    if candidate_path.read_bytes() != CANDIDATE_HEADER:
        raise SupersessionError("superseded aggregate candidate family is not exactly header-only")

    conditional_root = run / "conditional"
    if conditional_root.is_symlink() or not conditional_root.is_dir():
        raise SupersessionError("superseded conditional namespace is missing or symlinked")
    canonical_entries = sorted(conditional_root.iterdir(), key=lambda item: item.name)
    expected_names = [f"locus_{index:04d}" for index in range(1, 11)]
    if [item.name for item in canonical_entries] != expected_names:
        raise SupersessionError("conditional READY namespace is not exactly loci 1..10")

    members: list[dict[str, object]] = []
    ready_targets: set[Path] = set()
    for index, bundle in enumerate(canonical_entries, start=1):
        target = _bundle_target(bundle, attempts, f"conditional locus {index}")
        if target in ready_targets:
            raise SupersessionError("two conditional READY links share one attempt")
        ready_targets.add(target)
        verified = _validate_bundle(target, "conditional", index, source, fingerprint)
        members.append({
            "locus_index": index,
            "receipt": verified["receipt"],
            "result": verified["result"],
        })

    attempt_entries = sorted(attempts.iterdir(), key=lambda item: item.name)
    attempt_names = [item.name for item in attempt_entries]
    conditional_attempts = {
        _direct_attempt(item, attempts, "conditional inventory")
        for item in attempt_entries if item.name.startswith("conditional_locus_")
    }
    if conditional_attempts != ready_targets:
        raise SupersessionError("conditional attempt inventory is not exactly the ten READY attempts")
    ready_files = {
        item / "READY" for item in attempt_entries
        if item.is_dir() and not item.is_symlink() and (item / "READY").is_file()
    }
    expected_ready_files = {aggregate_target / "READY"} | {item / "READY" for item in ready_targets}
    if ready_files != expected_ready_files:
        raise SupersessionError("superseded run contains an unregistered READY attempt")
    for forbidden_phase in ("finalize", "terminal-qc"):
        phase_path = run / forbidden_phase
        if phase_path.exists() or phase_path.is_symlink():
            raise SupersessionError(f"superseded partial run unexpectedly has {forbidden_phase}")

    recipe = {
        "recipe_id": "SHA256_CANONICAL_SORTED_JSON_RECEIPT_RESULT_IDENTITIES_V1",
        "input": "the members array in ascending locus_index order",
        "member_schema": {
            "locus_index": "integer",
            "receipt": "repository-relative path, byte count, SHA-256",
            "result": "repository-relative path, byte count, SHA-256",
        },
        "canonicalization": (
            "json.dumps(members,sort_keys=True,separators=(',',':')).encode('utf-8'); "
            "no trailing newline"
        ),
        "hash": "SHA-256 lowercase hexadecimal",
    }
    frozen_family_digest = _digest(members)
    legacy_digest = record["conditional_partial_family"]["receipt_result_family_sha256"]
    return {
        "schema_version": LOCK_SCHEMA,
        "state": "FINAL_LOCKED",
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source,
        "superseded_continuation_fingerprint": fingerprint,
        "successor_continuation_fingerprint": record["successor_continuation_fingerprint"],
        "immutability_scope": {
            "effective_at": "FIRST_SUCCESSFUL_NO_REPLACE_PUBLICATION_OF_THIS_LOCK",
            "immutable": True,
            "covers": (
                "THE_EXACT_SUPERSESSION_JSON_BYTES_AND_FROZEN_NAMESPACE_IDENTITIES_"
                "LISTED_IN_THIS_LOCK_FROM_PUBLICATION_FORWARD"
            ),
            "does_not_retroactively_cover": (
                "ANY_PRE_FREEZE_DRAFT_VERSION_OR_EDIT_OF_THE_SUPERSESSION_JSON"
            ),
            "old_namespace_may_be_resumed": False,
        },
        "historical_disclosure": {
            "supersession_json_was_edited_during_drafting": True,
            "pre_freeze_draft_versions_attested": False,
            "legacy_json_immutable_claim_is_not_the_final_lock": True,
            "legacy_json_registry_rule_is_not_proof_of_prior_no_replace_publication": True,
        },
        "supersession_record": _identity(record_path),
        "legacy_pre_freeze_digest": {
            "field": "conditional_partial_family.receipt_result_family_sha256",
            "value_preserved_verbatim": legacy_digest,
            "recipe_available": False,
            "verification_status": "UNVERIFIABLE_PRE_FREEZE_DRAFT_VALUE_NOT_USED_AS_EVIDENCE",
        },
        "verified_snapshot": {
            "lineage": _identity(lineage_path),
            "benchmark": _identity(benchmark_path),
            "aggregate": {
                "ready": aggregate_verified["ready"],
                "receipt": aggregate_verified["receipt"],
                "result": aggregate_verified["result"],
                "candidate_attestation": _identity(candidate_path),
                "candidate_family": "HEADER_ONLY_ZERO_CANDIDATES",
            },
            "conditional_partial_family": {
                "complete": False,
                "ready_count": 10,
                "ready_indices": list(range(1, 11)),
                "no_other_ready_attempts": True,
                "digest_recipe": recipe,
                "members": members,
                "receipt_result_family_sha256": frozen_family_digest,
            },
            "attempt_directory_names": attempt_names,
        },
        "validator": _identity(Path(__file__).resolve()),
    }


def _family_digest_recipe() -> dict[str, object]:
    return {
        "recipe_id": "SHA256_CANONICAL_SORTED_JSON_RECEIPT_RESULT_IDENTITIES_V1",
        "input": "the members array in ascending locus_index order",
        "member_schema": {
            "locus_index": "integer",
            "receipt": "repository-relative path, byte count, SHA-256",
            "result": "repository-relative path, byte count, SHA-256",
        },
        "canonicalization": (
            "json.dumps(members,sort_keys=True,separators=(',',':')).encode('utf-8'); "
            "no trailing newline"
        ),
        "hash": "SHA-256 lowercase hexadecimal",
    }


def _direct_benchmark_snapshot(
    source: str, fingerprint: str, expected_ready_through: int,
) -> dict[str, object]:
    directory = ROOT / "results/track_b/lava_continuations/benchmarks" / fingerprint
    expected_names = {
        "RAM_BENCHMARK.tsv", "RAM_BENCHMARK.namespace.json", "RAM_BENCHMARK.provenance.json",
    }
    if directory.is_symlink() or not directory.is_dir():
        raise SupersessionError("continuation benchmark namespace is missing or symlinked")
    if {item.name for item in directory.iterdir()} != expected_names:
        raise SupersessionError("continuation benchmark namespace inventory drifted")
    table = directory / "RAM_BENCHMARK.tsv"
    namespace_path = directory / "RAM_BENCHMARK.namespace.json"
    provenance_path = directory / "RAM_BENCHMARK.provenance.json"
    namespace = _load_json(namespace_path, "benchmark namespace")
    provenance = _load_json(provenance_path, "benchmark provenance")
    if (
        namespace.get("source_discovery_fingerprint") != source
        or namespace.get("continuation_execution_fingerprint") != fingerprint
        or provenance.get("source_discovery_fingerprint") != source
        or provenance.get("continuation_execution_fingerprint") != fingerprint
        or provenance.get("namespace_sha256") != _identity(namespace_path)["sha256"]
        or provenance.get("table_sha256") != _identity(table)["sha256"]
        or provenance.get("row_count") != expected_ready_through + 1
    ):
        raise SupersessionError("continuation benchmark identity drifted")
    try:
        with table.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    except (OSError, UnicodeError, csv.Error) as error:
        raise SupersessionError(f"continuation benchmark table is unreadable: {error}") from error
    if len(rows) != expected_ready_through + 1:
        raise SupersessionError("continuation benchmark row count differs from the READY family")
    aggregate = [row for row in rows if row.get("phase") == "aggregate-discovery"]
    conditional = [row for row in rows if row.get("phase") == "conditional"]
    try:
        conditional_indices = [int(row["locus_index"]) for row in conditional]
    except (KeyError, TypeError, ValueError) as error:
        raise SupersessionError("continuation benchmark locus indices are malformed") from error
    if (
        len(aggregate) != 1 or aggregate[0].get("locus_index") != "NA"
        or conditional_indices != list(range(1, expected_ready_through + 1))
        or any(row.get("source_discovery_fingerprint") != source for row in rows)
        or any(row.get("continuation_execution_fingerprint") != fingerprint for row in rows)
    ):
        raise SupersessionError("continuation benchmark topology differs from the READY family")
    return {
        "table": _identity(table),
        "namespace": _identity(namespace_path),
        "provenance": _identity(provenance_path),
        "row_count": len(rows),
    }


def build_direct_lock_payload(
    fingerprint: str, source: str, expected_ready_through: int, *,
    stop_reason: str, performance_rationale: str, stop_boundary: str,
) -> dict[str, object]:
    """Build a prospective lock directly from one exact partial run namespace."""

    if (
        not FINGERPRINT_RE.fullmatch(fingerprint)
        or not FINGERPRINT_RE.fullmatch(source)
        or fingerprint == source
    ):
        raise SupersessionError("invalid direct supersession fingerprint identity")
    if (
        isinstance(expected_ready_through, bool)
        or not isinstance(expected_ready_through, int)
        or not 1 <= expected_ready_through < 2495
    ):
        raise SupersessionError("direct supersession requires a partial contiguous READY family")
    for label, value in (
        ("stop reason", stop_reason),
        ("performance rationale", performance_rationale),
        ("stop boundary", stop_boundary),
    ):
        if not isinstance(value, str) or not value or any(character in value for character in "\r\n\t"):
            raise SupersessionError(f"invalid {label}")

    run = RUN_ROOT / fingerprint
    if run.is_symlink() or not run.is_dir():
        raise SupersessionError("superseded continuation namespace is missing or symlinked")
    if {item.name for item in run.iterdir()} != {
        ".attempts", "aggregate-discovery", "conditional", "lineage.json",
    }:
        raise SupersessionError("superseded run namespace inventory is not exactly partial")
    attempts = run / ".attempts"
    if attempts.is_symlink() or not attempts.is_dir():
        raise SupersessionError("superseded attempt namespace is missing or symlinked")
    if (run / "failed_attempts").exists() or (run / "failed_attempts").is_symlink():
        raise SupersessionError("superseded run contains failed-attempt artifacts")
    for forbidden_phase in ("finalize", "terminal-qc"):
        phase_path = run / forbidden_phase
        if phase_path.exists() or phase_path.is_symlink():
            raise SupersessionError(f"superseded partial run unexpectedly has {forbidden_phase}")

    lineage_path = run / "lineage.json"
    lineage = _load_json(lineage_path, "lineage")
    if (
        lineage.get("source_discovery_fingerprint") != source
        or lineage.get("continuation_execution_fingerprint") != fingerprint
        or lineage.get("contract_payload_sha256") != fingerprint
    ):
        raise SupersessionError("superseded lineage identity drifted")
    source_lock_record = lineage.get("source_family_lock")
    source_lock_path = _verify_record_identity(source_lock_record, "source-family lock")

    aggregate_bundle = run / "aggregate-discovery/unit_all"
    aggregate_target = _bundle_target(aggregate_bundle, attempts, "aggregate")
    aggregate_verified = _validate_bundle(
        aggregate_target, "aggregate-discovery", None, source, fingerprint,
    )
    aggregate_receipt = _load_json(aggregate_target / "receipt.json", "aggregate receipt")
    candidate_receipt = next(
        (
            item for item in aggregate_receipt.get("artifacts", [])
            if isinstance(item, dict) and item.get("path") == "conditional_candidates.tsv"
        ),
        None,
    )
    if candidate_receipt is None:
        raise SupersessionError("aggregate receipt lacks its candidate attestation")
    candidate_path = aggregate_target / "conditional_candidates.tsv"
    candidate_identity = _identity(candidate_path)
    if (
        candidate_identity["bytes"] != candidate_receipt.get("bytes")
        or candidate_identity["sha256"] != candidate_receipt.get("sha256")
        or candidate_path.read_bytes() != CANDIDATE_HEADER
    ):
        raise SupersessionError("aggregate candidate seal is not exactly header-only")

    conditional_root = run / "conditional"
    if conditional_root.is_symlink() or not conditional_root.is_dir():
        raise SupersessionError("superseded conditional namespace is missing or symlinked")
    canonical_entries = sorted(conditional_root.iterdir(), key=lambda item: item.name)
    expected_names = [
        f"locus_{index:04d}" for index in range(1, expected_ready_through + 1)
    ]
    if [item.name for item in canonical_entries] != expected_names:
        raise SupersessionError(
            "conditional READY namespace is not the declared contiguous partial family"
        )

    members: list[dict[str, object]] = []
    ready_targets: set[Path] = set()
    for index, bundle in enumerate(canonical_entries, start=1):
        target = _bundle_target(bundle, attempts, f"conditional locus {index}")
        if target in ready_targets:
            raise SupersessionError("two conditional READY links share one attempt")
        ready_targets.add(target)
        verified = _validate_bundle(target, "conditional", index, source, fingerprint)
        members.append({
            "locus_index": index,
            "receipt": verified["receipt"],
            "result": verified["result"],
        })

    attempt_entries = sorted(attempts.iterdir(), key=lambda item: item.name)
    if any(item.is_symlink() or not item.is_dir() for item in attempt_entries):
        raise SupersessionError("attempt inventory contains a non-directory or symlink")
    aggregate_attempts = {
        _direct_attempt(item, attempts, "aggregate inventory")
        for item in attempt_entries if item.name.startswith("aggregate-discovery_")
    }
    conditional_attempts = {
        _direct_attempt(item, attempts, "conditional inventory")
        for item in attempt_entries if item.name.startswith("conditional_locus_")
    }
    if aggregate_attempts != {aggregate_target} or conditional_attempts != ready_targets:
        raise SupersessionError("attempt inventory differs from the exact READY family")
    if len(attempt_entries) != expected_ready_through + 1:
        raise SupersessionError("superseded namespace contains an extra attempt")
    ready_files = {
        item / "READY" for item in attempt_entries if (item / "READY").is_file()
    }
    expected_ready_files = {aggregate_target / "READY"} | {
        item / "READY" for item in ready_targets
    }
    if ready_files != expected_ready_files:
        raise SupersessionError("superseded run contains an unregistered READY attempt")

    benchmark = _direct_benchmark_snapshot(source, fingerprint, expected_ready_through)
    family_digest = _digest(members)
    return {
        "schema_version": DIRECT_LOCK_SCHEMA,
        "state": "FINAL_LOCKED",
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": source,
        "superseded_continuation_fingerprint": fingerprint,
        "status": "SUPERSEDED_PARTIAL_PERFORMANCE_OPTIMIZATION",
        "stop_reason": stop_reason,
        "performance_rationale": performance_rationale,
        "stop_boundary": stop_boundary,
        "scientific_failure": False,
        "scientific_results_relabelled": False,
        "successor_binding": {
            "mode": "SUCCESSOR_CONTRACT_BINDS_THIS_LOCK",
            "successor_fingerprint": None,
        },
        "immutability_scope": {
            "effective_at": "FIRST_SUCCESSFUL_NO_REPLACE_PUBLICATION_OF_THIS_LOCK",
            "immutable": True,
            "covers": "THIS_LOCK_AND_THE_EXACT_LISTED_RUN_STATE_FROM_PUBLICATION_FORWARD",
            "does_not_retroactively_cover": "ANY_STATE_BEFORE_THIS_DIRECT_FINAL_PUBLICATION",
            "old_namespace_may_be_resumed": False,
        },
        "verified_snapshot": {
            "source_family_lock": _identity(source_lock_path),
            "lineage": _identity(lineage_path),
            "benchmark": benchmark,
            "aggregate": {
                "ready": aggregate_verified["ready"],
                "receipt": aggregate_verified["receipt"],
                "result": aggregate_verified["result"],
                "candidate_attestation": candidate_identity,
                "candidate_family": "HEADER_ONLY_ZERO_CANDIDATES",
            },
            "conditional_partial_family": {
                "complete": False,
                "ready_count": expected_ready_through,
                "ready_indices": list(range(1, expected_ready_through + 1)),
                "no_other_ready_attempts": True,
                "digest_recipe": _family_digest_recipe(),
                "members": members,
                "receipt_result_family_sha256": family_digest,
            },
            "attempt_directory_names": [item.name for item in attempt_entries],
            "failed_attempt_count": 0,
            "finalize_present": False,
            "terminal_qc_present": False,
        },
        "validator": _identity(Path(__file__).resolve()),
    }


def lock_path(fingerprint: str) -> Path:
    if not FINGERPRINT_RE.fullmatch(fingerprint):
        raise SupersessionError("invalid superseded continuation fingerprint")
    return FINAL_LOCK_ROOT / f"{fingerprint}.lock.json"


def _lock_bytes(payload: dict[str, object]) -> bytes:
    return json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def publish_no_replace(path: Path, content: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise SupersessionError(f"refusing to replace existing final lock: {path}") from error
        _fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def freeze(fingerprint: str) -> Path:
    path = lock_path(fingerprint)
    payload = build_lock_payload(fingerprint)
    publish_no_replace(path, _lock_bytes(payload))
    verify(fingerprint)
    return path


def freeze_direct(
    fingerprint: str, source: str, expected_ready_through: int, *,
    stop_reason: str, performance_rationale: str, stop_boundary: str,
) -> Path:
    path = lock_path(fingerprint)
    payload = build_direct_lock_payload(
        fingerprint, source, expected_ready_through,
        stop_reason=stop_reason,
        performance_rationale=performance_rationale,
        stop_boundary=stop_boundary,
    )
    publish_no_replace(path, _lock_bytes(payload))
    verify(fingerprint)
    return path


def verify(fingerprint: str) -> dict[str, object]:
    path = lock_path(fingerprint)
    try:
        observed_bytes = path.read_bytes()
        observed = json.loads(observed_bytes)
    except (OSError, json.JSONDecodeError) as error:
        raise SupersessionError(f"final supersession lock is unreadable: {error}") from error
    if not isinstance(observed, dict):
        raise SupersessionError("final supersession lock is not a JSON object")
    if observed.get("schema_version") == LOCK_SCHEMA:
        expected = build_lock_payload(fingerprint)
        # Schema-.1 bound the exact validator bytes that first published it.
        # Preserve those already-frozen bytes while allowing this tool to add
        # prospective schema-.2 support; never rewrite the existing lock.
        expected["validator"] = V1_VALIDATOR_IDENTITY
    elif observed.get("schema_version") == DIRECT_LOCK_SCHEMA:
        snapshot = observed.get("verified_snapshot")
        if not isinstance(snapshot, dict):
            raise SupersessionError("direct final lock snapshot is malformed")
        family = snapshot.get("conditional_partial_family")
        if not isinstance(family, dict):
            raise SupersessionError("direct final lock conditional family is malformed")
        expected = build_direct_lock_payload(
            fingerprint,
            str(observed.get("source_discovery_fingerprint", "")),
            family.get("ready_count"),
            stop_reason=observed.get("stop_reason"),
            performance_rationale=observed.get("performance_rationale"),
            stop_boundary=observed.get("stop_boundary"),
        )
    else:
        raise SupersessionError("unsupported final supersession lock schema")
    if observed != expected or observed_bytes != _lock_bytes(expected):
        raise SupersessionError("final supersession lock differs from the frozen current state")
    return observed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--freeze", action="store_true")
    mode.add_argument("--freeze-direct", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument("--continuation-fingerprint", required=True)
    parser.add_argument("--source-fingerprint")
    parser.add_argument("--expected-ready-through", type=int)
    parser.add_argument("--stop-reason")
    parser.add_argument("--performance-rationale")
    parser.add_argument("--stop-boundary")
    args = parser.parse_args()
    try:
        if args.freeze:
            if any(value is not None for value in (
                args.source_fingerprint, args.expected_ready_through, args.stop_reason,
                args.performance_rationale, args.stop_boundary,
            )):
                parser.error("legacy --freeze does not accept direct-lock arguments")
            path = freeze(args.continuation_fingerprint)
            print(
                "TRACK_B_LAVA_SUPERSESSION_FINAL_LOCKED "
                f"continuation={args.continuation_fingerprint} lock={path}"
            )
        elif args.freeze_direct:
            required = (
                args.source_fingerprint, args.expected_ready_through, args.stop_reason,
                args.performance_rationale, args.stop_boundary,
            )
            if any(value is None for value in required):
                parser.error("--freeze-direct requires source, READY bound, reason, rationale, and boundary")
            path = freeze_direct(
                args.continuation_fingerprint, args.source_fingerprint,
                args.expected_ready_through, stop_reason=args.stop_reason,
                performance_rationale=args.performance_rationale,
                stop_boundary=args.stop_boundary,
            )
            print(
                "TRACK_B_LAVA_SUPERSESSION_DIRECT_FINAL_LOCKED "
                f"continuation={args.continuation_fingerprint} lock={path}"
            )
        else:
            if any(value is not None for value in (
                args.source_fingerprint, args.expected_ready_through, args.stop_reason,
                args.performance_rationale, args.stop_boundary,
            )):
                parser.error("--verify derives its exact scope from the no-replace lock")
            payload = verify(args.continuation_fingerprint)
            family = payload["verified_snapshot"]["conditional_partial_family"]
            print(
                "TRACK_B_LAVA_SUPERSESSION_FINAL_LOCK_VERIFIED "
                f"continuation={args.continuation_fingerprint} "
                f"ready={family['ready_count']} digest={family['receipt_result_family_sha256']}"
            )
    except SupersessionError as error:
        raise SystemExit(f"ERROR: {error}") from error


if __name__ == "__main__":
    main()
