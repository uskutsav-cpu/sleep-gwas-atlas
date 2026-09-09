#!/usr/bin/env python3
"""Verify and publish the final Track B RAM federation.

This is a reporting-only boundary.  It does not execute or reinterpret any
scientific model.  A candidate is assembled only after every contributing
execution family has reached, and still verifies at, its immutable terminal
boundary.  ``--preflight`` and ``--verify`` never write in the repository.
Only ``--execute`` may create the versioned package and promote the four
historical top-level RAM artifacts; the provenance file is the last promotion
and therefore the commit marker.
"""

from __future__ import annotations

import argparse
import base64
import collections
import contextlib
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
import shutil
import stat
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_REL = Path("scripts/153_finalize_track_b_ram_aware_execution.py")

V1_FINGERPRINT = "caeb6b5a1188560f27609cad77801715e0bedfb998660a531c5491afad1177c7"
CONTINUATION_FINGERPRINT = "c139e2368cae8f4b88c227a3f0a541b3ad988d941ce881a4b1c1aef9e28b8fc7"
SUPERSEDED_CONTINUATIONS = (
    "09cd41888ba405ae2b54b8584c5bc15c7d0213e81e7c1cd8f661136bca34a152",
    "ceb36dad6534e5f3256c32e53e2a5f0d8df3a56ad68cab93ebf96843a5fa1461",
)
EXPECTED_DISCOVERY_COUNT = 2495
EXPECTED_HISTORICAL_ROWS = 2496
EXPECTED_STALE_DECLARED_ROWS = 1376

HISTORY_SCRIPT_REL = Path("scripts/152_freeze_track_b_ram_v1_history.py")
CONTINUATION_CONTRACT_REL = Path("scripts/134_track_b_lava_continuation_contract.py")
CONTINUATION_RUNNER_REL = Path("scripts/136_run_track_b_lava_continuation.py")
CONTINUATION_RESULTS_REL = Path("scripts/138_validate_track_b_lava_results_v2.py")
SUPERSESSION_REL = Path("scripts/144_freeze_track_b_lava_supersession.py")
PLACO_REL = Path("scripts/143_run_track_b_placo_sequential_v2.py")
COLLATION_REL = Path("scripts/147_build_track_b_pleiotropy_results_v2.py")
FINEMAP_REL = Path("scripts/151_run_track_b_finemapping_sequential.py")
FINEMAP_ZERO_STATE = "NOT_APPLICABLE_ZERO_ELIGIBLE_LOCUS_FAMILY"
FINEMAP_VERIFIER_SHA256 = "a4f433830151f1fe80b707d7ce48f10bce432cbef2afb52c8e578244edd870c6"
EXPECTED_PLACO_PAIR_ORDER = ("A", "B", "CONTROL")

BENCHMARK_FIELDS = [
    "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
    "runtime_sec", "exit_status", "output_hash",
]
COLLISION_KEY = ["analysis", "pair", "locus", "chromosome"]
TOP_LEVEL = {
    "benchmark": Path("results/track_b/RAM_BENCHMARK.tsv"),
    "namespace": Path("results/track_b/RAM_BENCHMARK.namespace.json"),
    "report": Path("results/track_b/RAM_AWARE_EXECUTION_REPORT.md"),
    "provenance": Path("results/track_b/RAM_BENCHMARK.provenance.json"),
}
PROMOTION_ORDER = ("benchmark", "namespace", "report", "provenance")
PACKAGE_ROOT_REL = Path("results/track_b/ram_federation")
TRANSACTION_ROOT_REL = PACKAGE_ROOT_REL / ".transactions"
PUBLICATION_LOCK_REL = PACKAGE_ROOT_REL / ".publication.lock"
PACKAGE_FILES = {
    "benchmark": "RAM_BENCHMARK.tsv",
    "namespace": "RAM_BENCHMARK.namespace.json",
    "report": "RAM_AWARE_EXECUTION_REPORT.md",
    "provenance": "RAM_BENCHMARK.provenance.json",
}
PACKAGE_MANIFEST = "package.manifest.json"
SCHEMA = "sleep-atlas-track-b-final-ram-federation.1"
NAMESPACE_SCHEMA = "sleep-atlas-track-b-final-ram-namespace.1"
PACKAGE_SCHEMA = "sleep-atlas-track-b-final-ram-package.1"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
INTEGER_RE = re.compile(r"^(?:0|[1-9][0-9]*)$")
GIB = 1024**3
EIGHT_GIB = 8 * GIB
RESERVE_BYTES = GIB


class FederationError(RuntimeError):
    """An input family is incomplete, unsafe, stale, or internally inconsistent."""


@contextlib.contextmanager
def _suppress_bytecode_writes():
    """Suppress bytecode for the complete verifier call graph, including lazy imports."""

    previous = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        yield
    finally:
        sys.dont_write_bytecode = previous


def canonical_json_bytes(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode("utf-8")


def digest_json(value: object) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _load_module(name: str, path: Path) -> Any:
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise FederationError(f"could not import required verifier: {path}")
    module = importlib.util.module_from_spec(specification)
    with _suppress_bytecode_writes():
        specification.loader.exec_module(module)
    return module


def safe_path(root: Path, relative: Path | str, label: str) -> Path:
    relative = Path(relative)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise FederationError(f"{label} is not a safe repository-relative path: {relative}")
    root_real = root.resolve(strict=True)
    current = root
    for part in relative.parts:
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise FederationError(f"could not inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise FederationError(f"{label} contains a symbolic link: {current}")
    try:
        (root / relative).resolve(strict=False).relative_to(root_real)
    except (OSError, RuntimeError, ValueError) as error:
        raise FederationError(f"{label} escapes the repository: {relative}") from error
    return root / relative


def stable_bytes(path: Path, label: str, *, allow_empty: bool = False) -> bytes:
    try:
        first = os.lstat(path)
        if (
            not stat.S_ISREG(first.st_mode)
            or first.st_size < 0
            or (first.st_size == 0 and not allow_empty)
        ):
            raise FederationError(f"{label} is not a stable regular file: {path}")
        with path.open("rb") as handle:
            opened = os.fstat(handle.fileno())
            content = handle.read()
            after = os.fstat(handle.fileno())
        final = os.lstat(path)
    except OSError as error:
        raise FederationError(f"could not read {label}: {path}: {error}") from error
    identities = {
        (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns)
        for item in (first, opened, after, final)
    }
    if len(identities) != 1 or len(content) != final.st_size:
        raise FederationError(f"{label} changed while it was read: {path}")
    return content


def portable_identity(root: Path, path: Path, label: str, *, allow_empty: bool = False) -> dict[str, object]:
    try:
        relative = path.relative_to(root)
    except ValueError as error:
        raise FederationError(f"{label} is outside the repository: {path}") from error
    path = safe_path(root, relative, label)
    content = stable_bytes(path, label, allow_empty=allow_empty)
    return {
        "path": str(relative), "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def identity_for_bytes(path: Path | str, content: bytes) -> dict[str, object]:
    return {
        "path": str(path), "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def _verify_identity_record(
    root: Path, record: object, label: str, *, allow_empty: bool = False,
) -> tuple[Path, dict[str, object]]:
    if not isinstance(record, Mapping) or set(record) != {"path", "bytes", "sha256"}:
        raise FederationError(f"malformed {label} identity")
    path = safe_path(root, str(record["path"]), label)
    observed = portable_identity(root, path, label, allow_empty=allow_empty)
    if observed != dict(record):
        raise FederationError(f"{label} differs from its immutable identity")
    return path, observed


def read_json(path: Path, label: str) -> tuple[dict[str, Any], bytes]:
    content = stable_bytes(path, label)
    try:
        value = json.loads(content.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise FederationError(f"{label} is not UTF-8 JSON: {error}") from error
    if not isinstance(value, dict):
        raise FederationError(f"{label} must contain one JSON object")
    return value, content


def read_tsv_bytes(content: bytes, fields: Sequence[str], label: str) -> list[dict[str, str]]:
    try:
        reader = csv.DictReader(io.StringIO(content.decode("utf-8")), delimiter="\t")
        if list(reader.fieldnames or []) != list(fields):
            raise FederationError(f"{label} does not have the exact required schema")
        rows = list(reader)
    except (UnicodeError, csv.Error) as error:
        raise FederationError(f"{label} is not a readable UTF-8 TSV: {error}") from error
    for number, row in enumerate(rows, start=2):
        if None in row or set(row) != set(fields) or any(row[field] in {None, ""} for field in fields):
            raise FederationError(f"{label} has a malformed row {number}")
    return rows


def read_tsv(path: Path, fields: Sequence[str], label: str) -> tuple[list[dict[str, str]], bytes]:
    content = stable_bytes(path, label)
    return read_tsv_bytes(content, fields, label), content


def tsv_bytes(fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(fields), delimiter="\t", lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row[field] for field in fields})
    return buffer.getvalue().encode("utf-8")


def validate_federation_row(row: Mapping[str, object], label: str) -> dict[str, str]:
    if set(row) != set(BENCHMARK_FIELDS):
        raise FederationError(f"{label} does not have exactly the nine federation fields")
    normalized = {field: str(row[field]) for field in BENCHMARK_FIELDS}
    if any(not value or any(character in value for character in "\r\n\t") for value in normalized.values()):
        raise FederationError(f"{label} contains an empty or unsafe field")
    if normalized["n_snps"] != "NA" and not INTEGER_RE.fullmatch(normalized["n_snps"]):
        raise FederationError(f"{label} has an invalid SNP count")
    try:
        peak = float(normalized["peak_ram_gb"])
        runtime = float(normalized["runtime_sec"])
    except (ValueError, OverflowError) as error:
        raise FederationError(f"{label} has invalid RAM/runtime values") from error
    if not math.isfinite(peak) or peak < 0 or not math.isfinite(runtime) or runtime < 0:
        raise FederationError(f"{label} has non-finite or negative RAM/runtime values")
    if not SHA256_RE.fullmatch(normalized["output_hash"]):
        raise FederationError(f"{label} has an invalid output SHA-256")
    return normalized


def merge_component_rows(components: Sequence[Mapping[str, Any]]) -> list[dict[str, str]]:
    """Merge without filtering; reject superseded namespaces and exact-key collisions."""

    merged: list[dict[str, str]] = []
    seen: dict[tuple[str, str, str, str], str] = {}
    for component in components:
        name = str(component.get("name", "unnamed component"))
        fingerprint = component.get("execution_fingerprint")
        if fingerprint in SUPERSEDED_CONTINUATIONS:
            raise FederationError(f"superseded continuation cannot enter federation: {fingerprint}")
        rows = component.get("rows")
        if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
            raise FederationError(f"{name} has no row family")
        for index, row in enumerate(rows, start=1):
            if not isinstance(row, Mapping):
                raise FederationError(f"{name} row {index} is not a mapping")
            normalized = validate_federation_row(row, f"{name} row {index}")
            key = tuple(normalized[field] for field in COLLISION_KEY)
            if key in seen:
                raise FederationError(
                    "duplicate federation collision key "
                    f"{key!r} in {seen[key]} and {name}"
                )
            seen[key] = name
            merged.append(normalized)
    return merged


def validate_stale_v1_audit(
    audit: object, actual_rows: int, *, expected_discovery_count: int = EXPECTED_DISCOVERY_COUNT,
) -> None:
    if not isinstance(audit, Mapping):
        raise FederationError("V1 history receipt lacks its stale-provenance audit")
    if (
        audit.get("fully_consistent") is not False
        or audit.get("row_count_match") is not False
        or audit.get("benchmark_bytes_match") is not False
        or audit.get("benchmark_sha256_match") is not False
        or audit.get("observed_benchmark_rows") != actual_rows
        or audit.get("declared_attempt_rows") == actual_rows
        or audit.get("handling")
        != "PRESERVED_AS_PROVISIONAL_HISTORY_ONLY;FINAL_FEDERATION_MUST_REVALIDATE_RAW_RECEIPTS"
    ):
        raise FederationError("V1 stale provenance was not explicitly quarantined")
    if actual_rows != expected_discovery_count + 1:
        raise FederationError("V1 history does not contain discovery plus one failed aggregate attempt")


def _historical_exit_status(receipt: Mapping[str, Any]) -> str:
    qc = str(receipt["marker"]["qc"])
    return "0" if qc == "PROCESSED" else f"0:{qc}"


def _validate_v1_checkpoint(
    root: Path, entry: Mapping[str, Any], fingerprint: str, index: int,
) -> tuple[dict[str, str], dict[str, str], dict[str, object]]:
    required = {
        "locus_index", "locus", "chromosome", "qc", "bundle_path", "attempt_path",
        "ready", "receipt", "result", "log", "semantic_validation",
    }
    if set(entry) != required or entry.get("locus_index") != index:
        raise FederationError(f"source lock checkpoint {index} has a malformed identity")
    expected_bundle_rel = (
        Path("results/track_b/checkpoints/lava") / fingerprint
        / "discovery" / f"locus_{index:04d}"
    )
    if str(entry["bundle_path"]) != str(expected_bundle_rel):
        raise FederationError(f"source checkpoint {index} names a noncanonical READY bundle")
    attempt_rel = Path(str(entry["attempt_path"]))
    attempts_rel = expected_bundle_rel.parent / ".attempts"
    if attempt_rel.parent != attempts_rel or not attempt_rel.name.startswith(
        f"locus_{index:04d}."
    ):
        raise FederationError(f"source checkpoint {index} names a noncanonical attempt")
    attempt = safe_path(root, attempt_rel, f"source checkpoint {index} attempt")
    if attempt.is_symlink() or not attempt.is_dir():
        raise FederationError(f"source checkpoint {index} attempt is absent or unsafe")
    bundle_parent = safe_path(
        root, expected_bundle_rel.parent, f"source checkpoint {index} bundle parent",
    )
    bundle = bundle_parent / expected_bundle_rel.name
    if not bundle.is_symlink():
        raise FederationError(f"source checkpoint {index} lacks its atomic READY link")
    try:
        bundle_target = bundle.resolve(strict=True)
        attempt_target = attempt.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise FederationError(f"source checkpoint {index} READY link is unreadable") from error
    if bundle_target != attempt_target:
        raise FederationError(f"source checkpoint {index} READY link names a different attempt")
    expected_member_paths = {
        "ready": attempt_rel / "READY", "receipt": attempt_rel / "receipt.json",
        "result": attempt_rel / "result.rds", "log": attempt_rel / "worker.log",
        "semantic_validation": attempt_rel / "semantic_validation.txt",
    }
    for role, expected_path in expected_member_paths.items():
        record = entry[role]
        if not isinstance(record, Mapping) or str(record.get("path", "")) != str(expected_path):
            raise FederationError(f"source checkpoint {index} identity path family is noncanonical")
    paths: dict[str, Path] = {}
    identities: dict[str, dict[str, object]] = {}
    for role in ("ready", "receipt", "result", "log", "semantic_validation"):
        paths[role], identities[role] = _verify_identity_record(
            root, entry[role], f"source checkpoint {index} {role}",
        )
    ready, _ = read_json(paths["ready"], f"source checkpoint {index} READY")
    receipt, _ = read_json(paths["receipt"], f"source checkpoint {index} receipt")
    if ready != {
        "schema_version": 1, "state": "READY",
        "receipt_sha256": identities["receipt"]["sha256"],
    }:
        raise FederationError(f"source checkpoint {index} READY/receipt binding drifted")
    marker = receipt.get("marker")
    if (
        receipt.get("schema_version") != 3
        or receipt.get("analysis") != "LAVA_DISCOVERY"
        or receipt.get("phase") != "discovery"
        or receipt.get("locus_index") != index
        or receipt.get("execution_fingerprint") != fingerprint
        or receipt.get("exit_status") != 0
        or not isinstance(marker, Mapping)
        or marker.get("index") != str(index)
        or marker.get("locus") != str(entry["locus"])
        or marker.get("chromosome") != str(entry["chromosome"])
        or marker.get("qc") != str(entry["qc"])
    ):
        raise FederationError(f"source checkpoint {index} receipt identity drifted")
    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, list):
        raise FederationError(f"source checkpoint {index} has no artifact receipt")
    result = next(
        (item for item in artifacts if isinstance(item, Mapping) and item.get("path") == "result.rds"),
        None,
    )
    semantic = next(
        (
            item for item in artifacts
            if isinstance(item, Mapping) and item.get("path") == "semantic_validation.txt"
        ),
        None,
    )
    expected_result = {
        "path": "result.rds", "bytes": identities["result"]["bytes"],
        "sha256": identities["result"]["sha256"],
    }
    expected_semantic = {
        "path": "semantic_validation.txt", "bytes": identities["semantic_validation"]["bytes"],
        "sha256": identities["semantic_validation"]["sha256"],
    }
    if result != expected_result or semantic != expected_semantic:
        raise FederationError(f"source checkpoint {index} result/semantic receipt drifted")
    if (
        receipt.get("log_bytes") != identities["log"]["bytes"]
        or receipt.get("log_sha256") != identities["log"]["sha256"]
    ):
        raise FederationError(f"source checkpoint {index} worker log receipt drifted")
    pair = str(marker.get("pair", ""))
    n_snps = str(marker.get("n_snps", ""))
    if not pair or (n_snps != "NA" and not INTEGER_RE.fullmatch(n_snps)):
        raise FederationError(f"source checkpoint {index} marker lacks pair/SNP identity")
    row = {
        "analysis": "LAVA_DISCOVERY", "pair": pair,
        "locus": str(entry["locus"]), "chromosome": str(entry["chromosome"]),
        "n_snps": n_snps,
        "peak_ram_gb": f"{int(receipt['peak_rss_bytes']) / GIB:.6f}",
        "runtime_sec": f"{float(receipt['runtime_sec']):.3f}",
        "exit_status": _historical_exit_status(receipt),
        "output_hash": str(identities["result"]["sha256"]),
    }
    validate_federation_row(row, f"source checkpoint {index} reconstructed row")
    locus = {
        "locus": str(entry["locus"]), "chromosome": str(entry["chromosome"]),
        "n_snps": n_snps, "pair": pair,
    }
    evidence = {
        "index": index, "bundle_path": str(expected_bundle_rel),
        "attempt_path": str(attempt_rel),
        "peak_rss_bytes": int(receipt["peak_rss_bytes"]),
        "ready_sha256": identities["ready"]["sha256"],
        "receipt_sha256": identities["receipt"]["sha256"],
        "result_sha256": identities["result"]["sha256"],
        "log_sha256": identities["log"]["sha256"],
        "semantic_validation_sha256": identities["semantic_validation"]["sha256"],
    }
    return row, locus, evidence


def validate_v1_snapshot_rows(
    root: Path, benchmark_path: Path, history_receipt: Mapping[str, Any],
    source_lock: Mapping[str, Any], *, expected_count: int = EXPECTED_DISCOVERY_COUNT,
) -> dict[str, Any]:
    """Ignore stale provenance rows and rebuild every discovery row from raw receipts."""

    rows, benchmark_content = read_tsv(benchmark_path, BENCHMARK_FIELDS, "V1 history benchmark")
    validate_stale_v1_audit(
        history_receipt.get("historical_provenance_audit"), len(rows),
        expected_discovery_count=expected_count,
    )
    if (
        source_lock.get("source_discovery_fingerprint") != V1_FINGERPRINT
        or source_lock.get("checkpoint_count") != expected_count
        or not isinstance(source_lock.get("checkpoints"), list)
        or len(source_lock["checkpoints"]) != expected_count
        or source_lock.get("checkpoint_family_sha256")
        != digest_json(source_lock["checkpoints"])
    ):
        raise FederationError("source discovery lock count/fingerprint/family digest drifted")
    discovery_rows = [row for row in rows if row["analysis"] == "LAVA_DISCOVERY"]
    failed_rows = [row for row in rows if row["analysis"] != "LAVA_DISCOVERY"]
    if len(discovery_rows) != expected_count or len(failed_rows) != 1:
        raise FederationError("V1 history topology is not 2,495 discovery rows plus one failure")
    by_hash: dict[str, dict[str, str]] = {}
    for row in discovery_rows:
        validate_federation_row(row, "V1 history discovery row")
        if row["output_hash"] in by_hash:
            raise FederationError("V1 history repeats a source result hash")
        by_hash[row["output_hash"]] = row
    reconstructed: list[dict[str, str]] = []
    loci: dict[int, dict[str, str]] = {}
    members: list[dict[str, object]] = []
    for index, entry in enumerate(source_lock["checkpoints"], start=1):
        if not isinstance(entry, Mapping):
            raise FederationError(f"source lock checkpoint {index} is not an object")
        expected_row, locus, member = _validate_v1_checkpoint(root, entry, V1_FINGERPRINT, index)
        observed = by_hash.pop(expected_row["output_hash"], None)
        if observed != expected_row:
            raise FederationError(
                f"V1 historical discovery row {index} differs from its exact raw READY receipt/result"
            )
        reconstructed.append(expected_row)
        loci[index] = locus
        members.append(member)
    if by_hash:
        raise FederationError("V1 history contains discovery rows outside the source lock")

    failed = validate_federation_row(failed_rows[0], "V1 failed aggregate row")
    if (
        failed["analysis"] != "LAVA_DISCOVERY_AGGREGATION"
        or failed["pair"] != "UNKNOWN"
        or failed["locus"] != "ALL"
        or failed["chromosome"] != "ALL"
        or failed["n_snps"] != "NA"
        or failed["exit_status"] == "0"
    ):
        raise FederationError("V1 failed aggregate row identity drifted")
    failure_root = safe_path(
        root,
        Path("results/track_b/checkpoints/lava") / V1_FINGERPRINT / "failed_attempts",
        "V1 failed-attempt directory",
    )
    matches: list[dict[str, object]] = []
    if failure_root.is_dir() and not failure_root.is_symlink():
        for log in sorted(failure_root.glob("aggregate-discovery_all_*.log")):
            identity = portable_identity(root, log, "V1 failed aggregate log")
            if identity["sha256"] == failed["output_hash"]:
                text = stable_bytes(log, "V1 failed aggregate log").decode("utf-8")
                if (
                    "--phase aggregate-discovery" not in text
                    or f"--fingerprint {V1_FINGERPRINT}" not in text
                ):
                    raise FederationError("V1 failed aggregate log command identity drifted")
                matches.append(identity)
    if len(matches) != 1:
        raise FederationError("V1 failed aggregate benchmark row lacks one exact historical log")
    reconstructed.append(failed)
    return {
        "rows": reconstructed,
        # The retained historical failed aggregate predates exact-byte receipts.
        # Keep it, but never manufacture byte precision that does not exist.
        "peak_rss_bytes": [int(item["peak_rss_bytes"]) for item in members] + [None],
        "loci": loci,
        "checkpoint_family_revalidation_sha256": digest_json(members),
        "failed_aggregate_log": matches[0],
        "benchmark_identity": identity_for_bytes(
            benchmark_path.relative_to(root), benchmark_content,
        ),
        "stale_declared_rows": history_receipt["historical_provenance_audit"]["declared_attempt_rows"],
        "actual_rows": len(rows),
    }


def require_complete_continuation_family(
    aggregate: Mapping[str, Any] | None,
    conditionals: Sequence[Mapping[str, Any]],
    finalize: Mapping[str, Any] | None,
    terminal: Mapping[str, Any] | None,
    *, expected_count: int = EXPECTED_DISCOVERY_COUNT,
) -> str:
    if aggregate is None:
        raise FederationError("production c139 continuation lacks aggregate READY receipt")
    indices = [item.get("receipt", {}).get("locus_index") for item in conditionals]
    if len(conditionals) != expected_count or indices != list(range(1, expected_count + 1)):
        raise FederationError(
            f"production c139 continuation is incomplete: conditional={len(conditionals)}/{expected_count}"
        )
    if (finalize is None) == (terminal is None):
        raise FederationError(
            "production c139 continuation must have exactly one terminal success/failure receipt"
        )
    return "COMPLETE_PUBLISHED" if finalize is not None else "TERMINAL_FAILED_QC"


def _continuation_bundle_target(engine: Any, phase: str, index: int | None) -> tuple[Path, Path]:
    bundle = engine.bundle_path(phase, index, CONTINUATION_FINGERPRINT)
    if not bundle.is_symlink():
        raise FederationError(f"c139 {phase}/{index} is not an atomic READY link")
    try:
        target = bundle.resolve(strict=True)
        attempts = (engine.run_root(CONTINUATION_FINGERPRINT) / ".attempts").resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise FederationError(f"c139 {phase}/{index} READY link is unreadable") from error
    if target.parent != attempts or target.is_symlink() or not target.is_dir():
        raise FederationError(f"c139 {phase}/{index} READY link escapes its attempt namespace")
    return bundle, target


def _active_content_identity(
    engine: Any, item: Mapping[str, Any], cache: dict[str, tuple[int, str]], label: str,
) -> dict[str, object]:
    required = {"path", "bytes", "sha256", "role", "pre_stat", "post_stat"}
    if set(item) != required or item.get("pre_stat") != item.get("post_stat"):
        raise FederationError(f"{label} has a malformed active-input record")
    path_text = str(item["path"])
    path = engine._safe_relative(path_text)
    current = engine.stable_stat(path)
    if current != item["pre_stat"] or current["bytes"] != item["bytes"]:
        raise FederationError(f"{label} active input changed since execution: {path_text}")
    if path_text not in cache:
        content = stable_bytes(path, f"{label} active input")
        cache[path_text] = (len(content), hashlib.sha256(content).hexdigest())
    if cache[path_text] != (int(item["bytes"]), str(item["sha256"])):
        raise FederationError(f"{label} active input hash differs from receipt: {path_text}")
    return {key: item[key] for key in ("path", "bytes", "sha256", "role")}


def _precompute_continuation_invariants(
    root: Path, engine: Any, contract: Any, active_cache: dict[str, tuple[int, str]],
) -> dict[str, str | None]:
    """Hash invariant receipt dependencies once and seed the shared content cache."""

    def remember(path: Path, label: str) -> str:
        path = Path(path)
        if not path.is_absolute():
            path = root / path
        content = stable_bytes(path, label)
        try:
            relative = str(path.relative_to(root))
        except ValueError as error:
            raise FederationError(f"{label} is outside the repository: {path}") from error
        identity = (len(content), hashlib.sha256(content).hexdigest())
        prior = active_cache.setdefault(relative, identity)
        if prior != identity:
            raise FederationError(f"{label} has conflicting invariant identities")
        return identity[1]

    return {
        "lineage_sha256": remember(
            contract.lineage_path(CONTINUATION_FINGERPRINT), "c139 lineage invariant",
        ),
        "source_family_lock_sha256": remember(
            contract.source_lock_path(V1_FINGERPRINT), "c139 source-family lock invariant",
        ),
        "semantic_validator_sha256": remember(
            engine.CHECKPOINT_VALIDATOR, "c139 semantic-validator invariant",
        ),
        "terminal_qc_validator_sha256": remember(
            engine.RESULT_VALIDATOR, "c139 terminal-validator invariant",
        ),
    }


@contextlib.contextmanager
def _cache_continuation_legacy_loci(engine: Any):
    """Keep script 136's checks while validating its immutable legacy inputs once."""

    contract = engine.CONTRACT
    original_factory = contract.legacy_supervisor
    legacy = original_factory()
    loci = legacy.load_loci()
    baselines: dict[str, Any] = {}
    if len(loci) != EXPECTED_DISCOVERY_COUNT:
        raise FederationError("c139 official marker locus family is not exactly 2,495 loci")

    class CachedLegacy:
        def load_loci(self) -> Any:
            return loci

        def validate_input_baseline(self, fingerprint: str) -> Any:
            if fingerprint not in baselines:
                baselines[fingerprint] = legacy.validate_input_baseline(fingerprint)
            return baselines[fingerprint]

        def __getattr__(self, name: str) -> Any:
            return getattr(legacy, name)

    cached = CachedLegacy()
    contract.legacy_supervisor = lambda: cached
    try:
        yield loci
    finally:
        contract.legacy_supervisor = original_factory


def _validate_continuation_bundle(
    root: Path, engine: Any, contract: Any, phase: str, index: int | None,
    active_cache: dict[str, tuple[int, str]], invariant_hashes: Mapping[str, str | None],
) -> dict[str, Any] | None:
    bundle = engine.bundle_path(phase, index, CONTINUATION_FINGERPRINT)
    if not bundle.exists() and not bundle.is_symlink():
        return None
    _, target = _continuation_bundle_target(engine, phase, index)
    required_names = {"READY", "receipt.json", "worker.log", "result.rds", "semantic_validation.txt"}
    entries = list(target.iterdir())
    for item in entries:
        status = os.lstat(item)
        if not stat.S_ISREG(status.st_mode) or stat.S_ISLNK(status.st_mode):
            raise FederationError(f"c139 {phase}/{index} bundle contains unsafe content: {item.name}")
    live = {item.name for item in entries}
    if not required_names.issubset(live):
        raise FederationError(f"c139 {phase}/{index} bundle is incomplete")
    receipt, receipt_content = read_json(target / "receipt.json", f"c139 {phase}/{index} receipt")
    ready, ready_content = read_json(target / "READY", f"c139 {phase}/{index} READY")
    receipt_identity = identity_for_bytes(
        (target / "receipt.json").relative_to(root), receipt_content,
    )
    if ready != {
        "schema_version": engine.READY_SCHEMA, "state": "READY",
        "receipt_sha256": receipt_identity["sha256"],
    }:
        raise FederationError(f"c139 {phase}/{index} READY marker drifted")
    analysis = {
        "aggregate-discovery": "LAVA_DISCOVERY_AGGREGATION_V2",
        "conditional": "LAVA_CONDITIONAL_V2",
        "finalize": "LAVA_FINALIZATION_V2",
        "terminal-qc": "LAVA_TERMINAL_QC_V2",
    }[phase]
    expected_scalar = {
        "schema_version": engine.RECEIPT_SCHEMA,
        "phase": phase,
        "locus_index": index,
        "analysis": analysis,
        "requested_phase": "finalize" if phase == "terminal-qc" else phase,
        "source_discovery_fingerprint": V1_FINGERPRINT,
        "continuation_execution_fingerprint": CONTINUATION_FINGERPRINT,
        "exit_status": 78 if phase == "terminal-qc" else 0,
        "bundle_path": str(bundle.relative_to(root)),
        "lineage_sha256": invariant_hashes["lineage_sha256"],
        "source_family_lock_sha256": invariant_hashes["source_family_lock_sha256"],
        "semantic_validator_sha256": invariant_hashes["semantic_validator_sha256"],
        "terminal_qc_validator_sha256": (
            invariant_hashes["terminal_qc_validator_sha256"] if phase == "terminal-qc" else None
        ),
    }
    if any(receipt.get(key) != value for key, value in expected_scalar.items()):
        raise FederationError(f"c139 {phase}/{index} receipt scalar identity drifted")
    required_receipt = {
        "schema_version", "analysis", "phase", "requested_phase", "locus_index",
        "source_discovery_fingerprint", "continuation_execution_fingerprint", "marker",
        "command", "runtime_sec", "peak_rss_bytes", "exit_status", "bundle_path",
        "log_bytes", "log_sha256", "semantic_validator_sha256",
        "terminal_qc_validator_sha256", "lineage_sha256", "source_family_lock_sha256",
        "active_inputs", "semantic_artifact_binding", "artifacts",
    }
    if set(receipt) != required_receipt:
        raise FederationError(f"c139 {phase}/{index} receipt schema drifted")
    marker = receipt.get("marker")
    if not isinstance(marker, dict):
        raise FederationError(f"c139 {phase}/{index} receipt lacks marker")
    try:
        engine.validate_marker(
            marker, phase, index, target / "result.rds",
            V1_FINGERPRINT, CONTINUATION_FINGERPRINT,
        )
    except Exception as error:
        raise FederationError(f"c139 {phase}/{index} marker validation failed: {error}") from error
    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise FederationError(f"c139 {phase}/{index} lacks artifact receipts")
    artifact_records: dict[str, dict[str, object]] = {}
    artifact_identities: dict[str, dict[str, object]] = {}
    for artifact in artifacts:
        if not isinstance(artifact, Mapping) or set(artifact) != {"path", "bytes", "sha256"}:
            raise FederationError(f"c139 {phase}/{index} has malformed artifact receipt")
        name = str(artifact["path"])
        if Path(name).name != name or name.startswith(".") or name in artifact_records:
            raise FederationError(f"c139 {phase}/{index} has unsafe/duplicate artifact")
        identity = portable_identity(root, target / name, f"c139 {phase}/{index} artifact {name}")
        expected = {"path": name, "bytes": identity["bytes"], "sha256": identity["sha256"]}
        if dict(artifact) != expected:
            raise FederationError(f"c139 {phase}/{index} artifact differs from receipt: {name}")
        artifact_records[name] = expected
        artifact_identities[name] = identity
    if live != set(artifact_records) | {"READY", "receipt.json", "worker.log"}:
        raise FederationError(f"c139 {phase}/{index} contains unreceipted files")
    if "result.rds" not in artifact_records or "semantic_validation.txt" not in artifact_records:
        raise FederationError(f"c139 {phase}/{index} lacks result/semantic evidence")
    worker_identity = portable_identity(root, target / "worker.log", f"c139 {phase}/{index} log")
    if (
        receipt.get("log_bytes") != worker_identity["bytes"]
        or receipt.get("log_sha256") != worker_identity["sha256"]
    ):
        raise FederationError(f"c139 {phase}/{index} worker log differs from receipt")
    binding = receipt.get("semantic_artifact_binding")
    try:
        snapshot = engine.semantic_snapshot_from_binding(binding, phase)
    except Exception as error:
        raise FederationError(f"c139 {phase}/{index} semantic binding is invalid: {error}") from error
    for item in snapshot:
        name = str(item["path"])
        if artifact_records.get(name) != {
            "path": name, "bytes": item["bytes"], "sha256": item["sha256"],
        } or engine.stable_stat(target / name) != item["stable_stat"]:
            raise FederationError(f"c139 {phase}/{index} semantic artifact binding drifted: {name}")
    active = receipt.get("active_inputs")
    if not isinstance(active, list) or not active:
        raise FederationError(f"c139 {phase}/{index} lacks active-input family")
    observed_specs: list[dict[str, object]] = []
    observed_paths: set[str] = set()
    for item in active:
        if not isinstance(item, Mapping):
            raise FederationError(f"c139 {phase}/{index} active input is malformed")
        spec = _active_content_identity(
            engine, item, active_cache, f"c139 {phase}/{index}",
        )
        if str(spec["path"]) in observed_paths:
            raise FederationError(f"c139 {phase}/{index} repeats an active input")
        observed_paths.add(str(spec["path"]))
        observed_specs.append(spec)
    return {
        "receipt": receipt, "target": target,
        "ready_identity": identity_for_bytes(
            (target / "READY").relative_to(root), ready_content,
        ),
        "receipt_identity": receipt_identity,
        "result_identity": artifact_identities["result.rds"],
        "artifact_records": artifact_records,
        "artifact_identities": artifact_identities,
        "semantic_snapshot": snapshot,
        "active_specs": sorted(observed_specs, key=lambda item: str(item["path"])),
    }


def _spec(path: Path, role: str, root: Path) -> dict[str, object]:
    identity = portable_identity(root, path, f"continuation expected {role}")
    return {**identity, "role": role}


def _source_checkpoint_specs(source_lock: Mapping[str, Any]) -> list[dict[str, object]]:
    """Materialize the 2,495 source-result specs once, not once per conditional."""

    return [
        {
            "path": str(entry["result"]["path"]), "bytes": int(entry["result"]["bytes"]),
            "sha256": str(entry["result"]["sha256"]), "role": "SOURCE_DISCOVERY_CHECKPOINT",
        }
        for entry in source_lock["checkpoints"]
    ]


def _expected_continuation_active_specs(
    root: Path, engine: Any, contract: Any,
    aggregate: Mapping[str, Any], conditionals: Sequence[Mapping[str, Any]],
    phase: str, index: int | None, candidate_loci: frozenset[int],
    cached: dict[tuple[str, int | None], list[dict[str, object]]],
    source_specs: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    base_key = (phase, None if phase != "conditional" else (index if index in candidate_loci else 0))
    if base_key not in cached:
        base = [
            _spec(contract.lineage_path(CONTINUATION_FINGERPRINT), "CONTINUATION_LINEAGE", root),
            _spec(contract.source_lock_path(V1_FINGERPRINT), "SOURCE_FAMILY_LOCK", root),
            _spec(contract.predecessor_supersession_lock_path(), "PREDECESSOR_SUPERSESSION_LOCK", root),
        ]
        try:
            base.extend(engine._execution_code_specs())
            base.extend(engine._scientific_input_specs(
                phase, index, V1_FINGERPRINT,
                include_conditional_chromosome=(
                    phase == "conditional" and index is not None and index in candidate_loci
                ),
            ))
        except Exception as error:
            raise FederationError(f"could not reconstruct c139 active-input contract: {error}") from error
        cached[base_key] = base
    specs = list(cached[base_key])
    aggregate_specs = []
    for name, role in (
        ("result.rds", "CONTINUATION_CHECKPOINT"),
        ("conditional_candidates.tsv", "CONDITIONAL_CANDIDATE_ATTESTATION"),
    ):
        record = aggregate["artifact_identities"][name]
        aggregate_specs.append({**record, "role": role})
    aggregate_specs.append({**aggregate["receipt_identity"], "role": "CONTINUATION_CHECKPOINT_RECEIPT"})
    if phase == "aggregate-discovery":
        specs.extend(dict(item) for item in source_specs)
    elif phase == "conditional":
        if index is None:
            raise FederationError("conditional expected-spec reconstruction lacks index")
        specs.append(dict(source_specs[index - 1]))
        specs.extend(aggregate_specs)
    else:
        specs.extend(dict(item) for item in source_specs)
        specs.extend(aggregate_specs)
        specs.extend(
            {**item["result_identity"], "role": "CONTINUATION_CHECKPOINT"}
            for item in conditionals
        )
    return sorted(specs, key=lambda item: str(item["path"]))


def _continuation_raw_benchmark_row(bundle: Mapping[str, Any], fields: Sequence[str]) -> dict[str, str]:
    receipt = bundle["receipt"]
    marker = receipt["marker"]
    result_hash = bundle["result_identity"]["sha256"]
    row = {
        "source_discovery_fingerprint": V1_FINGERPRINT,
        "continuation_execution_fingerprint": CONTINUATION_FINGERPRINT,
        "analysis": str(receipt["analysis"]), "phase": str(receipt["phase"]),
        "locus_index": "NA" if receipt["locus_index"] is None else str(receipt["locus_index"]),
        "locus": str(marker["locus"]), "chromosome": str(marker["chromosome"]),
        "pair": str(marker["pair"]), "qc": str(marker["qc"]),
        "peak_ram_gb": f"{int(receipt['peak_rss_bytes']) / GIB:.6f}",
        "runtime_sec": f"{float(receipt['runtime_sec']):.3f}",
        "exit_status": str(receipt["exit_status"]), "output_sha256": str(result_hash),
    }
    if set(row) != set(fields):
        raise FederationError("c139 benchmark schema differs from its receipt reconstruction")
    return {field: row[field] for field in fields}


def _require_no_unbenchmarkable_failed_attempts(root: Path, failure_root: Path) -> None:
    """Refuse a completeness claim when script 136 retained only unmetered failure logs."""

    failed_attempts: list[dict[str, object]] = []
    if failure_root.exists() or failure_root.is_symlink():
        if failure_root.is_symlink() or not failure_root.is_dir():
            raise FederationError("c139 failed-attempt namespace is unsafe")
        for path in sorted(failure_root.iterdir()):
            if not path.is_file() or path.is_symlink():
                raise FederationError("c139 failed-attempt namespace contains unsafe content")
            failed_attempts.append(portable_identity(root, path, "c139 failed attempt"))
    if failed_attempts:
        raise FederationError(
            "c139 contains failed attempts without benchmarkable immutable resource receipts; "
            "publishing a table that claims full failure retention is forbidden"
        )


def _continuation_federation_row(
    bundle: Mapping[str, Any], source_loci: Mapping[int, Mapping[str, str]],
) -> dict[str, str]:
    receipt = bundle["receipt"]
    marker = receipt["marker"]
    phase = str(receipt["phase"])
    index = receipt["locus_index"]
    if phase == "conditional":
        source = source_loci.get(int(index))
        if source is None:
            raise FederationError(f"c139 conditional locus {index} lacks exact source discovery locus")
        if (
            str(marker["locus"]) != source["locus"]
            or str(marker["chromosome"]) != source["chromosome"]
        ):
            raise FederationError(f"c139 conditional locus {index} differs from source discovery")
        n_snps = source["n_snps"]
        if str(marker.get("n_snps")) not in {n_snps, "NA"}:
            raise FederationError(f"c139 conditional locus {index} SNP count contradicts source")
    else:
        n_snps = str(marker.get("n_snps", "NA"))
        if n_snps != "NA" and not INTEGER_RE.fullmatch(n_snps):
            n_snps = "NA"
    qc = str(marker["qc"])
    success_qc = {
        "aggregate-discovery": {"FULL_FAMILY_BH_COMPLETE"},
        "conditional": {"CONDITIONAL_LOCUS_PROCESSED"},
        "finalize": {"STAGED_SEMANTIC_VALIDATION_COMPLETE"},
        "terminal-qc": set(),
    }[phase]
    code = str(receipt["exit_status"])
    exit_status = code if qc in success_qc else f"{code}:{qc}"
    row = {
        "analysis": str(receipt["analysis"]), "pair": str(marker["pair"]),
        "locus": str(marker["locus"]), "chromosome": str(marker["chromosome"]),
        "n_snps": n_snps,
        "peak_ram_gb": f"{int(receipt['peak_rss_bytes']) / GIB:.6f}",
        "runtime_sec": f"{float(receipt['runtime_sec']):.3f}",
        "exit_status": exit_status,
        "output_hash": str(bundle["result_identity"]["sha256"]),
    }
    return validate_federation_row(row, f"c139 {phase}/{index} federation row")


def _validate_terminal_qc_attestation(
    root: Path, terminal: Mapping[str, Any], aggregate: Mapping[str, Any],
    conditionals: Sequence[Mapping[str, Any]], source_lock: Mapping[str, Any],
    contract: Any, results: Any,
) -> dict[str, Any]:
    path = terminal["target"] / "terminal_qc.json"
    payload, _ = read_json(path, "c139 terminal-QC attestation")
    try:
        expected_payload = results._terminal_payload(
            terminal["target"], V1_FINGERPRINT, CONTINUATION_FINGERPRINT,
        )
    except Exception as error:
        raise FederationError(f"c139 terminal-QC deep verifier failed: {error}") from error
    if payload != expected_payload:
        raise FederationError("c139 terminal-QC attestation differs from its deep reconstruction")
    compact = [
        {
            "locus_index": item["receipt"]["locus_index"],
            "receipt_sha256": item["receipt_identity"]["sha256"],
            "result_sha256": item["result_identity"]["sha256"],
            "qc": str(item["receipt"]["marker"]["qc"]),
        }
        for item in conditionals
    ]
    aggregate_summary = {
        "phase": "aggregate-discovery", "locus_index": None,
        "receipt": aggregate["receipt_identity"],
        "result": aggregate["result_identity"],
        "qc": str(aggregate["receipt"]["marker"]["qc"]),
    }
    upstream = {
        "aggregate": aggregate_summary,
        "conditional_checkpoint_count": len(conditionals),
        "conditional_checkpoint_family_sha256": digest_json(compact),
        "complete": True,
    }
    qc = payload.get("qc")
    if not isinstance(qc, Mapping) or not isinstance(qc.get("reasons"), list) or not qc["reasons"]:
        raise FederationError("c139 terminal-QC attestation lacks explicit frozen-gate failures")
    expected_scalar = {
        "schema_version": results.TERMINAL_QC_SCHEMA,
        "state": "TERMINAL_FAILED_QC", "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": V1_FINGERPRINT,
        "continuation_execution_fingerprint": CONTINUATION_FINGERPRINT,
        "full_family_complete": True, "scientific_validation_passed": False,
        "canonical_publication_allowed": False, "terminal": True,
        "upstream_continuation_family": upstream,
        "source_family_lock": portable_identity(
            root, contract.source_lock_path(V1_FINGERPRINT), "c139 source-family lock",
        ),
        "source_checkpoint_family_sha256": source_lock["checkpoint_family_sha256"],
        "lineage": portable_identity(
            root, contract.lineage_path(CONTINUATION_FINGERPRINT), "c139 lineage",
        ),
    }
    if any(payload.get(key) != value for key, value in expected_scalar.items()):
        raise FederationError("c139 terminal-QC attestation/family binding drifted")
    validators = payload.get("validators")
    expected_validators = {
        "legacy_scientific_validator": portable_identity(
            root, results.LEGACY_VALIDATOR_PATH, "legacy scientific validator",
        ),
        "continuation_result_validator": portable_identity(
            root, root / CONTINUATION_RESULTS_REL,
            "continuation result validator",
        ),
    }
    if validators != expected_validators:
        raise FederationError("c139 terminal-QC validator identities drifted")
    staged = payload.get("staged_low_level_results")
    if not isinstance(staged, list) or len(staged) != 4:
        raise FederationError("c139 terminal-QC staged result identity family drifted")
    for record in staged:
        _verify_identity_record(root, record, "c139 terminal-QC staged result")
    return {
        "state": "TERMINAL_FAILED_QC",
        "attestation": portable_identity(root, path, "c139 terminal-QC attestation"),
        "qc": qc,
    }


def validate_continuation_component(
    root: Path, source_lock: Mapping[str, Any], source_loci: Mapping[int, Mapping[str, str]],
    *, expected_count: int = EXPECTED_DISCOVERY_COUNT,
) -> dict[str, Any]:
    engine = _load_module("track_b_ram_federation_runner_136", root / CONTINUATION_RUNNER_REL)
    # Use script 136's exact loaded contract instance; importing 134 a second
    # time would repeat its source-lock and official-locus reconstruction.
    contract = engine.CONTRACT
    results = _load_module("track_b_ram_federation_results_138", root / CONTINUATION_RESULTS_REL)
    failure_root = engine.run_root(CONTINUATION_FINGERPRINT) / "failed_attempts"
    _require_no_unbenchmarkable_failed_attempts(root, failure_root)
    active_cache: dict[str, tuple[int, str]] = {}
    with _cache_continuation_legacy_loci(engine):
        try:
            contract.validate_lineage(V1_FINGERPRINT, CONTINUATION_FINGERPRINT, deep=False)
            engine.quick_contract(V1_FINGERPRINT, CONTINUATION_FINGERPRINT)
        except Exception as error:
            raise FederationError(
                f"production c139 lineage/contract failed verification: {error}"
            ) from error
        invariant_hashes = _precompute_continuation_invariants(
            root, engine, contract, active_cache,
        )
        aggregate = _validate_continuation_bundle(
            root, engine, contract, "aggregate-discovery", None, active_cache,
            invariant_hashes,
        )
        conditionals: list[dict[str, Any]] = []
        for index in range(1, expected_count + 1):
            member = _validate_continuation_bundle(
                root, engine, contract, "conditional", index, active_cache,
                invariant_hashes,
            )
            if member is not None:
                conditionals.append(member)
        finalize = _validate_continuation_bundle(
            root, engine, contract, "finalize", None, active_cache, invariant_hashes,
        )
        terminal = _validate_continuation_bundle(
            root, engine, contract, "terminal-qc", None, active_cache, invariant_hashes,
        )
        state = require_complete_continuation_family(
            aggregate, conditionals, finalize, terminal, expected_count=expected_count,
        )
        if aggregate is None:  # guarded by require_complete_continuation_family
            raise FederationError("production c139 aggregate vanished during verification")
        try:
            candidate_loci = engine._candidate_loci(
                aggregate["target"] / "conditional_candidates.tsv"
            )
        except Exception as error:
            raise FederationError(f"c139 candidate attestation failed validation: {error}") from error
        expected_cache: dict[tuple[str, int | None], list[dict[str, object]]] = {}
        source_specs = _source_checkpoint_specs(source_lock)
        all_bundles: list[dict[str, Any]] = [aggregate, *conditionals]
        all_bundles.append(finalize if finalize is not None else terminal)  # type: ignore[arg-type]
        for member in all_bundles:
            phase = str(member["receipt"]["phase"])
            index = member["receipt"]["locus_index"]
            expected = _expected_continuation_active_specs(
                root, engine, contract, aggregate, conditionals,
                phase, index, candidate_loci, expected_cache, source_specs,
            )
            if member["active_specs"] != expected:
                raise FederationError(f"c139 {phase}/{index} active-input family drifted")

    terminal_evidence: dict[str, Any]
    if finalize is not None:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                results.verify_results(V1_FINGERPRINT, CONTINUATION_FINGERPRINT)
        except Exception as error:
            raise FederationError(f"c139 canonical terminal result verification failed: {error}") from error
        result_lock = results.LEGACY.contract.RESULT_LOCK
        terminal_evidence = {
            "state": state,
            "result_lock": portable_identity(root, result_lock, "c139 terminal result lock"),
        }
    else:
        if terminal is None:  # guarded by require_complete_continuation_family
            raise FederationError("production c139 terminal receipt vanished during verification")
        terminal_evidence = _validate_terminal_qc_attestation(
            root, terminal, aggregate, conditionals, source_lock, contract, results,
        )

    benchmark_path, namespace_path, provenance_path = engine.benchmark_paths(
        CONTINUATION_FINGERPRINT
    )
    raw_rows, benchmark_content = read_tsv(
        benchmark_path, engine.BENCHMARK_FIELDS, "c139 RAM benchmark",
    )
    expected_raw = [
        _continuation_raw_benchmark_row(member, engine.BENCHMARK_FIELDS)
        for member in all_bundles
    ]
    expected_raw.sort(key=lambda row: (
        row["phase"], int(row["locus_index"]) if row["locus_index"] != "NA" else 0,
    ))
    if raw_rows != expected_raw:
        raise FederationError("c139 RAM benchmark differs from the complete READY receipt family")
    namespace, namespace_content = read_json(namespace_path, "c139 RAM namespace")
    expected_namespace = {
        "schema_version": "track-b-lava-continuation-ram-namespace.1",
        "analysis_id": "track-b-v1.0-local",
        "source_discovery_fingerprint": V1_FINGERPRINT,
        "continuation_execution_fingerprint": CONTINUATION_FINGERPRINT,
    }
    if namespace != expected_namespace:
        raise FederationError("c139 RAM namespace drifted")
    provenance, provenance_content = read_json(provenance_path, "c139 RAM provenance")
    expected_provenance = {
        "schema_version": "track-b-lava-continuation-ram-provenance.1",
        "source_discovery_fingerprint": V1_FINGERPRINT,
        "continuation_execution_fingerprint": CONTINUATION_FINGERPRINT,
        "namespace_sha256": hashlib.sha256(namespace_content).hexdigest(),
        "row_count": len(raw_rows),
        "table_sha256": hashlib.sha256(benchmark_content).hexdigest(),
    }
    if provenance != expected_provenance:
        raise FederationError("c139 RAM provenance is stale or inconsistent")
    federation_rows = [
        _continuation_federation_row(member, source_loci) for member in all_bundles
    ]
    compact = [
        {
            "index": member["receipt"]["locus_index"],
            "ready_sha256": member["ready_identity"]["sha256"],
            "receipt_sha256": member["receipt_identity"]["sha256"],
            "result_sha256": member["result_identity"]["sha256"],
        }
        for member in conditionals
    ]
    return {
        "rows": federation_rows,
        "peak_rss_bytes": [int(member["receipt"]["peak_rss_bytes"]) for member in all_bundles],
        "execution_fingerprint": CONTINUATION_FINGERPRINT,
        "evidence": {
            "source_fingerprint": V1_FINGERPRINT,
            "continuation_fingerprint": CONTINUATION_FINGERPRINT,
            "lineage": portable_identity(
                root, contract.lineage_path(CONTINUATION_FINGERPRINT), "c139 lineage",
            ),
            "aggregate_ready": aggregate["ready_identity"],
            "aggregate_receipt": aggregate["receipt_identity"],
            "aggregate_result": aggregate["result_identity"],
            "conditional_count": len(conditionals),
            "conditional_ready_receipt_result_family_sha256": digest_json(compact),
            "terminal_ready": all_bundles[-1]["ready_identity"],
            "terminal_receipt": all_bundles[-1]["receipt_identity"],
            "terminal_result": all_bundles[-1]["result_identity"],
            "terminal": terminal_evidence,
            "benchmark": identity_for_bytes(benchmark_path.relative_to(root), benchmark_content),
            "namespace": identity_for_bytes(namespace_path.relative_to(root), namespace_content),
            "provenance": identity_for_bytes(provenance_path.relative_to(root), provenance_content),
            "failed_attempts": [],
        },
    }


def validate_superseded_exclusions(root: Path) -> dict[str, dict[str, object]]:
    module = _load_module("track_b_ram_federation_supersession_144", root / SUPERSESSION_REL)
    output: dict[str, dict[str, object]] = {}
    for fingerprint in SUPERSEDED_CONTINUATIONS:
        try:
            payload = module.verify(fingerprint)
        except Exception as error:
            raise FederationError(f"superseded continuation {fingerprint} failed lock verification: {error}") from error
        if payload.get("superseded_continuation_fingerprint") != fingerprint:
            raise FederationError("supersession lock names the wrong continuation")
        output[fingerprint] = portable_identity(
            root, module.lock_path(fingerprint), f"supersession lock {fingerprint}",
        )
    return output


def require_placo_terminal_family(
    completions: Mapping[str, object],
) -> None:
    missing = [pair for pair in EXPECTED_PLACO_PAIR_ORDER if pair not in completions]
    extra = sorted(set(completions) - set(EXPECTED_PLACO_PAIR_ORDER))
    if missing or extra:
        raise FederationError(
            "PLACO full-family cleanup is incomplete "
            f"(missing={','.join(missing) or 'NONE'} extra={','.join(extra) or 'NONE'})"
        )


def _placo_completion_expected(module: Any, context: Mapping[str, Any], bundle: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": module.CLEANUP_COMPLETE_SCHEMA,
        "analysis_id": context["policy"]["analysis_id"],
        "pair_id": context["pair_id"],
        "run_fingerprint": context["fingerprint"],
        "cleanup_bundle_sha256": module.stable_identity(module.cleanup_bundle_path(context))["sha256"],
        "removed_artifact_count": len(bundle["cleanup_candidates"]),
        "removed_logical_bytes": sum(int(item["bytes"]) for item in bundle["cleanup_candidates"]),
        "reclaim_accounting_at_authorization": bundle["cleanup_reclaim_accounting"],
        "canonical_ledger_retained": bundle["canonical"]["ledger"],
        "canonical_provenance_retained": bundle["canonical"]["provenance"],
        "source_or_canonical_deletion": False,
        "ld_extracted": False,
        "locus_claim_status": "BLOCKED_PENDING_FULL_P_FAMILY_AND_SEPARATE_LD_COLLATOR",
        "coordinator_sha256": module.coordinator_sha256(),
    }


def validate_placo_execution_contract(module: Any, root: Path) -> dict[str, Any]:
    """Validate script 143's exact sealed contract without its legacy O(N*locus) gate walk.

    ``script 143.verify_execution_contract`` rebuilds the upstream terminal gate,
    whose historical implementation rehashes the multi-gigabyte V1 baseline
    separately for every discovery locus.  The final federation has already
    revalidated every member of that immutable V1 source lock.  Here we still
    use script 143's own payload constructor: the exact sealed terminal-gate
    identity, reference-provenance identity, all code identities, constants,
    RAM rules, pair order, cleanup rules, and every other contract field must
    reconstruct byte-for-byte.
    """

    if tuple(getattr(module, "PAIR_ORDER", ())) != EXPECTED_PLACO_PAIR_ORDER:
        raise FederationError("PLACO module pair order is not exactly A/B/CONTROL")
    if tuple(getattr(module.ENGINE, "PAIR_IDS", ())) != EXPECTED_PLACO_PAIR_ORDER:
        raise FederationError("PLACO engine pair scope is not exactly A/B/CONTROL")
    contract_path = root / module.EXECUTION_CONTRACT
    observed = module.read_json(contract_path)
    if observed.get("pair_order") != list(EXPECTED_PLACO_PAIR_ORDER):
        raise FederationError("PLACO sealed contract pair order is not exactly A/B/CONTROL")
    guard = observed.get("lava_reference_archive_state_guard")
    required_guard = {
        "accepted_states", "guard_implementation", "reference_provenance",
        "archive_family_sha256", "extracted_payload_count",
        "extracted_family_sha256", "verification",
    }
    if not isinstance(guard, Mapping) or set(guard) != required_guard:
        raise FederationError("PLACO sealed contract has a malformed reference-state invariant")
    reference_record = guard.get("reference_provenance")
    _verify_identity_record(root, reference_record, "PLACO sealed reference provenance")
    gate_lock = module.read_json(root / module.GATE.TERMINAL_GATE_LOCK)
    readiness_rel = module.GATE.READINESS_GATE
    if gate_lock.get("readiness_gate") != str(readiness_rel):
        raise FederationError("PLACO terminal gate names the wrong readiness artifact")
    readiness = stable_bytes(
        safe_path(root, readiness_rel, "PLACO readiness gate"), "PLACO readiness gate",
    )
    if hashlib.sha256(readiness).hexdigest() != gate_lock.get("readiness_gate_sha256"):
        raise FederationError("PLACO readiness gate differs from its terminal lock")
    try:
        reference_state = module.verify_reference_archive_state(gate_lock)
    except Exception as error:
        raise FederationError(f"PLACO reference archive-state invariant failed: {error}") from error
    expected = module.execution_contract_payload(gate_lock, reference_state)
    if observed != expected:
        raise FederationError("PLACO sequential execution contract is missing or drifted")
    return observed


def derive_placo_context_from_verified_contract(
    module: Any, root: Path, pair: str, gate_lock: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Derive script 140's pair context from the exact gate bound by the contract."""

    bridge = module.BRIDGE
    engine = module.ENGINE
    gate = module.GATE
    evidence = gate_lock.get("lava_scientific_evidence", {})
    if (
        not isinstance(evidence, Mapping)
        or evidence.get("status") != "FAILED_QC_NOT_CONSUMED"
        or evidence.get("scientific_validation_passed") is not False
        or evidence.get("passed_lava_results_used_by_placo") is not False
        or evidence.get("consumption") != "ORDERING_TERMINALITY_ATTESTATION_ONLY"
    ):
        raise FederationError("PLACO sealed gate lost FAILED_QC_NOT_CONSUMED semantics")
    if gate_lock.get("pair_scope") != ["A", "B", "CONTROL"] or pair not in engine.PAIR_IDS:
        raise FederationError("PLACO sealed gate/context scope is not exactly A/B/CONTROL")

    original_assert = engine.assert_upstream_ready
    original_input_lock = engine.INPUT_LOCK
    original_runner = engine.RUNNER
    engine.assert_upstream_ready = bridge.inherited_v1_upstream
    engine.INPUT_LOCK = gate.V1_INPUT_LOCK
    engine.RUNNER = bridge.V2_RUNNER
    try:
        base = engine.production_context(root, pair)
    except SystemExit as error:
        engine.RUNNER = original_runner
        raise FederationError(str(error)) from error
    except BaseException:
        engine.RUNNER = original_runner
        raise
    finally:
        engine.assert_upstream_ready = original_assert
        engine.INPUT_LOCK = original_input_lock

    # Match the terminal V2 context used by script 140 after its V1-compatible
    # source discovery.  The execution contract has already revalidated the
    # current identities of script 140 and every code/input lock named below.
    engine.INPUT_LOCK = gate.TERMINAL_GATE_LOCK
    gate_hash = engine.sha256(root / gate.TERMINAL_GATE_LOCK)
    payload = {
        "schema_version": "sleep-atlas-track-b-placo-execution-fingerprint.2",
        "analysis_id": base["policy"]["analysis_id"],
        "pair_id": pair,
        "traits": [base["trait1"], base["trait2"]],
        "family_role": base["family_role"],
        "policy_sha256": engine.sha256(root / engine.POLICY),
        "v1_contract_lock_sha256": engine.sha256(root / engine.CONTRACT_LOCK),
        "v1_input_gate_lock_sha256": engine.sha256(root / gate.V1_INPUT_LOCK),
        "v2_terminal_gate_lock_sha256": gate_hash,
        "v1_materialization_validation_engine_sha256": engine.sha256(bridge.ENGINE_PATH),
        "v2_materializer_bridge_sha256": engine.sha256(Path(bridge.__file__).resolve()),
        "v1_runner_sha256": engine.sha256(root / bridge.V1_RUNNER),
        "v2_worker_bound_runner_sha256": engine.sha256(root / bridge.V2_RUNNER),
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
    fingerprint = hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
    context = {
        **base, "input_lock": dict(gate_lock), "fingerprint": fingerprint,
        "fingerprint_payload": payload,
    }
    try:
        bridge.validate_v2_paths(context)
    except Exception as error:
        raise FederationError(f"PLACO {pair} V2 path validation failed: {error}") from error
    return context, dict(gate_lock)


def validate_placo_component(root: Path) -> dict[str, Any]:
    module = _load_module("track_b_ram_federation_placo_143", root / PLACO_REL)
    try:
        contract = validate_placo_execution_contract(module, root)
    except Exception as error:
        raise FederationError(f"PLACO sequential contract failed verification: {error}") from error
    gate_lock = module.read_json(root / module.GATE.TERMINAL_GATE_LOCK)
    contexts: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {}
    completions: dict[str, dict[str, Any]] = {}
    cleanup_identities: dict[str, dict[str, object]] = {}
    cleanup_bundle_identities: dict[str, dict[str, object]] = {}
    projections: dict[str, dict[str, str]] = {}
    for pair in EXPECTED_PLACO_PAIR_ORDER:
        try:
            context, gate = derive_placo_context_from_verified_contract(
                module, root, pair, gate_lock,
            )
            bundle = module.validate_cleanup_bundle(context, contract)
            completion_path = module.cleanup_complete_path(context)
            completion = module.read_json(completion_path)
        except Exception as error:
            raise FederationError(f"PLACO {pair} terminal cleanup verification failed: {error}") from error
        if completion != _placo_completion_expected(module, context, bundle):
            raise FederationError(f"PLACO {pair} cleanup completion receipt drifted")
        for record in bundle["cleanup_candidates"]:
            candidate = root / str(record["path"])
            if candidate.exists() or candidate.is_symlink():
                raise FederationError(f"PLACO {pair} cleanup is not terminal")
        work_root = module.ENGINE.materialized_paths(context)["base"]
        if work_root.exists() or work_root.is_symlink():
            raise FederationError(f"PLACO {pair} retired work root still exists")
        contexts[pair] = (context, gate)
        completions[pair] = completion
        cleanup_identities[pair] = portable_identity(
            root, completion_path, f"PLACO {pair} cleanup completion",
        )
        cleanup_bundle_identities[pair] = portable_identity(
            root, module.cleanup_bundle_path(context), f"PLACO {pair} cleanup bundle",
        )
        try:
            embedded = base64.b64decode(
                bundle["embedded_reproducibility_artifacts"]["benchmark.tsv"]["base64"],
                validate=True,
            ).decode("utf-8")
            row = next(csv.DictReader(io.StringIO(embedded), delimiter="\t"))
            projections[pair] = {
                "projected_full_family_seconds": row["projected_full_family_seconds"],
                "measured_runner_elapsed_seconds": row["measured_runner_elapsed_seconds"],
                "global_nuisance_elapsed_seconds": row["global_nuisance_elapsed_seconds"],
            }
        except (KeyError, ValueError, UnicodeError, StopIteration, csv.Error) as error:
            raise FederationError(f"PLACO {pair} embedded runtime projection is unreadable") from error
    require_placo_terminal_family(completions)

    record_root = root / module.RECORD_ROOT
    rows_with_paths: list[tuple[dict[str, str], Path, Path, int]] = []
    fingerprints_by_pair: dict[str, str] = {}
    if record_root.is_symlink() or not record_root.is_dir():
        raise FederationError("PLACO immutable measurement namespace is absent or unsafe")
    for fingerprint_root in sorted(record_root.iterdir()):
        if (
            fingerprint_root.is_symlink() or not fingerprint_root.is_dir()
            or not SHA256_RE.fullmatch(fingerprint_root.name)
        ):
            raise FederationError("PLACO measurement fingerprint namespace is unsafe")
        pair_dirs = sorted(
            path for path in fingerprint_root.iterdir()
            if path.name in EXPECTED_PLACO_PAIR_ORDER
        )
        unknown = [
            path.name for path in fingerprint_root.iterdir()
            if path.name not in EXPECTED_PLACO_PAIR_ORDER
        ]
        if unknown or len(pair_dirs) != 1:
            raise FederationError("PLACO measurement fingerprint must contain exactly one pair")
        pair = pair_dirs[0].name
        if pair in fingerprints_by_pair:
            raise FederationError(f"PLACO measurement namespace repeats pair {pair}")
        fingerprints_by_pair[pair] = fingerprint_root.name
        ram_root = pair_dirs[0] / "ram"
        if ram_root.is_symlink() or not ram_root.is_dir():
            raise FederationError(f"PLACO {pair} RAM evidence namespace is absent or unsafe")
        entries = sorted(ram_root.iterdir())
        if any(
            path.is_symlink() or not path.is_file()
            or not (path.name.endswith(".tsv") or path.name.endswith(".provenance.json"))
            for path in entries
        ):
            raise FederationError(f"PLACO {pair} RAM evidence contains an unsafe artifact")
        provenance_paths = [path for path in entries if path.name.endswith(".provenance.json")]
        for provenance_path in provenance_paths:
            phase = provenance_path.name.removesuffix(".provenance.json")
            receipt = provenance_path.with_name(f"{phase}.tsv")
            if not receipt.is_file() or receipt.is_symlink():
                raise FederationError(f"PLACO {pair} measurement provenance lacks its TSV")
            try:
                provenance = module.read_json(provenance_path)
                row, expected_content = module.validate_ram_measurement_provenance(
                    provenance, receipt, context=contexts[pair][0], require_live_task=False,
                )
            except Exception as error:
                raise FederationError(f"PLACO {pair}/{phase} RAM provenance failed: {error}") from error
            if stable_bytes(receipt, f"PLACO {pair}/{phase} RAM TSV") != expected_content:
                raise FederationError(f"PLACO {pair}/{phase} RAM TSV differs from provenance")
            rows_with_paths.append((
                validate_federation_row(row, f"PLACO {pair}/{phase} row"),
                receipt, provenance_path, int(provenance["peak_process_tree_rss_bytes"]),
            ))
        orphan = [
            path for path in entries if path.name.endswith(".tsv")
            and not path.with_suffix(".provenance.json").is_file()
        ]
        if orphan:
            raise FederationError(f"PLACO {pair} RAM receipt lacks provenance")
    if fingerprints_by_pair != {
        pair: contexts[pair][0]["fingerprint"] for pair in EXPECTED_PLACO_PAIR_ORDER
    }:
        raise FederationError("PLACO measurement namespaces differ from terminal pair contexts")
    order = {pair: index for index, pair in enumerate(EXPECTED_PLACO_PAIR_ORDER)}
    rows_with_paths.sort(key=lambda item: (order[item[0]["pair"]], item[1].name))
    rows = [item[0] for item in rows_with_paths]
    receipt_paths = [item[1] for item in rows_with_paths]
    provenance_paths = [item[2] for item in rows_with_paths]
    peak_rss_bytes = [item[3] for item in rows_with_paths]
    aggregate_content = module.tsv_bytes(module.RAM_FIELDS, rows)
    aggregate_path = root / module.RAM_AGGREGATE
    if stable_bytes(aggregate_path, "PLACO aggregate RAM table") != aggregate_content:
        raise FederationError("PLACO aggregate RAM table is stale")
    family_fingerprint = hashlib.sha256(module.canonical_json_bytes({
        "pair_order": list(EXPECTED_PLACO_PAIR_ORDER),
        "pair_fingerprints": fingerprints_by_pair,
    })).hexdigest()
    expected_provenance = {
        "schema_version": "sleep-atlas-track-b-placo-ram-benchmark-derived.1",
        "analysis_id": "track-b-v1.0-pleiotropy",
        "row_count": len(rows), "pair_order": list(EXPECTED_PLACO_PAIR_ORDER),
        "execution_fingerprints_by_pair": fingerprints_by_pair,
        "execution_family_fingerprint": family_fingerprint,
        "sequential_contract": module.stable_identity(root / module.EXECUTION_CONTRACT),
        "source_receipts": [module.stable_identity(path) for path in receipt_paths],
        "source_measurement_provenance": [
            module.stable_identity(path) for path in provenance_paths
        ],
        "aggregate": {
            "path": str(module.RAM_AGGREGATE), "bytes": len(aggregate_content),
            "sha256": hashlib.sha256(aggregate_content).hexdigest(),
            "schema": module.RAM_FIELDS,
        },
        "note": "Additive PLACO measurements; federation into the source-fingerprint-bound top-level table is separate.",
    }
    observed_provenance = module.read_json(root / module.RAM_AGGREGATE_PROVENANCE)
    if observed_provenance != expected_provenance:
        raise FederationError("PLACO aggregate RAM provenance is stale")
    return {
        "rows": rows,
        "peak_rss_bytes": peak_rss_bytes,
        "execution_fingerprint": family_fingerprint,
        "evidence": {
            "execution_family_fingerprint": family_fingerprint,
            "execution_contract": portable_identity(
                root, root / module.EXECUTION_CONTRACT, "PLACO execution contract",
            ),
            "aggregate": portable_identity(root, aggregate_path, "PLACO aggregate RAM table"),
            "aggregate_provenance": portable_identity(
                root, root / module.RAM_AGGREGATE_PROVENANCE, "PLACO aggregate RAM provenance",
            ),
            "cleanup_bundles": cleanup_bundle_identities,
            "cleanup_completions": cleanup_identities,
            "measurement_count": len(rows),
            "runtime_projections": projections,
            "global_nuisance_decomposable": False,
            "complete_pair_bh_decomposable": False,
        },
    }


def validate_collation_component(root: Path) -> dict[str, Any]:
    module = _load_module("track_b_ram_federation_collation_147", root / COLLATION_REL)
    try:
        verified = module.verify_results(root)
    except Exception as error:
        raise FederationError(f"script 147 deep result/component verification failed: {error}") from error
    complete = verified.get("complete_family")
    if (
        not isinstance(complete, Mapping)
        or complete.get("order") != list(EXPECTED_PLACO_PAIR_ORDER)
        or not isinstance(complete.get("pairs"), Mapping)
        or set(complete["pairs"]) != set(EXPECTED_PLACO_PAIR_ORDER)
    ):
        raise FederationError("script 147 verifier did not return exact A/B/CONTROL collation")
    rows, _ = read_tsv(root / module.RAM_COMPONENT_PATH, BENCHMARK_FIELDS, "script 147 RAM component")
    if len(rows) != 1:
        raise FederationError("script 147 RAM component must have exactly one row")
    rows = [validate_federation_row(rows[0], "script 147 RAM component row")]
    result_provenance, _ = read_json(
        root / module.RESULT_PROVENANCE_PATH, "script 147 result provenance",
    )
    try:
        peak_rss_bytes = int(result_provenance["resource_metrics"]["peak_ram_bytes"])
    except (KeyError, TypeError, ValueError) as error:
        raise FederationError("script 147 result provenance lacks exact peak RAM bytes") from error
    return {
        "rows": rows,
        "peak_rss_bytes": [peak_rss_bytes],
        "execution_fingerprint": str(verified["fingerprint"]),
        "evidence": {
            "deep_verifier_state": verified["state"],
            "collation_fingerprint": verified["fingerprint"],
            "component": portable_identity(
                root, root / module.RAM_COMPONENT_PATH, "script 147 RAM component",
            ),
            "namespace": portable_identity(
                root, root / module.RAM_NAMESPACE_PATH, "script 147 RAM namespace",
            ),
            "component_report": portable_identity(
                root, root / module.RAM_REPORT_PATH, "script 147 RAM report",
            ),
            "result_provenance": portable_identity(
                root, root / module.RESULT_PROVENANCE_PATH, "script 147 result provenance",
            ),
            "method_availability": verified["method_availability"],
        },
    }


def require_finemap_terminal_state(state: Mapping[str, Any], zero_state: str) -> str:
    observed = state.get("state")
    if observed not in {"COMPLETE_CANONICAL_FAMILY", zero_state}:
        raise FederationError(f"fine-mapping family is incomplete: {observed}")
    return str(observed)


def validate_finemap_component(root: Path) -> dict[str, Any]:
    verifier_identity = portable_identity(
        root, root / FINEMAP_REL, "fine-mapping public verifier",
    )
    if verifier_identity["sha256"] != FINEMAP_VERIFIER_SHA256:
        raise FederationError("fine-mapping public verifier is not the reconciled terminal version")
    module = _load_module("track_b_ram_federation_finemap_151", root / FINEMAP_REL)
    try:
        state = module.verify_current_state(root)
    except Exception as error:
        raise FederationError(f"fine-mapping deepest-state verification failed: {error}") from error
    accepted = require_finemap_terminal_state(state, FINEMAP_ZERO_STATE)
    if accepted == FINEMAP_ZERO_STATE:
        return {
            "rows": [], "execution_fingerprint": FINEMAP_ZERO_STATE,
            "peak_rss_bytes": [],
            "evidence": {
                "state": FINEMAP_ZERO_STATE, "component_row_count": 0,
                "verifier": verifier_identity, "genuine_larger_host_stages": [],
            },
        }
    canonical = state.get("canonical")
    if not isinstance(canonical, Mapping) or canonical.get("state") != "COMPLETE_CANONICAL_FAMILY":
        raise FederationError("fine-mapping canonical verifier result is malformed")
    component_hint = canonical.get("component_report")
    if not isinstance(component_hint, Mapping):
        raise FederationError("fine-mapping canonical family lacks component RAM verification")
    report_identity = component_hint.get("report_identity")
    if not isinstance(report_identity, Mapping) or not isinstance(report_identity.get("path"), str):
        raise FederationError("fine-mapping canonical component report identity is absent")
    report_rel = Path(str(report_identity["path"]))
    try:
        component = module.verify_component_ram_report(root, report_rel)
    except Exception as error:
        raise FederationError(f"fine-mapping component RAM verifier failed: {error}") from error
    if component.get("state") != "VERIFIED_COMPONENT_RAM_REPORT":
        raise FederationError("fine-mapping component verifier did not return its terminal state")
    if component.get("collision_key") != COLLISION_KEY:
        raise FederationError("fine-mapping component verifier returned a different collision key")
    for key in ("report_identity", "namespace_identity", "provenance_identity"):
        if component.get(key) != component_hint.get(key):
            raise FederationError("fine-mapping component/canonical identity binding drifted")
    rows = [
        validate_federation_row(row, f"fine-mapping component row {index}")
        for index, row in enumerate(component["rows"], start=1)
    ]
    blocked = [row for row in rows if row["exit_status"].startswith("BLOCKED_BY_DATA:")]
    other_incomplete = [
        row for row in rows
        if row["exit_status"] != "COMPLETE_SEALED"
        and not row["exit_status"].startswith("BLOCKED_BY_DATA:")
    ]
    if other_incomplete:
        raise FederationError("fine-mapping canonical RAM component retains nonterminal mapped rows")
    if any(
        row["n_snps"] != "0" or row["peak_ram_gb"] != "0" or row["runtime_sec"] != "0"
        for row in blocked
    ):
        raise FederationError("fine-mapping unavailable locus row is not an explicit data block")
    canonical_runs = canonical.get("runs")
    if not isinstance(canonical_runs, list):
        raise FederationError("fine-mapping canonical verifier omitted its exact run family")
    runs_by_locus: dict[str, Mapping[str, Any]] = {}
    for run in canonical_runs:
        if not isinstance(run, Mapping) or not isinstance(run.get("row"), Mapping):
            raise FederationError("fine-mapping canonical verifier returned a malformed run")
        locus = str(run["row"].get("locus_entry_id", ""))
        if not locus or locus in runs_by_locus:
            raise FederationError("fine-mapping canonical verifier repeats a run locus")
        runs_by_locus[locus] = run
    peak_rss_bytes: list[int] = []
    larger_stages: list[dict[str, object]] = []
    for row in rows:
        if row["exit_status"].startswith("BLOCKED_BY_DATA:"):
            peak_rss_bytes.append(0)
            continue
        run = runs_by_locus.pop(row["locus"], None)
        if run is None:
            raise FederationError("fine-mapping component row lacks its verified canonical run")
        run_row = run["row"]
        metric = run.get("resource_row")
        provenance = run.get("provenance")
        if (
            not isinstance(metric, Mapping) or not isinstance(provenance, Mapping)
            or str(run_row.get("pair_id")) != row["pair"]
            or str(run_row.get("CHR")) != row["chromosome"]
        ):
            raise FederationError("fine-mapping component/run identity binding drifted")
        try:
            exact_peak = int(metric["peak_rss_bytes"])
        except (KeyError, TypeError, ValueError) as error:
            raise FederationError("fine-mapping run lacks exact peak RSS bytes") from error
        if format(exact_peak / GIB, ".12g") != row["peak_ram_gb"]:
            raise FederationError("fine-mapping component peak differs from exact run receipt")
        peak_rss_bytes.append(exact_peak)
        materialization = run.get("materialization")
        materialization_receipt = (
            materialization.get("receipt") if isinstance(materialization, Mapping) else None
        )
        if not isinstance(materialization_receipt, Mapping):
            raise FederationError("fine-mapping run lacks its verified materialization receipt")
        try:
            materialization_peak = int(materialization_receipt["peak_rss_bytes"])
        except (KeyError, TypeError, ValueError) as error:
            raise FederationError("fine-mapping materialization lacks exact peak RSS bytes") from error
        if materialization_peak < 0 or materialization_peak > exact_peak:
            raise FederationError("fine-mapping materialization peak contradicts the component maximum")
        if materialization_receipt.get("larger_host_continuation") is True:
            larger_stages.append({
                "analysis": row["analysis"], "pair": row["pair"],
                "locus": row["locus"], "chromosome": row["chromosome"],
                "stage": "MATERIALIZATION", "peak_rss_bytes": materialization_peak,
                "predecessor_count": None,
                "lineage_verifier": "SCRIPT_151_VERIFY_MATERIALIZATION_MEASUREMENT",
            })
        if provenance.get("larger_host_continuation") is True:
            predecessors = provenance.get("engine_continuation_predecessors")
            if not isinstance(predecessors, list) or not predecessors:
                raise FederationError("fine-mapping larger-host run lacks predecessor lineage")
            larger_stages.append({
                "analysis": row["analysis"], "pair": row["pair"],
                "locus": row["locus"], "chromosome": row["chromosome"],
                "stage": "R_ENGINE", "peak_rss_bytes": exact_peak,
                "predecessor_count": len(predecessors),
                "lineage_verifier": "SCRIPT_151_VERIFY_SEALED_RUN",
            })
    if runs_by_locus:
        raise FederationError("fine-mapping canonical run family exceeds its component rows")
    canonical_provenance = canonical.get("provenance_identity")
    if (
        not isinstance(canonical_provenance, Mapping)
        or not SHA256_RE.fullmatch(str(canonical_provenance.get("sha256", "")))
    ):
        raise FederationError("fine-mapping canonical result provenance identity is absent")
    return {
        "rows": rows,
        "peak_rss_bytes": peak_rss_bytes,
        "execution_fingerprint": str(canonical_provenance["sha256"]),
        "evidence": {
            "state": accepted, "component_row_count": len(rows),
            "blocked_by_data_rows": len(blocked),
            "report": component["report_identity"],
            "namespace": component["namespace_identity"],
            "provenance": component["provenance_identity"],
            "canonical_result_provenance": canonical_provenance,
            "verifier": verifier_identity,
            "genuine_larger_host_stages": larger_stages,
            "component_verified_only_via_script_151_public_api": True,
        },
    }


def _is_success(row: Mapping[str, str]) -> bool:
    return row["exit_status"] in {"0", "COMPLETE_SEALED"}


def _component_peak_rss_bytes(component: Mapping[str, Any]) -> list[int | None]:
    values = component.get("peak_rss_bytes")
    rows = component.get("rows")
    if (
        not isinstance(values, list) or not isinstance(rows, Sequence)
        or isinstance(rows, (str, bytes)) or len(values) != len(rows)
    ):
        raise FederationError(f"{component.get('name', 'component')} lacks one peak-byte value per row")
    normalized: list[int | None] = []
    for value in values:
        if value is None:
            normalized.append(None)
        elif isinstance(value, int) and not isinstance(value, bool) and value >= 0:
            normalized.append(value)
        else:
            raise FederationError("component peak-byte evidence is malformed")
    return normalized


def _merge_peak_rss_bytes(components: Sequence[Mapping[str, Any]]) -> list[int | None]:
    return [value for component in components for value in _component_peak_rss_bytes(component)]


def _analysis_statistics(
    rows: Sequence[Mapping[str, str]], peak_rss_bytes: Sequence[int | None],
) -> list[dict[str, object]]:
    if len(rows) != len(peak_rss_bytes):
        raise FederationError("federation row/peak-byte families differ in length")
    grouped: dict[str, list[tuple[Mapping[str, str], int | None]]] = collections.defaultdict(list)
    for row, peak_bytes in zip(rows, peak_rss_bytes, strict=True):
        grouped[row["analysis"]].append((row, peak_bytes))
    output: list[dict[str, object]] = []
    for analysis in sorted(grouped):
        family = grouped[analysis]
        family_rows = [item[0] for item in family]
        exact_peaks = [item[1] for item in family]
        runtimes = [float(row["runtime_sec"]) for row in family_rows]
        peaks = [float(row["peak_ram_gb"]) for row in family_rows]
        exact_peak = max((value for value in exact_peaks if value is not None), default=None)
        all_exact = all(value is not None for value in exact_peaks)
        fits = (
            exact_peak + RESERVE_BYTES <= EIGHT_GIB
            if all_exact and exact_peak is not None else None
        )
        output.append({
            "analysis": analysis, "units": len(family),
            "successful": sum(_is_success(row) for row in family_rows),
            "retained_non_success": sum(not _is_success(row) for row in family_rows),
            "peak_gib": max(peaks, default=0.0),
            "peak_rss_bytes": exact_peak,
            "observed_seconds": sum(runtimes),
            "fits": fits,
        })
    return output


def build_report(
    fingerprint: str, rows: Sequence[Mapping[str, str]], components: Sequence[Mapping[str, Any]],
    peak_rss_bytes: Sequence[int | None] | None = None,
) -> bytes:
    exact_peaks = list(peak_rss_bytes) if peak_rss_bytes is not None else _merge_peak_rss_bytes(components)
    statistics_rows = _analysis_statistics(rows, exact_peaks)
    failures = [row for row in rows if not _is_success(row)]
    components_by_name = {str(item["name"]): item for item in components}
    placo_projection = components_by_name.get("placo_143", {}).get("evidence", {}).get(
        "runtime_projections", {}
    )
    larger_measured = [item for item in statistics_rows if item["fits"] is False]
    unknown_fit = [item for item in statistics_rows if item["fits"] is None]
    larger_lineages = [
        stage
        for component in components
        for stage in component.get("evidence", {}).get("genuine_larger_host_stages", [])
    ]
    lines = [
        "# Final Track B RAM-aware execution report",
        "",
        f"Federation fingerprint: `{fingerprint}`.",
        "",
        "This final table is measurement/provenance only. It changes no scientific model, threshold, family denominator, phenotype pairing, or claim rule. The provisional V1 provenance was intentionally not trusted: all 2,495 discovery rows were reconstructed from the exact source-lock READY receipts and result identities, and the historical failed aggregate log was rehashed.",
        "",
        "## Verified component boundaries",
        "",
        f"- Frozen LAVA V1 discovery: `{V1_FINGERPRINT}` (2,495 READY results plus its retained failed aggregate attempt).",
        f"- Production LAVA continuation: `{CONTINUATION_FINGERPRINT}` (the 09cd and ceb continuations are verified superseded and excluded).",
        "- PLACO+: all A, B, and CONTROL canonical ledgers, cleanup plans/completion receipts, and every immutable RAM measurement provenance validate under script 143.",
        "- Post-PLACO collation: script 147 deep result verification completed before its component row was admitted.",
        "- Fine-mapping/trait-coloc: script 151 canonical/terminal verification and its versioned component RAM verifier completed before rows were admitted; unavailable loci remain explicit `BLOCKED_BY_DATA` rows.",
        "",
        "## Safe decomposition",
        "",
        "LAVA discovery and conditional work are decomposed by official locus and executed sequentially; whole-family aggregation/finalization remains a separate sealed boundary. PLACO+ is decomposed only by the frozen A, B, and CONTROL pair families and those pairs run sequentially. Post-PLACO ledgers are streamed and PLINK commands run one pair at a time. Fine-mapping and trait-coloc run one frozen locus in a fresh process at a time; a locus itself is not subdivided.",
        "",
        "## Non-decomposable PLACO stages",
        "",
        "The PLACO global nuisance variance/correlation fit uses all valid variants for one pair and is not decomposable. The complete within-pair PLACO p-value family and its BH denominator are also indivisible. The table retains both the freshly matched nuisance benchmark and the measured complete scan; neither is represented by a subset benchmark.",
        "",
        "## Measured peak RAM and sequential 8-GiB feasibility",
        "",
        "The admission rule is sequential execution on an 8-GiB host with a 1-GiB reserve, so a measured stage passes only when its peak is at most 7 GiB.",
        "",
        "| Analysis | Units | Successful | Retained failures/skips | Measured peak GiB | Exact peak bytes | Observed cumulative runtime s | Exact peak + 1 GiB <= 8 GiB |",
        "|---|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for item in statistics_rows:
        lines.append(
            f"| {item['analysis']} | {item['units']} | {item['successful']} | "
            f"{item['retained_non_success']} | {item['peak_gib']:.6f} | "
            f"{item['peak_rss_bytes'] if item['peak_rss_bytes'] is not None else 'UNKNOWN'} | "
            f"{item['observed_seconds']:.3f} | "
            f"{ {True: 'YES', False: 'NO', None: 'UNKNOWN'}[item['fits']] } |"
        )
    lines.extend([
        "",
        "The table reports observed cumulative runtime only; it does not relabel observations as an independent expectation. The independently sealed fresh PLACO matched-worker projections are retained below for comparison:",
        "",
        "| PLACO pair | Fresh benchmark measured s | Projected complete scan s | Global nuisance s |",
        "|---|---:|---:|---:|",
    ])
    if isinstance(placo_projection, Mapping) and placo_projection:
        for pair in EXPECTED_PLACO_PAIR_ORDER:
            item = placo_projection.get(pair, {})
            lines.append(
                f"| {pair} | {item.get('measured_runner_elapsed_seconds', 'NA')} | "
                f"{item.get('projected_full_family_seconds', 'NA')} | "
                f"{item.get('global_nuisance_elapsed_seconds', 'NA')} |"
            )
    else:
        lines.append("| No PLACO projection | NA | NA | NA |")
    lines.extend([
        "",
        "## Genuine larger-server stages",
        "",
    ])
    if larger_measured:
        for item in larger_measured:
            lines.append(
                f"- `{item['analysis']}` measured {item['peak_gib']:.6f} GiB and does not fit the 8-GiB/1-GiB-reserve envelope."
            )
    elif not larger_lineages:
        lines.append("- No completed measured stage in this federation requires a larger server under the frozen sequential rule.")
    for item in unknown_fit:
        lines.append(
            f"- `{item['analysis']}` has retained display-only historical RAM evidence; its exact-byte 8-GiB fit is `UNKNOWN`."
        )
    for stage in larger_lineages:
        predecessor_text = (
            f"{stage['predecessor_count']} sealed resource predecessor(s)"
            if stage.get("predecessor_count") is not None
            else "a nonempty failed-resource predecessor family verified inside script 151"
        )
        lines.append(
            f"- `{stage['analysis']}` pair `{stage['pair']}` locus `{stage['locus']}` "
            f"stage `{stage.get('stage', 'UNSPECIFIED')}` is a verified genuine larger-host "
            f"continuation with {predecessor_text}, independent of rounded display RAM."
        )
    lines.append(
        "- The independently blocked MATLAB conjFDR method remains a genuine >=16-GiB RAM / >=20-GiB free-storage stage under its sealed policy. It is not silently substituted, inferred, or counted as completed PLACO work."
    )
    lines.extend([
        "",
        "## Full retained failures and skips",
        "",
        "Every non-success row remains in `RAM_BENCHMARK.tsv`; none was filtered before these counts were computed.",
        "",
        "| Analysis | Pair | Locus | Chromosome | Exit status | Output hash |",
        "|---|---|---|---|---|---|",
    ])
    if failures:
        for row in failures:
            lines.append(
                f"| {row['analysis']} | {row['pair']} | {row['locus']} | "
                f"{row['chromosome']} | {row['exit_status']} | `{row['output_hash']}` |"
            )
    else:
        lines.append("| None | NONE | NA | NA | NONE | NA |")
    lines.extend([
        "",
        "## Publication and recovery",
        "",
        "The exact nine-column table, namespace, report, and provenance are sealed in one versioned package. For each top-level role, the frozen V1 inode is first atomically withdrawn to a retained transaction path without replacement, then the mode-0444 packaged inode is linked into the now-empty destination without replacement. Promotion is resumable in the fixed order table, namespace, report, provenance; provenance is last and is the commit marker. A partial prefix is accepted only when every promoted byte equals this package and every remaining or withdrawn source equals its frozen V1 identity.",
        "",
        "Preflight deliberately retains script 138's independent terminal-family reconstruction. Shared c139 identities are cached in the outer pass, but the frozen terminal safeguard may reread its complete active-input receipts; this verification cost is not presented as analysis runtime.",
        "",
    ])
    return ("\n".join(lines)).encode("utf-8")


def _component_binding(component: Mapping[str, Any]) -> dict[str, Any]:
    peak_rss_bytes = _component_peak_rss_bytes(component)
    return {
        "name": component["name"],
        "execution_fingerprint": component.get("execution_fingerprint"),
        "row_count": len(component["rows"]),
        "rows_sha256": hashlib.sha256(tsv_bytes(BENCHMARK_FIELDS, component["rows"])).hexdigest(),
        "exact_peak_rss_bytes": peak_rss_bytes,
        "evidence": component["evidence"],
    }


def build_candidate(
    root: Path, components: Sequence[Mapping[str, Any]], history_receipt: Mapping[str, Any],
) -> dict[str, Any]:
    rows = merge_component_rows(components)
    peak_rss_bytes = _merge_peak_rss_bytes(components)
    benchmark_content = tsv_bytes(BENCHMARK_FIELDS, rows)
    script_identity = portable_identity(root, root / SCRIPT_REL, "RAM federation script")
    bindings = [_component_binding(component) for component in components]
    fingerprint = digest_json({
        "schema_version": SCHEMA,
        "source_fingerprint": V1_FINGERPRINT,
        "continuation_fingerprint": CONTINUATION_FINGERPRINT,
        "superseded_excluded": list(SUPERSEDED_CONTINUATIONS),
        "benchmark_sha256": hashlib.sha256(benchmark_content).hexdigest(),
        "component_bindings": bindings,
        "script": script_identity,
    })
    report_content = build_report(fingerprint, rows, components, peak_rss_bytes)
    namespace = {
        "schema_version": NAMESPACE_SCHEMA,
        "analysis_id": "track-b-v1.0-final-ram-federation",
        "federation_fingerprint": fingerprint,
        "table": str(TOP_LEVEL["benchmark"]),
        "exact_fields": BENCHMARK_FIELDS,
        "collision_key": COLLISION_KEY,
        "row_count": len(rows),
        "source_discovery_fingerprint": V1_FINGERPRINT,
        "production_continuation_fingerprint": CONTINUATION_FINGERPRINT,
        "superseded_continuations_excluded": list(SUPERSEDED_CONTINUATIONS),
        "component_row_counts": {
            str(component["name"]): len(component["rows"]) for component in components
        },
    }
    namespace_content = json.dumps(namespace, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    package_rel = PACKAGE_ROOT_REL / fingerprint
    old_sources = history_receipt.get("source_identities_at_freeze")
    if not isinstance(old_sources, Mapping) or set(old_sources) != set(TOP_LEVEL):
        raise FederationError("V1 history receipt lacks the exact four old top-level identities")
    expected_old: dict[str, dict[str, object]] = {}
    for role in TOP_LEVEL:
        record = old_sources[role]
        if (
            not isinstance(record, Mapping)
            or record.get("path") != str(TOP_LEVEL[role])
            or not isinstance(record.get("bytes"), int)
            or not SHA256_RE.fullmatch(str(record.get("sha256", "")))
        ):
            raise FederationError(f"V1 old top-level source identity is malformed: {role}")
        expected_old[role] = dict(record)
    output_identities = {
        "benchmark": identity_for_bytes(TOP_LEVEL["benchmark"], benchmark_content),
        "namespace": identity_for_bytes(TOP_LEVEL["namespace"], namespace_content),
        "report": identity_for_bytes(TOP_LEVEL["report"], report_content),
    }
    provenance = {
        "schema_version": SCHEMA,
        "analysis_id": "track-b-v1.0-final-ram-federation",
        "artifact_role": "FINAL_TERMINAL_COMPONENT_RAM_FEDERATION",
        "federation_fingerprint": fingerprint,
        "exact_fields": BENCHMARK_FIELDS,
        "collision_key": COLLISION_KEY,
        "collision_key_unique": True,
        "row_count": len(rows),
        "terminal_status_counts": dict(sorted(collections.Counter(
            row["exit_status"] for row in rows
        ).items())),
        "outputs_before_commit_marker": output_identities,
        "component_bindings": bindings,
        "historical_v1": {
            "fingerprint": V1_FINGERPRINT,
            "source_identities_at_freeze": expected_old,
            "stale_provenance_trusted": False,
            "raw_source_ready_receipts_revalidated": EXPECTED_DISCOVERY_COUNT,
        },
        "production_continuation_fingerprint": CONTINUATION_FINGERPRINT,
        "superseded_continuations_excluded": list(SUPERSEDED_CONTINUATIONS),
        "memory_envelope": {
            "sequential_host_bytes": EIGHT_GIB,
            "reserve_bytes": RESERVE_BYTES,
            "maximum_admitted_measured_peak_bytes": EIGHT_GIB - RESERVE_BYTES,
            "fit_decision_source": "EXACT_PEAK_RSS_BYTES_OR_UNKNOWN",
        },
        "scientific_contract_changed": False,
        "package": {
            "path": str(package_rel), "manifest": str(package_rel / PACKAGE_MANIFEST),
            "promotion_order": [str(TOP_LEVEL[role]) for role in PROMOTION_ORDER],
            "commit_marker": str(TOP_LEVEL["provenance"]),
            "promotion": "PRESERVATION_FIRST_WITHDRAW_THEN_KERNEL_NO_REPLACE_LINK",
            "withdrawal_recovery_root": str(TRANSACTION_ROOT_REL / fingerprint),
            "withdrawn_sources_retained": True,
            "final_artifact_mode": "0444",
        },
        "generator": script_identity,
    }
    provenance_content = json.dumps(provenance, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    contents = {
        "benchmark": benchmark_content, "namespace": namespace_content,
        "report": report_content, "provenance": provenance_content,
    }
    return {
        "fingerprint": fingerprint, "rows": rows, "peak_rss_bytes": peak_rss_bytes,
        "contents": contents,
        "old_sources": expected_old, "package_rel": package_rel,
        "namespace": namespace, "provenance": provenance,
    }


def collect_federation(root: Path = ROOT) -> dict[str, Any]:
    """Verify inputs and build the deterministic candidate without repository writes.

    The history freezer is deliberately the first imported component verifier,
    and its immutable snapshot is verified before any live family is inspected.
    """

    history_module = _load_module("track_b_ram_federation_history_152", root / HISTORY_SCRIPT_REL)
    try:
        history_receipt = history_module.verify_snapshot(root, V1_FINGERPRINT)
    except Exception as error:
        raise FederationError(f"immutable V1 RAM history failed verification: {error}") from error
    if (
        history_receipt.get("benchmark_row_count") != EXPECTED_HISTORICAL_ROWS
        or history_receipt.get("historical_provenance_audit", {}).get("declared_attempt_rows")
        != EXPECTED_STALE_DECLARED_ROWS
    ):
        raise FederationError("immutable V1 stale-provenance accounting is not the known 2496/1376 state")

    contract = _load_module("track_b_ram_federation_source_contract_134", root / CONTINUATION_CONTRACT_REL)
    try:
        source_lock = contract.validate_source_lock(
            # The contract validates the immutable lock itself.  Its legacy
            # ``validate_bundles=True`` path rehashes the entire multi-gigabyte
            # active-input baseline separately for each of 2,495 loci.  The
            # federation instead performs the stricter, bounded loop below:
            # exact canonical READY link plus every lock-recorded READY,
            # receipt, result, worker-log, and semantic-attestation identity.
            V1_FINGERPRINT, validate_bundles=False,
            revalidate_semantics=False, require_live_fingerprint=False,
        )
    except Exception as error:
        raise FederationError(f"V1 source lock/raw READY family failed verification: {error}") from error
    benchmark_record = history_receipt["snapshots"]["benchmark"]
    benchmark_path = safe_path(root, benchmark_record["path"], "V1 history benchmark")
    v1 = validate_v1_snapshot_rows(root, benchmark_path, history_receipt, source_lock)
    superseded = validate_superseded_exclusions(root)
    continuation = validate_continuation_component(root, source_lock, v1["loci"])
    placo = validate_placo_component(root)
    collation = validate_collation_component(root)
    finemap = validate_finemap_component(root)
    components = [
        {
            "name": "lava_v1_history", "rows": v1["rows"],
            "peak_rss_bytes": v1["peak_rss_bytes"],
            "execution_fingerprint": V1_FINGERPRINT,
            "evidence": {
                "history_receipt": portable_identity(
                    root,
                    safe_path(
                        root,
                        Path(history_module.HISTORY_ROOT) / V1_FINGERPRINT / history_module.RECEIPT_NAME,
                        "V1 history receipt",
                    ),
                    "V1 history receipt",
                ),
                "source_lock": portable_identity(
                    root, contract.source_lock_path(V1_FINGERPRINT), "V1 source lock",
                ),
                "source_checkpoint_family_sha256": source_lock["checkpoint_family_sha256"],
                "raw_revalidation_family_sha256": v1["checkpoint_family_revalidation_sha256"],
                "stale_declared_rows": v1["stale_declared_rows"],
                "actual_rows": v1["actual_rows"],
                "failed_aggregate_log": v1["failed_aggregate_log"],
            },
        },
        {
            "name": "lava_superseded_exclusions", "rows": [],
            "peak_rss_bytes": [],
            "execution_fingerprint": digest_json(superseded),
            "evidence": {
                "excluded_fingerprints": list(SUPERSEDED_CONTINUATIONS),
                "verified_final_locks": superseded,
            },
        },
        {"name": "lava_c139_continuation", **continuation},
        {"name": "placo_143", **placo},
        {"name": "post_placo_147", **collation},
        {"name": "finemapping_151", **finemap},
    ]
    candidate = build_candidate(root, components, history_receipt)
    return candidate


def _package_manifest(candidate: Mapping[str, Any]) -> tuple[dict[str, Any], bytes]:
    artifacts = {
        PACKAGE_FILES[role]: {
            "bytes": len(candidate["contents"][role]),
            "sha256": hashlib.sha256(candidate["contents"][role]).hexdigest(),
            "top_level_path": str(TOP_LEVEL[role]),
        }
        for role in PROMOTION_ORDER
    }
    payload = {
        "schema_version": PACKAGE_SCHEMA,
        "federation_fingerprint": candidate["fingerprint"],
        "artifacts": artifacts,
        "promotion_order": [PACKAGE_FILES[role] for role in PROMOTION_ORDER],
        "commit_marker": PACKAGE_FILES["provenance"],
        "top_level_commit_marker": str(TOP_LEVEL["provenance"]),
        "promotion": "PRESERVATION_FIRST_WITHDRAW_THEN_KERNEL_NO_REPLACE_LINK",
        "withdrawn_sources_recoverable": True,
        "partial_prefix_recovery": True,
    }
    return payload, json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_new_file(path: Path, content: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    try:
        os.fchmod(descriptor, 0o444)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        try:
            path.unlink()
        except FileNotFoundError:
            pass
        raise


def verify_package(root: Path, candidate: Mapping[str, Any]) -> dict[str, Any]:
    package = safe_path(root, candidate["package_rel"], "RAM federation package")
    if package.is_symlink() or not package.is_dir():
        raise FederationError("versioned RAM federation package is absent or unsafe")
    expected_manifest, manifest_content = _package_manifest(candidate)
    names = {item.name for item in package.iterdir()}
    if names != set(PACKAGE_FILES.values()) | {PACKAGE_MANIFEST}:
        raise FederationError("versioned RAM federation package inventory drifted")
    observed_manifest, observed_content = read_json(
        package / PACKAGE_MANIFEST, "RAM federation package manifest",
    )
    if observed_manifest != expected_manifest or observed_content != manifest_content:
        raise FederationError("versioned RAM federation package manifest drifted")
    if stat.S_IMODE(os.lstat(package / PACKAGE_MANIFEST).st_mode) != 0o444:
        raise FederationError("versioned RAM federation package manifest is not mode 0444")
    identities: dict[str, dict[str, object]] = {}
    for role in PROMOTION_ORDER:
        path = package / PACKAGE_FILES[role]
        content = stable_bytes(path, f"packaged {role}")
        if content != candidate["contents"][role]:
            raise FederationError(f"versioned RAM federation package {role} drifted")
        if stat.S_IMODE(os.lstat(path).st_mode) != 0o444:
            raise FederationError(f"versioned RAM federation package {role} is not mode 0444")
        identities[role] = portable_identity(root, path, f"packaged {role}")
    return {
        "path": str(candidate["package_rel"]),
        "manifest": identity_for_bytes(
            candidate["package_rel"] / PACKAGE_MANIFEST, manifest_content,
        ),
        "artifacts": identities,
    }


def _create_package(root: Path, candidate: Mapping[str, Any]) -> dict[str, Any]:
    package = safe_path(root, candidate["package_rel"], "RAM federation package destination")
    if package.exists() or package.is_symlink():
        return verify_package(root, candidate)
    package_root = safe_path(root, PACKAGE_ROOT_REL, "RAM federation package root")
    package_root.mkdir(parents=True, exist_ok=True)
    if package_root.is_symlink() or not package_root.is_dir():
        raise FederationError("RAM federation package root is unsafe")
    stage = Path(tempfile.mkdtemp(prefix=f".{candidate['fingerprint']}.", suffix=".staging", dir=package_root))
    try:
        for role in PROMOTION_ORDER:
            _write_new_file(stage / PACKAGE_FILES[role], candidate["contents"][role])
        _, manifest_content = _package_manifest(candidate)
        _write_new_file(stage / PACKAGE_MANIFEST, manifest_content)
        _fsync_directory(stage)
        try:
            os.rename(stage, package)
        except OSError as error:
            if error.errno not in {errno.EEXIST, errno.ENOTEMPTY}:
                raise
        _fsync_directory(package_root)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return verify_package(root, candidate)


def _transaction_backup(root: Path, candidate: Mapping[str, Any], role: str) -> Path:
    fingerprint = str(candidate["fingerprint"])
    if not SHA256_RE.fullmatch(fingerprint) or role not in PROMOTION_ORDER:
        raise FederationError("invalid RAM federation transaction identity")
    return safe_path(
        root, TRANSACTION_ROOT_REL / fingerprint / f"{role}.withdrawn",
        f"withdrawn top-level {role}",
    )


def _matches_content_identity(content: bytes, identity: Mapping[str, Any]) -> bool:
    return (
        len(content) == identity.get("bytes")
        and hashlib.sha256(content).hexdigest() == identity.get("sha256")
    )


def publication_state(root: Path, candidate: Mapping[str, Any]) -> dict[str, Any]:
    states: list[str] = []
    for role in PROMOTION_ORDER:
        target = safe_path(root, TOP_LEVEL[role], f"top-level {role}")
        backup = _transaction_backup(root, candidate, role)
        if target.is_symlink() or backup.is_symlink():
            raise FederationError(f"top-level RAM {role} transaction contains a symlink")
        old = candidate["old_sources"][role]
        backup_present = backup.exists()
        if backup_present:
            withdrawn = stable_bytes(backup, f"withdrawn top-level RAM {role}")
            if not _matches_content_identity(withdrawn, old):
                raise FederationError(f"withdrawn top-level RAM {role} is not the frozen V1 source")
        if target.exists():
            content = stable_bytes(target, f"top-level RAM {role}")
            if content == candidate["contents"][role]:
                if stat.S_IMODE(os.lstat(target).st_mode) != 0o444:
                    raise FederationError(f"promoted top-level RAM {role} is not mode 0444")
                states.append("NEW")
            elif _matches_content_identity(content, old):
                if backup_present:
                    raise FederationError(f"top-level RAM {role} duplicates its withdrawn source")
                states.append("OLD")
            else:
                raise FederationError(f"top-level RAM {role} matches neither frozen V1 nor candidate")
        elif backup_present:
            states.append("WITHDRAWN")
        else:
            raise FederationError(f"top-level RAM {role} is missing without a recoverable withdrawal")
    prefix = 0
    while prefix < len(states) and states[prefix] == "NEW":
        prefix += 1
    withdrawn_role: str | None = None
    remainder = states[prefix:]
    if remainder and remainder[0] == "WITHDRAWN":
        withdrawn_role = PROMOTION_ORDER[prefix]
        remainder = remainder[1:]
    if any(state != "OLD" for state in remainder):
        raise FederationError("top-level RAM publication is not a valid recoverable prefix")
    if withdrawn_role is not None:
        state = f"WITHDRAWN_{withdrawn_role.upper()}_AFTER_PREFIX_{prefix}_OF_{len(states)}"
    elif prefix == 0:
        state = "UNPUBLISHED_FROZEN_V1"
    elif prefix == len(states):
        state = "COMMITTED"
    else:
        state = f"PARTIAL_PREFIX_{prefix}_OF_{len(states)}"
    return {
        "state": state, "promoted_prefix": prefix, "withdrawn_role": withdrawn_role,
        "members": dict(zip(PROMOTION_ORDER, states)),
    }


def _rename_no_replace(source: Path, destination: Path) -> None:
    """Atomically withdraw a path without ever replacing the recovery destination."""

    import ctypes
    import ctypes.util

    try:
        libc = ctypes.CDLL(ctypes.util.find_library("c") or None, use_errno=True)
    except OSError as error:
        raise FederationError(f"could not load libc for preservation-first promotion: {error}") from error
    source_raw, destination_raw = os.fsencode(source), os.fsencode(destination)
    if sys.platform == "darwin":
        try:
            rename = libc.renamex_np
        except AttributeError as error:
            raise FederationError("renamex_np is unavailable for preservation-first promotion") from error
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
        raise FederationError("kernel no-replace rename is unavailable")
    if result:
        number = ctypes.get_errno()
        if number in {errno.EEXIST, errno.ENOTEMPTY}:
            raise FederationError("RAM federation recovery destination appeared concurrently")
        raise FederationError(f"RAM federation no-replace rename failed: {os.strerror(number)}")


def _promote_file_preserving(
    root: Path, candidate: Mapping[str, Any], package: Path, role: str,
) -> None:
    target = safe_path(root, TOP_LEVEL[role], f"top-level {role}")
    backup = _transaction_backup(root, candidate, role)
    source = package / PACKAGE_FILES[role]
    backup.parent.mkdir(parents=True, exist_ok=True)
    if backup.parent.is_symlink() or not backup.parent.is_dir():
        raise FederationError("RAM federation transaction directory is unsafe")
    old = candidate["old_sources"][role]

    if backup.exists():
        withdrawn = stable_bytes(backup, f"withdrawn top-level RAM {role}")
        if not _matches_content_identity(withdrawn, old):
            if not target.exists() and not target.is_symlink():
                _rename_no_replace(backup, target)
                _fsync_directory(target.parent)
                if backup.parent != target.parent:
                    _fsync_directory(backup.parent)
            raise FederationError(f"withdrawn top-level RAM {role} is not recoverable frozen V1")
        if target.exists():
            if stable_bytes(target, f"top-level RAM {role}") == candidate["contents"][role]:
                return
            raise FederationError(f"top-level RAM {role} appeared beside its withdrawal")
    else:
        if not target.exists() or target.is_symlink():
            raise FederationError(f"top-level RAM {role} cannot be withdrawn safely")
        _rename_no_replace(target, backup)
        _fsync_directory(target.parent)
        if backup.parent != target.parent:
            _fsync_directory(backup.parent)
        withdrawn = stable_bytes(backup, f"withdrawn top-level RAM {role}")
        if not _matches_content_identity(withdrawn, old):
            if not target.exists() and not target.is_symlink():
                _rename_no_replace(backup, target)
                _fsync_directory(target.parent)
            raise FederationError(
                f"top-level RAM {role} changed concurrently; unexpected bytes were preserved"
            )

    try:
        os.link(source, target)
    except FileExistsError as error:
        raise FederationError(
            f"top-level RAM {role} appeared concurrently; no file was overwritten"
        ) from error
    _fsync_directory(target.parent)
    observed = stable_bytes(target, f"promoted top-level {role}")
    if observed != candidate["contents"][role]:
        raise FederationError(f"top-level RAM {role} promotion did not commit exact bytes")
    if stat.S_IMODE(os.lstat(target).st_mode) != 0o444:
        raise FederationError(f"promoted top-level RAM {role} is not mode 0444")


def _restore_interrupted_unexpected_withdrawals(
    root: Path, candidate: Mapping[str, Any],
) -> None:
    """Restore raced bytes if a process died between withdrawal and mismatch handling."""

    for role in PROMOTION_ORDER:
        target = safe_path(root, TOP_LEVEL[role], f"top-level {role}")
        backup = _transaction_backup(root, candidate, role)
        if backup.is_symlink() or target.is_symlink() or not backup.exists():
            continue
        withdrawn = stable_bytes(backup, f"withdrawn top-level RAM {role}")
        if _matches_content_identity(withdrawn, candidate["old_sources"][role]):
            continue
        if target.exists():
            raise FederationError(
                f"withdrawn top-level RAM {role} contains unexpected bytes beside a live target"
            )
        _rename_no_replace(backup, target)
        _fsync_directory(target.parent)
        if backup.parent != target.parent:
            _fsync_directory(backup.parent)
        raise FederationError(
            f"interrupted top-level RAM {role} race was restored without overwriting it"
        )


def publish_candidate(root: Path, candidate: Mapping[str, Any]) -> dict[str, Any]:
    # Use one stable lock inode that is never one of the four promoted paths.
    package_root = safe_path(root, PACKAGE_ROOT_REL, "RAM federation package root")
    package_root.mkdir(parents=True, exist_ok=True)
    lock_path = safe_path(root, PUBLICATION_LOCK_REL, "RAM federation publication lock")
    descriptor = os.open(lock_path, os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0), 0o600)
    with os.fdopen(descriptor, "rb+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        _restore_interrupted_unexpected_withdrawals(root, candidate)
        state = publication_state(root, candidate)
        package_evidence = _create_package(root, candidate)
        package = root / str(package_evidence["path"])
        for position, role in enumerate(PROMOTION_ORDER, start=1):
            if position <= state["promoted_prefix"]:
                continue
            # Revalidate the complete recoverable prefix immediately before
            # every mutation.  This also proves all unpromoted targets retain
            # the exact frozen V1 bytes.
            current = publication_state(root, candidate)
            if (
                current["promoted_prefix"] != position - 1
                or current["withdrawn_role"] not in {None, role}
            ):
                raise FederationError("top-level RAM publication changed concurrently")
            _promote_file_preserving(root, candidate, package, role)
        final = publication_state(root, candidate)
        if final["state"] != "COMMITTED":
            raise FederationError("top-level RAM promotion did not reach its commit marker")
    return {"package": package_evidence, "publication": final}


def preflight(root: Path = ROOT, *, collector: Callable[[Path], dict[str, Any]] = collect_federation) -> dict[str, Any]:
    with _suppress_bytecode_writes():
        candidate = collector(root)
        state = publication_state(root, candidate)
    return {
        "status": "PREFLIGHT_PASS", "federation_fingerprint": candidate["fingerprint"],
        "row_count": len(candidate["rows"]), "publication_state": state["state"],
        "writes_performed": False,
    }


def verify(root: Path = ROOT, *, collector: Callable[[Path], dict[str, Any]] = collect_federation) -> dict[str, Any]:
    with _suppress_bytecode_writes():
        candidate = collector(root)
        package = verify_package(root, candidate)
        state = publication_state(root, candidate)
    if state["state"] != "COMMITTED":
        raise FederationError(f"final RAM federation is not committed: {state['state']}")
    return {
        "status": "VERIFIED", "federation_fingerprint": candidate["fingerprint"],
        "row_count": len(candidate["rows"]), "package": package,
        "publication_state": state["state"], "writes_performed": False,
    }


def execute(root: Path = ROOT, *, collector: Callable[[Path], dict[str, Any]] = collect_federation) -> dict[str, Any]:
    with _suppress_bytecode_writes():
        candidate = collector(root)
        result = publish_candidate(root, candidate)
    return {
        "status": "PUBLISHED", "federation_fingerprint": candidate["fingerprint"],
        "row_count": len(candidate["rows"]), **result,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=str(ROOT))
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--preflight", action="store_true")
    action.add_argument("--verify", action="store_true")
    action.add_argument("--execute", action="store_true")
    args = parser.parse_args(argv)
    root = Path(args.root).resolve(strict=True)
    if args.preflight:
        result = preflight(root)
    elif args.verify:
        result = verify(root)
    else:
        result = execute(root)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FederationError, OSError) as error:
        raise SystemExit(f"ERROR: {error}") from error
