#!/usr/bin/env python3
"""Audit and, only on explicit request, evict redundant LAVA acquisition ZIPs.

The seven UKB v1.1 ZIP archives are acquisition containers.  The scientific LD
reference is the 44 extracted chromosome ``.info``/``.bcor`` payloads.  This
module preserves the distinction with a content-addressed, crash-resumable
plan -> intent -> final-receipt protocol.

Downstream code should call :func:`verify_reference_state`.  It accepts exactly
two usable states:

* ``ARCHIVES_PRESENT_FULLY_VERIFIED``; or
* ``ARCHIVES_INTENTIONALLY_EVICTED_WITH_FINAL_RECEIPT_AND_ALL_44_EXTRACTED_PAYLOADS_REHASHED``.

An interrupted eviction is never presented as a valid reference.  It can only
be resumed by this script against the already-fsynced intent record.
"""

from __future__ import annotations

import argparse
import ast
import csv
import fcntl
import hashlib
import importlib.util
import io
import json
import os
import re
import stat
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator


ROOT = Path(__file__).resolve().parents[1]
IMPLEMENTATION = Path("scripts/146_manage_lava_reference_archives.py")
POLICY = Path("config/lava_analysis_policy.json")
SOURCE_REGISTRY = Path("config/lava_reference_sources.tsv")
REFERENCE_ROOT = Path("ref/lava/ukb_v1.1")
REFERENCE_PROVENANCE = REFERENCE_ROOT / "reference.provenance.json"
DOWNLOAD_MANIFEST = REFERENCE_ROOT / "download_manifest.tsv"
EXTRACTED_MANIFEST = REFERENCE_ROOT / "extracted_manifest.tsv"
LEGACY_CONTRACT = Path("scripts/lava_contract.py")
SOURCE_LOCK_ROOT = Path("results/track_b/checkpoints/lava_continuations/source_locks")
CONTINUATION_RUN_ROOT = Path("results/track_b/checkpoints/lava_continuations/runs")
SUPERSESSION_LOCK_ROOT = Path("results/track_b/lava_continuations/supersessions/final_locks")
PLACO_TERMINAL_GATE = Path(
    "results/track_b/pleiotropy/continuations/post_lava_terminal_v2/terminal_gate.lock.json"
)

EVICTION_ROOT = Path("results/track_b/reference_archive_eviction")
PLAN_ROOT = EVICTION_ROOT / "plans"
INTENT_ROOT = EVICTION_ROOT / "intents"
RECEIPT_ROOT = EVICTION_ROOT / "receipts"
MUTEX = EVICTION_ROOT / ".archive_eviction.lock"

PLAN_SCHEMA = "track-b-lava-reference-archive-eviction-plan.1"
INTENT_SCHEMA = "track-b-lava-reference-archive-eviction-intent.1"
RECEIPT_SCHEMA = "track-b-lava-reference-archive-eviction-receipt.1"
SOURCE_LOCK_SCHEMA = "track-b-lava-discovery-family.1"
SOURCE_BASELINE_SCHEMA = 1
PINNED_SOURCE_FINGERPRINT = "caeb6b5a1188560f27609cad77801715e0bedfb998660a531c5491afad1177c7"
PINNED_CONTINUATION_FINGERPRINT = "c139e2368cae8f4b88c227a3f0a541b3ad988d941ce881a4b1c1aef9e28b8fc7"

ARCHIVES_PRESENT = "ARCHIVES_PRESENT_FULLY_VERIFIED"
ARCHIVES_EVICTED = (
    "ARCHIVES_INTENTIONALLY_EVICTED_WITH_FINAL_RECEIPT_AND_ALL_44_"
    "EXTRACTED_PAYLOADS_REHASHED"
)
PARTIAL_STATE = "PARTIAL_ARCHIVE_EVICTION_RESTART_REQUIRED"
CONFIRMATION = "EVICT_EXACTLY_SEVEN_VERIFIED_LAVA_ACQUISITION_ARCHIVES"

SOURCE_FIELDS = ["archive_id", "chromosomes", "url", "archive_bytes", "archive_filename"]
DOWNLOAD_FIELDS = SOURCE_FIELDS + ["sha256"]
EXTRACTED_FIELDS = ["chromosome", "file_type", "path", "bytes", "sha256"]
HEX64 = re.compile(r"^[0-9a-f]{64}$")
EXPECTED_ARCHIVE_COUNT = 7
EXPECTED_ARCHIVE_BYTES = 14_110_596_095
EXPECTED_EXTRACTED_COUNT = 44
MAX_STRUCTURED_ARTIFACT_BYTES = 64 * 1024**2
EXPECTED_REFERENCE_ROLES = {
    "COMMON_INPUT": 12,
    "REFERENCE_CHROMOSOME": 44,
    "CHROMOSOME_SHARD": 176,
    "CHROMOSOME_INPUT_INFO": 22,
}


class ArchiveEvictionError(RuntimeError):
    """A fail-closed archive/reference validation or lifecycle error."""


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest_json(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def _require_plan_id(value: str) -> str:
    if not HEX64.fullmatch(value):
        raise ArchiveEvictionError("plan ID must be one lowercase SHA-256 digest")
    return value


def plan_path(plan_id: str) -> Path:
    return PLAN_ROOT / f"{_require_plan_id(plan_id)}.json"


def intent_path(plan_id: str) -> Path:
    return INTENT_ROOT / f"{_require_plan_id(plan_id)}.json"


def receipt_path(plan_id: str) -> Path:
    return RECEIPT_ROOT / f"{_require_plan_id(plan_id)}.json"


def _relative(root: Path, path: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError as error:
        raise ArchiveEvictionError(f"path escapes repository: {path}") from error


def safe_path(root: Path, relative: Path | str, label: str) -> Path:
    relative = Path(relative)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ArchiveEvictionError(f"unsafe {label} path: {relative}")
    try:
        root_real = root.resolve(strict=True)
    except OSError as error:
        raise ArchiveEvictionError(f"repository root is unavailable: {error}") from error
    path = root / relative
    current = root
    for offset, part in enumerate(relative.parts):
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise ArchiveEvictionError(f"could not inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise ArchiveEvictionError(f"{label} contains a symbolic link: {current}")
        if offset < len(relative.parts) - 1 and not stat.S_ISDIR(observed.st_mode):
            raise ArchiveEvictionError(f"{label} has a non-directory ancestor: {current}")
    try:
        path.resolve(strict=False).relative_to(root_real)
    except (OSError, RuntimeError, ValueError) as error:
        raise ArchiveEvictionError(f"{label} escapes repository: {path}") from error
    return path


def stable_bytes(
    root: Path, relative: Path | str, *, allow_empty: bool = False,
) -> tuple[bytes, dict[str, object]]:
    path = safe_path(root, relative, "artifact")
    digest = hashlib.sha256()
    chunks: list[bytes] = []
    try:
        path_lstat = os.lstat(path)
        if (
            not stat.S_ISREG(path_lstat.st_mode)
            or path_lstat.st_size < 0
            or (path_lstat.st_size == 0 and not allow_empty)
        ):
            raise ArchiveEvictionError(f"artifact is not a real non-empty regular file: {path}")
        if path_lstat.st_size > MAX_STRUCTURED_ARTIFACT_BYTES:
            raise ArchiveEvictionError(f"structured artifact exceeds safe in-memory bound: {path}")
        flags = os.O_RDONLY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path, flags)
        try:
            before = os.fstat(descriptor)
            while True:
                block = os.read(descriptor, 4 * 1024 * 1024)
                if not block:
                    break
                chunks.append(block)
                digest.update(block)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        current = os.lstat(path)
    except ArchiveEvictionError:
        raise
    except OSError as error:
        raise ArchiveEvictionError(f"could not hash artifact {path}: {error}") from error
    identities = [
        (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns)
        for item in (path_lstat, before, after, current)
    ]
    if len(set(identities)) != 1 or not stat.S_ISREG(current.st_mode):
        raise ArchiveEvictionError(f"artifact changed while hashing: {path}")
    return b"".join(chunks), {
        "path": _relative(root, path),
        "bytes": int(current.st_size),
        "sha256": digest.hexdigest(),
    }


def stable_identity(root: Path, relative: Path | str, *, allow_empty: bool = False) -> dict[str, object]:
    # Avoid retaining multi-gigabyte payloads merely to compute their identity.
    path = safe_path(root, relative, "artifact")
    digest = hashlib.sha256()
    try:
        path_lstat = os.lstat(path)
        if (
            not stat.S_ISREG(path_lstat.st_mode)
            or path_lstat.st_size < 0
            or (path_lstat.st_size == 0 and not allow_empty)
        ):
            raise ArchiveEvictionError(f"artifact is not a real non-empty regular file: {path}")
        flags = os.O_RDONLY
        if hasattr(os, "O_CLOEXEC"):
            flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        descriptor = os.open(path, flags)
        try:
            before = os.fstat(descriptor)
            while True:
                block = os.read(descriptor, 4 * 1024 * 1024)
                if not block:
                    break
                digest.update(block)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        current = os.lstat(path)
    except ArchiveEvictionError:
        raise
    except OSError as error:
        raise ArchiveEvictionError(f"could not hash artifact {path}: {error}") from error
    identities = [
        (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns)
        for item in (path_lstat, before, after, current)
    ]
    if len(set(identities)) != 1 or not stat.S_ISREG(current.st_mode):
        raise ArchiveEvictionError(f"artifact changed while hashing: {path}")
    return {
        "path": _relative(root, path),
        "bytes": int(current.st_size),
        "sha256": digest.hexdigest(),
    }


def identity_matches(observed: object, expected: object) -> bool:
    return (
        isinstance(observed, dict)
        and isinstance(expected, dict)
        and all(observed.get(key) == expected.get(key) for key in ("path", "bytes", "sha256"))
    )


def read_json(root: Path, relative: Path | str) -> dict[str, Any]:
    encoded, _ = stable_bytes(root, relative)
    try:
        value = json.loads(encoded.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ArchiveEvictionError(f"unreadable JSON artifact {relative}: {error}") from error
    if not isinstance(value, dict):
        raise ArchiveEvictionError(f"JSON artifact must contain one object: {relative}")
    return value


def read_tsv(root: Path, relative: Path | str, expected_fields: list[str]) -> list[dict[str, str]]:
    encoded, _ = stable_bytes(root, relative)
    try:
        text = encoded.decode("utf-8")
        reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    except (UnicodeError, csv.Error) as error:
        raise ArchiveEvictionError(f"unreadable TSV artifact {relative}: {error}") from error
    if (
        fields != expected_fields
        or any(None in row or any(value is None for value in row.values()) for row in rows)
    ):
        raise ArchiveEvictionError(f"tabular schema or row shape drifted: {relative}")
    return rows


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def ensure_durable_directory(root: Path, relative: Path) -> Path:
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ArchiveEvictionError(f"unsafe durable-directory path: {relative}")
    current = root
    for part in relative.parts:
        candidate = current / part
        try:
            os.mkdir(candidate)
        except FileExistsError:
            try:
                observed = os.lstat(candidate)
            except OSError as error:
                raise ArchiveEvictionError(
                    f"could not inspect durable-directory component: {candidate}"
                ) from error
            if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
                raise ArchiveEvictionError(
                    f"durable-directory component is not a real directory: {candidate}"
                )
        except OSError as error:
            raise ArchiveEvictionError(
                f"could not create durable-directory component: {candidate}"
            ) from error
        else:
            # Persist every newly created directory entry, not only the final
            # plans/intents/receipts directory, before any archive can be unlinked.
            fsync_directory(current)
        current = candidate
    fsync_directory(current)
    return safe_path(root, relative, "durable directory")


def publish_no_replace(root: Path, relative: Path, payload: dict[str, object]) -> dict[str, object]:
    ensure_durable_directory(root, relative.parent)
    path = safe_path(root, relative, "no-replace publication")
    encoded = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise ArchiveEvictionError(f"refusing to replace existing artifact: {relative}") from error
        fsync_directory(path.parent)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
    return stable_identity(root, relative)


@contextmanager
def eviction_mutex(root: Path) -> Iterator[None]:
    ensure_durable_directory(root, MUTEX.parent)
    path = safe_path(root, MUTEX, "archive-eviction mutex")
    flags = os.O_RDWR | os.O_CREAT
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    descriptor = os.open(path, flags, 0o600)
    try:
        opened = os.fstat(descriptor)
        current = os.lstat(path)
        if (
            not stat.S_ISREG(opened.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (opened.st_dev, opened.st_ino) != (current.st_dev, current.st_ino)
            or opened.st_nlink != 1
            or current.st_nlink != 1
        ):
            raise ArchiveEvictionError("archive-eviction mutex is not one private regular file")
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise ArchiveEvictionError("another archive-eviction action holds the mutex") from error
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def _plain_identity(record: dict[str, object]) -> dict[str, object]:
    return {key: record[key] for key in ("path", "bytes", "sha256")}


def validate_reference_metadata(
    root: Path, *, require_archives: bool,
) -> dict[str, Any]:
    """Rebuild the sealed reference family and fully rehash scientific payloads."""

    policy = read_json(root, POLICY)
    if (
        policy.get("analysis_id") != "atlas-v1.0-lava"
        or policy.get("lava_version") != "0.1.5"
        or policy.get("reference_expected_archives") != EXPECTED_ARCHIVE_COUNT
        or policy.get("reference_expected_archive_bytes") != EXPECTED_ARCHIVE_BYTES
        or policy.get("reference_expected_extracted_files") != EXPECTED_EXTRACTED_COUNT
        or policy.get("reference_source_registry") != str(SOURCE_REGISTRY)
        or policy.get("reference_download_manifest") != str(DOWNLOAD_MANIFEST)
        or policy.get("reference_extracted_manifest") != str(EXTRACTED_MANIFEST)
        or policy.get("reference_payload_root") != str(REFERENCE_ROOT)
        or policy.get("reference_provenance") != str(REFERENCE_PROVENANCE)
    ):
        raise ArchiveEvictionError("LAVA reference policy identity or exact family scope drifted")

    source_identity = stable_identity(root, SOURCE_REGISTRY)
    if source_identity["sha256"] != policy.get("reference_source_registry_sha256"):
        raise ArchiveEvictionError("LAVA reference source registry differs from policy")
    sources = read_tsv(root, SOURCE_REGISTRY, SOURCE_FIELDS)
    downloads = read_tsv(root, DOWNLOAD_MANIFEST, DOWNLOAD_FIELDS)
    extracted_rows = read_tsv(root, EXTRACTED_MANIFEST, EXTRACTED_FIELDS)
    if (
        len(sources) != EXPECTED_ARCHIVE_COUNT
        or len(downloads) != EXPECTED_ARCHIVE_COUNT
        or len({row["archive_id"] for row in sources}) != EXPECTED_ARCHIVE_COUNT
    ):
        raise ArchiveEvictionError("source/download manifests do not define exactly seven archives")

    archive_records: list[dict[str, object]] = []
    recovery_records: list[dict[str, object]] = []
    for source, download in zip(sources, downloads, strict=True):
        if any(source[field] != download[field] for field in SOURCE_FIELDS):
            raise ArchiveEvictionError("download manifest differs from the ordered source registry")
        filename = Path(download["archive_filename"])
        if len(filename.parts) != 1 or filename.suffix != ".zip":
            raise ArchiveEvictionError(f"unsafe LAVA archive filename: {filename}")
        if not download["url"].startswith("https://") or not HEX64.fullmatch(download["sha256"]):
            raise ArchiveEvictionError(f"invalid recovery URL or SHA-256: {download['archive_id']}")
        try:
            size = int(download["archive_bytes"])
        except ValueError as error:
            raise ArchiveEvictionError("invalid archive byte count") from error
        relative = REFERENCE_ROOT / "archives" / filename
        record = {"path": str(relative), "bytes": size, "sha256": download["sha256"]}
        archive_records.append(record)
        recovery_records.append({
            "archive_id": download["archive_id"],
            "chromosomes": download["chromosomes"],
            "url": download["url"],
            **record,
        })
    if sum(int(record["bytes"]) for record in archive_records) != EXPECTED_ARCHIVE_BYTES:
        raise ArchiveEvictionError("seven-archive byte total drifted")

    expected_extracted = {
        (str(chromosome), suffix): f"lava-ukb-v1.1_chr{chromosome}.{suffix}"
        for chromosome in range(1, 23) for suffix in ("info", "bcor")
    }
    extracted_records: list[dict[str, object]] = []
    keys: set[tuple[str, str]] = set()
    for row in extracted_rows:
        key = (row["chromosome"], row["file_type"])
        if key in keys or expected_extracted.get(key) != row["path"]:
            raise ArchiveEvictionError("extracted manifest differs from the exact 22-chromosome pair family")
        keys.add(key)
        if not HEX64.fullmatch(row["sha256"]):
            raise ArchiveEvictionError(f"invalid extracted payload SHA-256: {row['path']}")
        try:
            size = int(row["bytes"])
        except ValueError as error:
            raise ArchiveEvictionError("invalid extracted payload byte count") from error
        extracted_records.append({
            "path": str(REFERENCE_ROOT / row["path"]),
            "bytes": size,
            "sha256": row["sha256"],
        })
    if keys != set(expected_extracted) or len(extracted_records) != EXPECTED_EXTRACTED_COUNT:
        raise ArchiveEvictionError("extracted manifest omits a required .info/.bcor payload")

    provenance = read_json(root, REFERENCE_PROVENANCE)
    expected_provenance = {
        "schema_version": "sleep-atlas-lava-reference.1",
        "analysis_id": policy["analysis_id"],
        "reference": policy["reference"],
        "reference_population": policy["reference_population"],
        "source_registry": str(SOURCE_REGISTRY),
        "source_registry_sha256": source_identity["sha256"],
        "download_manifest": str(DOWNLOAD_MANIFEST),
        "download_manifest_sha256": stable_identity(root, DOWNLOAD_MANIFEST)["sha256"],
        "extracted_manifest": str(EXTRACTED_MANIFEST),
        "extracted_manifest_sha256": stable_identity(root, EXTRACTED_MANIFEST)["sha256"],
        "archive_count": EXPECTED_ARCHIVE_COUNT,
        "archive_bytes": EXPECTED_ARCHIVE_BYTES,
        "extracted_file_count": EXPECTED_EXTRACTED_COUNT,
        "extracted_bytes": sum(int(record["bytes"]) for record in extracted_records),
        "archives": archive_records,
        "extracted_files": extracted_records,
        "verification": "SHA-256 verified after official HTTPS acquisition",
        "contract_script_sha256": stable_identity(root, LEGACY_CONTRACT)["sha256"],
    }
    if provenance != expected_provenance:
        raise ArchiveEvictionError("immutable reference provenance differs from manifests or contract")

    # These 44 hashes are always recomputed, including in post-eviction mode.
    live_payloads: list[dict[str, object]] = []
    for expected in extracted_records:
        observed = stable_identity(root, str(expected["path"]))
        if not identity_matches(observed, expected):
            raise ArchiveEvictionError(f"extracted reference payload differs: {expected['path']}")
        live_payloads.append(observed)

    live_archives: list[dict[str, object]] = []
    archive_presence: list[bool] = []
    for expected in archive_records:
        path = safe_path(root, str(expected["path"]), "archive")
        present = path.exists() or path.is_symlink()
        archive_presence.append(present)
        if present:
            observed = stable_identity(root, str(expected["path"]))
            if not identity_matches(observed, expected):
                raise ArchiveEvictionError(f"acquisition archive differs: {expected['path']}")
            live_archives.append(observed)
    if require_archives and not all(archive_presence):
        raise ArchiveEvictionError("all seven exact acquisition archives are required for this audit")

    return {
        "policy": stable_identity(root, POLICY),
        "source_registry": source_identity,
        "download_manifest": stable_identity(root, DOWNLOAD_MANIFEST),
        "extracted_manifest": stable_identity(root, EXTRACTED_MANIFEST),
        "reference_provenance_identity": stable_identity(root, REFERENCE_PROVENANCE),
        "reference_provenance": provenance,
        "archives": archive_records,
        "archive_presence": archive_presence,
        "live_archives": live_archives,
        "recovery": recovery_records,
        "extracted_payloads": live_payloads,
        "archive_family_sha256": digest_json(archive_records),
        "extracted_family_sha256": digest_json(live_payloads),
    }


def _all_strings(value: object) -> Iterator[str]:
    if isinstance(value, str):
        yield value
    elif isinstance(value, dict):
        for key, item in value.items():
            yield str(key)
            yield from _all_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from _all_strings(item)


def _archive_string(value: str, archive_paths: set[str]) -> bool:
    return value in archive_paths or "/archives/" in value or value.endswith(".zip")


def _superseded_fingerprints(root: Path) -> set[str]:
    directory = safe_path(root, SUPERSESSION_LOCK_ROOT, "supersession-lock root")
    if not directory.exists():
        return set()
    fingerprints: set[str] = set()
    for path in directory.glob("*.lock.json"):
        relative = Path(_relative(root, path))
        payload = read_json(root, relative)
        value = payload.get("superseded_continuation_fingerprint")
        if isinstance(value, str) and HEX64.fullmatch(value):
            fingerprints.add(value)
    return fingerprints


def _current_lineage(root: Path, source_lock_identity: dict[str, object]) -> tuple[str, dict[str, Any], dict[str, object]]:
    directory = safe_path(root, CONTINUATION_RUN_ROOT, "continuation run root")
    superseded = _superseded_fingerprints(root)
    candidates: list[tuple[str, dict[str, Any], dict[str, object]]] = []
    for lineage_path in directory.glob("*/lineage.json"):
        fingerprint = lineage_path.parent.name
        if not HEX64.fullmatch(fingerprint) or fingerprint in superseded:
            continue
        relative = Path(_relative(root, lineage_path))
        payload = read_json(root, relative)
        if (
            payload.get("analysis_id") != "track-b-v1.0-local"
            or payload.get("source_discovery_fingerprint") != PINNED_SOURCE_FINGERPRINT
            or payload.get("continuation_execution_fingerprint") != fingerprint
            or payload.get("contract_payload_sha256") != fingerprint
            or not identity_matches(payload.get("source_family_lock", {}), source_lock_identity)
        ):
            raise ArchiveEvictionError(f"continuation lineage identity drifted: {relative}")
        candidates.append((fingerprint, payload, stable_identity(root, relative)))
    if len(candidates) != 1:
        raise ArchiveEvictionError(
            f"expected one unsuperseded continuation lineage; observed {len(candidates)}"
        )
    if candidates[0][0] != PINNED_CONTINUATION_FINGERPRINT:
        raise ArchiveEvictionError(
            "unsuperseded continuation is not the pinned c139 production lineage"
        )
    return candidates[0]


def _published_continuation_receipts(root: Path, fingerprint: str) -> list[dict[str, object]]:
    attempts = safe_path(
        root, CONTINUATION_RUN_ROOT / fingerprint / ".attempts", "continuation attempt root",
    )
    if not attempts.is_dir():
        raise ArchiveEvictionError("current continuation attempt namespace is missing")
    result: list[dict[str, object]] = []
    for receipt_path in sorted(attempts.glob("*/receipt.json")):
        ready_path = receipt_path.parent / "READY"
        if not ready_path.is_file():
            continue
        receipt_relative = Path(_relative(root, receipt_path))
        ready_relative = Path(_relative(root, ready_path))
        receipt_identity = stable_identity(root, receipt_relative)
        ready = read_json(root, ready_relative)
        if (
            ready.get("state") != "READY"
            or ready.get("receipt_sha256") != receipt_identity["sha256"]
        ):
            raise ArchiveEvictionError(f"continuation READY/receipt mismatch: {receipt_relative}")
        receipt = read_json(root, receipt_relative)
        active = receipt.get("active_inputs")
        if not isinstance(active, list) or not active:
            raise ArchiveEvictionError(f"continuation receipt lacks active inputs: {receipt_relative}")
        result.append({
            "phase": receipt.get("phase"),
            "locus_index": receipt.get("locus_index"),
            "receipt": receipt_identity,
            "ready": stable_identity(root, ready_relative),
            "active_input_count": len(active),
            "active_input_paths": [item.get("path") for item in active if isinstance(item, dict)],
        })
    return result


def validate_lava_dependency_baseline(root: Path, reference: dict[str, Any]) -> dict[str, Any]:
    source_lock_relative = SOURCE_LOCK_ROOT / f"{PINNED_SOURCE_FINGERPRINT}.json"
    source_lock_identity = stable_identity(root, source_lock_relative)
    source_lock = read_json(root, source_lock_relative)
    checkpoints = source_lock.get("checkpoints")
    if (
        source_lock.get("schema_version") != SOURCE_LOCK_SCHEMA
        or source_lock.get("analysis_id") != "track-b-v1.0-local"
        or source_lock.get("source_discovery_fingerprint") != PINNED_SOURCE_FINGERPRINT
        or source_lock.get("checkpoint_count") != 2495
        or not isinstance(checkpoints, list)
        or len(checkpoints) != 2495
        or source_lock.get("checkpoint_family_sha256") != digest_json(checkpoints)
        or [item.get("locus_index") if isinstance(item, dict) else None for item in checkpoints]
        != list(range(1, 2496))
    ):
        raise ArchiveEvictionError("immutable LAVA source-family lock drifted")
    baseline_record = source_lock.get("active_input_baseline")
    if not isinstance(baseline_record, dict):
        raise ArchiveEvictionError("source lock lacks its active-input baseline identity")
    baseline_identity = stable_identity(root, str(baseline_record.get("path", "")))
    if not identity_matches(baseline_identity, baseline_record):
        raise ArchiveEvictionError("active-input baseline differs from source lock")
    baseline = read_json(root, str(baseline_record["path"]))
    entries = baseline.get("entries")
    if (
        baseline.get("schema_version") != SOURCE_BASELINE_SCHEMA
        or baseline.get("analysis_id") != "track-b-v1.0-local"
        or baseline.get("execution_fingerprint") != PINNED_SOURCE_FINGERPRINT
        or baseline.get("content_verification") != "FULL_SHA256_PREFLIGHT_BEFORE_BASELINE_PUBLICATION"
        or not isinstance(entries, list)
        or len(entries) != sum(EXPECTED_REFERENCE_ROLES.values())
    ):
        raise ArchiveEvictionError("source active-input baseline scope drifted")
    roles: dict[str, int] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            raise ArchiveEvictionError("source active-input baseline entry is malformed")
        role = str(entry.get("role", ""))
        roles[role] = roles.get(role, 0) + 1
    if roles != EXPECTED_REFERENCE_ROLES:
        raise ArchiveEvictionError(f"source active-input role family drifted: {roles}")

    archive_paths = {str(record["path"]) for record in reference["archives"]}
    offenders = sorted({value for value in _all_strings(baseline) if _archive_string(value, archive_paths)})
    offenders += sorted({value for value in _all_strings(source_lock) if _archive_string(value, archive_paths)})
    if offenders:
        raise ArchiveEvictionError(f"LAVA source baseline depends on acquisition archives: {offenders[:3]}")
    reference_entries = {
        (str(item.get("path")), int(item.get("bytes", -1)), str(item.get("sha256")))
        for item in entries if item.get("role") == "REFERENCE_CHROMOSOME"
    }
    extracted_entries = {
        (str(item["path"]), int(item["bytes"]), str(item["sha256"]))
        for item in reference["extracted_payloads"]
    }
    if reference_entries != extracted_entries:
        raise ArchiveEvictionError("source baseline is not bound to the exact 44 extracted payloads")
    if not any(
        item.get("role") == "COMMON_INPUT"
        and item.get("path") == str(REFERENCE_PROVENANCE)
        and int(item.get("bytes", -1)) == reference["reference_provenance_identity"]["bytes"]
        and item.get("sha256") == reference["reference_provenance_identity"]["sha256"]
        for item in entries
    ):
        raise ArchiveEvictionError("source baseline lacks the exact reference provenance identity")

    fingerprint, lineage, lineage_identity = _current_lineage(root, source_lock_identity)
    receipts = _published_continuation_receipts(root, fingerprint)
    receipt_offenders: list[str] = []
    for receipt in receipts:
        for value in receipt["active_input_paths"]:
            if isinstance(value, str) and _archive_string(value, archive_paths):
                receipt_offenders.append(value)
    if receipt_offenders:
        raise ArchiveEvictionError(
            f"current continuation receipt depends on acquisition archives: {receipt_offenders[:3]}"
        )
    return {
        "source_discovery_fingerprint": PINNED_SOURCE_FINGERPRINT,
        "source_family_lock": source_lock_identity,
        "active_input_baseline": baseline_identity,
        "baseline_entry_count": len(entries),
        "baseline_role_counts": roles,
        "archive_dependency_count": 0,
        "extracted_reference_dependency_count": len(reference_entries),
        "current_continuation_fingerprint": fingerprint,
        "current_lineage": lineage_identity,
        # Receipt count is deliberately omitted: a plan created while the
        # continuation is running must not become stale as more immutable
        # checkpoints appear.  Every currently published receipt was inspected;
        # terminal identity is captured later in the durable intent.
        "published_receipt_archive_dependency_count": 0,
        "proof": (
            "THE_CONTENT_LOCKED_SOURCE_BASELINE_AND_EVERY_CURRENT_PUBLISHED_"
            "CONTINUATION_ACTIVE_INPUT_RECEIPT_CONTAIN_ZERO_ARCHIVE_PATHS_AND_"
            "BIND_ALL_44_EXTRACTED_REFERENCE_PAYLOADS"
        ),
    }


def _normalized_additional_consumers(
    values: list[Path | str] | tuple[Path | str, ...],
) -> list[Path]:
    normalized: list[Path] = []
    for value in values:
        relative = Path(value)
        if (
            relative.is_absolute()
            or not relative.parts
            or ".." in relative.parts
            or relative.parts[0] != "scripts"
            or relative.suffix not in {".py", ".R"}
        ):
            raise ArchiveEvictionError(f"unsafe additional consumer path: {relative}")
        if relative not in normalized:
            normalized.append(relative)
    return sorted(normalized, key=str)


def _has_executable_guard_integration(relative: Path, source: str) -> bool:
    if relative.suffix == ".py":
        try:
            tree = ast.parse(source, filename=str(relative))
        except SyntaxError as error:
            raise ArchiveEvictionError(
                f"could not parse Python consumer while auditing guard integration: {relative}"
            ) from error
        guard_call = any(
            isinstance(node, ast.Call)
            and (
                isinstance(node.func, ast.Name)
                and node.func.id == "verify_reference_state"
                or isinstance(node.func, ast.Attribute)
                and node.func.attr == "verify_reference_state"
            )
            for node in ast.walk(tree)
        )
        literal_strings = {
            node.value for node in ast.walk(tree)
            if isinstance(node, ast.Constant) and isinstance(node.value, str)
        }
        return (
            guard_call
            and any("146_manage_lava_reference_archives.py" in value for value in literal_strings)
            and ARCHIVES_PRESENT in literal_strings
            and ARCHIVES_EVICTED in literal_strings
        )
    # No R consumer currently owns this state transition.  This fallback keeps
    # the extension API usable if one is introduced, while requiring an actual
    # function-call spelling rather than accepting comments alone.
    uncommented = "\n".join(line.split("#", 1)[0] for line in source.splitlines())
    return (
        "verify_reference_state(" in uncommented
        and "146_manage_lava_reference_archives.py" in uncommented
        and ARCHIVES_PRESENT in uncommented
        and ARCHIVES_EVICTED in uncommented
    )


def audit_downstream_code(
    root: Path, additional_consumers: list[Path | str] | tuple[Path | str, ...] = (),
) -> dict[str, Any]:
    paths = {
        "legacy_reference_contract": Path("scripts/lava_contract.py"),
        "track_b_legacy_contract": Path("scripts/119_track_b_lava_contract.py"),
        "lava_v2_contract": Path("scripts/134_track_b_lava_continuation_contract.py"),
        "placo_terminal_gate": Path("scripts/139_build_track_b_placo_terminal_gate_v2.py"),
        "placo_bridge": Path("scripts/140_materialize_track_b_placo_pair_v2.py"),
        "placo_worker": Path("scripts/141_run_track_b_placo_pair_v2.R"),
        "placo_coordinator": Path("scripts/143_run_track_b_placo_sequential_v2.py"),
        "finemapping_gate": Path("scripts/145_build_track_b_finemapping_continuation_gate.py"),
        "finemapping_contract": Path("scripts/fine_mapping_contract.py"),
        "finemapping_materializer": Path("scripts/58_materialize_finemapping_locus.py"),
    }
    text: dict[str, str] = {}
    identities: dict[str, dict[str, object]] = {}
    for name, relative in paths.items():
        encoded, identity = stable_bytes(root, relative)
        try:
            text[name] = encoded.decode("utf-8")
        except UnicodeError as error:
            raise ArchiveEvictionError(f"downstream script is not UTF-8: {relative}") from error
        identities[name] = identity
    required = {
        "legacy_reference_contract": ["require_live=True", "archives + payloads"],
        "track_b_legacy_contract": ["lava_contract.validate_reference"],
        "lava_v2_contract": ["_require_live_source_fingerprint", "legacy_contract().run_fingerprint()"],
        "placo_terminal_gate": ["validate_source_bundles=True", "def deep_lava_terminal"],
        "placo_bridge": ["GATE.verify_gate()"],
        "placo_coordinator": ["GATE.verify_gate()"],
        "finemapping_contract": ["lava_contract.validate_reference"],
        "finemapping_materializer": [".info", ".bcor"],
    }
    for name, snippets in required.items():
        if any(snippet not in text[name] for snippet in snippets):
            raise ArchiveEvictionError(f"downstream archive-dependency audit assumption drifted: {name}")
    direct_placo_archive_literals = {
        name: [
            line.strip() for line in text[name].splitlines()
            if "/archives/" in line or ".zip" in line.lower()
        ]
        for name in ("placo_bridge", "placo_worker", "placo_coordinator")
    }
    if any(direct_placo_archive_literals.values()):
        raise ArchiveEvictionError("a PLACO execution script directly references the LAVA ZIP family")
    accepted_state_tokens = (ARCHIVES_PRESENT, ARCHIVES_EVICTED)
    integration = {
        name: _has_executable_guard_integration(paths[name], text[name])
        for name in ("placo_terminal_gate", "placo_coordinator", "finemapping_gate")
    }
    fixed_paths = set(paths.values())
    additional_paths = _normalized_additional_consumers(additional_consumers)
    if any(relative in fixed_paths or relative == IMPLEMENTATION for relative in additional_paths):
        raise ArchiveEvictionError("an additional consumer duplicates a core audited script")
    additional_identities: dict[str, dict[str, object]] = {}
    additional_integration: dict[str, bool] = {}
    for relative in additional_paths:
        encoded, identity = stable_bytes(root, relative)
        try:
            consumer_text = encoded.decode("utf-8")
        except UnicodeError as error:
            raise ArchiveEvictionError(f"additional consumer is not UTF-8: {relative}") from error
        if any(
            "/archives/" in line or ".zip" in line.lower()
            for line in consumer_text.splitlines()
        ):
            raise ArchiveEvictionError(
                f"additional scientific consumer directly references LAVA ZIPs: {relative}"
            )
        key = str(relative)
        additional_identities[key] = identity
        additional_integration[key] = _has_executable_guard_integration(
            relative, consumer_text,
        )
    compatible = all(integration.values()) and all(additional_integration.values())
    return {
        "script_identities": identities,
        "additional_consumer_audit": {
            "paths": [str(path) for path in additional_paths],
            "script_identities": additional_identities,
            "integrations": additional_integration,
            "extension_interface": (
                "PASS_EACH_NEW_POST_LAVA_PLACO_OR_FINEMAPPING_CONSUMER_VIA_"
                "--additional-consumer_BEFORE_PUBLISHING_THE_REAL_PLAN"
            ),
        },
        "placo_scientific_input": {
            "archive_bytes_required": False,
            "status": "NO_DIRECT_ZIP_DEPENDENCY_IN_FULL_P_MATERIALIZATION_OR_TESTING",
        },
        "placo_current_operational_verifier": {
            "archive_bytes_required": not (
                integration["placo_terminal_gate"] and integration["placo_coordinator"]
            ),
            "status": (
                "EVICTION_RECEIPT_API_INTEGRATED"
                if integration["placo_terminal_gate"] and integration["placo_coordinator"]
                else "BLOCKED_AFTER_EVICTION_UNTIL_SCRIPT_139_DEEP_GATE_ACCEPTS_THE_"
                "FINAL_EVICTION_RECEIPT;SCRIPTS_140_AND_143_CALL_THAT_GATE_REPEATEDLY"
            ),
        },
        "finemapping_scientific_input": {
            "archive_bytes_required": False,
            "required_payloads": "FULL_LOCUS_INFO_AND_BCOR_ONLY",
        },
        "finemapping_current_operational_verifier": {
            "archive_bytes_required": not integration["finemapping_gate"],
            "status": (
                "EVICTION_RECEIPT_API_INTEGRATED_IN_ADDITIVE_GATE"
                if integration["finemapping_gate"]
                else "BLOCKED_AFTER_EVICTION_UNTIL_FINE_MAPPING_REFERENCE_VALIDATION_"
                "ACCEPTS_THE_FINAL_EVICTION_RECEIPT"
            ),
        },
        "eviction_compatibility": {
            "state": (
                "READY_BOTH_REFERENCE_STATES_ACCEPTED_BY_ACTIVE_CONSUMER_GATES"
                if compatible else "BLOCKED_ACTIVE_CONSUMER_GATES_DO_NOT_ACCEPT_FINAL_RECEIPT"
            ),
            "integrations": integration,
            "additional_consumer_integrations": additional_integration,
            "required_module": "scripts/146_manage_lava_reference_archives.py",
            "required_api": "verify_reference_state(root: pathlib.Path) -> dict",
            "accepted_states": list(accepted_state_tokens),
        },
        "legacy_verifier_after_eviction": "WILL_FAIL_BY_DESIGN_BECAUSE_REQUIRE_LIVE_INCLUDES_ARCHIVES",
    }


def plan_payload(
    root: Path, additional_consumers: list[Path | str] | tuple[Path | str, ...] = (),
) -> dict[str, Any]:
    reference = validate_reference_metadata(root, require_archives=True)
    baseline = validate_lava_dependency_baseline(root, reference)
    downstream = audit_downstream_code(root, additional_consumers)
    return {
        "schema_version": PLAN_SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "state": "AUDITED_NO_FILES_REMOVED",
        "scope": {
            "archive_count": EXPECTED_ARCHIVE_COUNT,
            "archive_bytes": EXPECTED_ARCHIVE_BYTES,
            "scientific_payload_count": EXPECTED_EXTRACTED_COUNT,
            "scientific_payload_bytes": sum(
                int(item["bytes"]) for item in reference["extracted_payloads"]
            ),
            "removal_allowlist": [item["path"] for item in reference["archives"]],
            "all_other_files_forbidden": True,
        },
        "implementation": stable_identity(root, IMPLEMENTATION),
        "immutable_reference": {
            key: reference[key]
            for key in (
                "policy", "source_registry", "download_manifest", "extracted_manifest",
                "reference_provenance_identity", "archive_family_sha256",
                "extracted_family_sha256", "archives", "extracted_payloads", "recovery",
            )
        },
        "lava_dependency_audit": baseline,
        "downstream_dependency_audit": downstream,
        "eviction_prerequisites": {
            "current_continuation_terminal_receipt_required": True,
            "placo_terminal_gate_deeply_frozen_before_eviction_required": True,
            "current_operational_verifier_breakage_acknowledgement_required": True,
            "durable_intent_before_first_unlink_required": True,
            "all_seven_archives_rehashed_immediately_before_intent": True,
            "all_44_payloads_rehashed_before_and_after_eviction": True,
        },
        "scientific_assessment": (
            "SCIENTIFIC_CONTENT_SAFE_BECAUSE_ARCHIVES_ARE_REDUNDANT_ACQUISITION_"
            "CONTAINERS;OPERATIONALLY_GATED_UNTIL_CONSUMERS_ACCEPT_FINAL_RECEIPT"
        ),
        "legacy_verifier_warning": (
            "scripts/lava_contract.py --verify-reference WILL FAIL AFTER INTENTIONAL_"
            "EVICTION;THIS IS EXPECTED AND MUST NOT BE MISREPORTED_AS_LD_PAYLOAD_LOSS"
        ),
        "recovery_protocol": (
            "DOWNLOAD_EACH_EXACT_HTTPS_URL_TO_A_NEW_PART_FILE;VERIFY_RECORDED_BYTE_"
            "COUNT_AND_SHA256;ATOMically_INSTALL_AT_THE_RECORDED_ARCHIVE_PATH;DO_NOT_"
            "REEXTRACT_OVER_THE_VERIFIED_INFO_BCOR_PAYLOADS"
        ),
    }


def create_or_verify_plan(
    root: Path = ROOT,
    additional_consumers: list[Path | str] | tuple[Path | str, ...] = (),
) -> dict[str, Any]:
    expected = plan_payload(root, additional_consumers)
    plan_id = digest_json(expected)
    relative = plan_path(plan_id)
    path = safe_path(root, relative, "archive-eviction plan")
    if path.exists() or path.is_symlink():
        observed = read_json(root, relative)
        if observed != expected:
            raise ArchiveEvictionError("existing no-replace archive-eviction plan differs")
        return observed
    publish_no_replace(root, relative, expected)
    observed = read_json(root, relative)
    if observed != expected:
        raise ArchiveEvictionError("published archive-eviction plan differs")
    return observed


def _deep_verify_placo_gate(root: Path, expected: dict[str, Any]) -> str:
    """Run the production gate's own verifier before the irreversible boundary."""

    if root != ROOT:
        # Test fixtures cannot import repository-global modules.  Their exact
        # structural and identity checks still exercise the deletion protocol.
        return "TEST_FIXTURE_STRUCTURAL_VALIDATION"
    gate_script = root / "scripts/139_build_track_b_placo_terminal_gate_v2.py"
    specification = importlib.util.spec_from_file_location(
        "track_b_placo_terminal_gate_for_archive_eviction", gate_script,
    )
    if specification is None or specification.loader is None:
        raise ArchiveEvictionError("could not load the PLACO terminal-gate verifier")
    module = importlib.util.module_from_spec(specification)
    try:
        specification.loader.exec_module(module)
        observed = module.verify_gate()
    except BaseException as error:
        raise ArchiveEvictionError(
            f"PLACO terminal gate did not pass its native deep verifier: {error}"
        ) from error
    if observed != expected:
        raise ArchiveEvictionError("PLACO native deep verifier returned different gate content")
    return "PRODUCTION_NATIVE_GATE_VERIFY_PASS_BEFORE_DURABLE_INTENT"


def _terminal_binding(root: Path, baseline: dict[str, Any]) -> dict[str, Any]:
    receipts = _published_continuation_receipts(
        root, str(baseline["current_continuation_fingerprint"]),
    )
    terminal_receipts = [item for item in receipts if item["phase"] == "terminal-qc"]
    if len(terminal_receipts) != 1:
        raise ArchiveEvictionError("current LAVA continuation is not terminal; eviction is forbidden")
    terminal = terminal_receipts[0]
    gate_path = safe_path(root, PLACO_TERMINAL_GATE, "PLACO terminal gate")
    if not gate_path.is_file():
        raise ArchiveEvictionError(
            "PLACO terminal gate must be deeply frozen while archives remain present"
        )
    gate = read_json(root, PLACO_TERMINAL_GATE)
    lava = gate.get("lava_terminal")
    scientific_evidence = gate.get("lava_scientific_evidence")
    if (
        gate.get("schema_version") != "sleep-atlas-track-b-pleiotropy-terminal-gate.2"
        or not isinstance(lava, dict)
        or not isinstance(scientific_evidence, dict)
        or lava.get("source_discovery_fingerprint") != PINNED_SOURCE_FINGERPRINT
        or lava.get("continuation_execution_fingerprint")
        != baseline["current_continuation_fingerprint"]
        or lava.get("receipt") != terminal["receipt"]
        or lava.get("ready") != terminal["ready"]
        or scientific_evidence.get("status") != "FAILED_QC_NOT_CONSUMED"
    ):
        raise ArchiveEvictionError("PLACO terminal gate does not bind the current terminal LAVA lineage")
    gate_identity_before = stable_identity(root, PLACO_TERMINAL_GATE)
    validation = _deep_verify_placo_gate(root, gate)
    gate_identity_after = stable_identity(root, PLACO_TERMINAL_GATE)
    if gate_identity_after != gate_identity_before:
        raise ArchiveEvictionError("PLACO terminal gate changed during native deep verification")
    return {
        "source_discovery_fingerprint": PINNED_SOURCE_FINGERPRINT,
        "continuation_execution_fingerprint": baseline["current_continuation_fingerprint"],
        "source_family_lock": baseline["source_family_lock"],
        "active_input_baseline": baseline["active_input_baseline"],
        "current_lineage": baseline["current_lineage"],
        "terminal_receipt": terminal["receipt"],
        "terminal_ready": terminal["ready"],
        "placo_terminal_gate": gate_identity_after,
        "placo_terminal_gate_validation": validation,
    }


def _intent_payload(
    root: Path, plan_id: str, plan: dict[str, Any], reference: dict[str, Any],
    baseline: dict[str, Any],
) -> dict[str, Any]:
    compatibility = plan.get("downstream_dependency_audit", {}).get(
        "eviction_compatibility", {}
    )
    if compatibility.get("state") != "READY_BOTH_REFERENCE_STATES_ACCEPTED_BY_ACTIVE_CONSUMER_GATES":
        raise ArchiveEvictionError(
            "active PLACO/fine-mapping consumer gates do not yet accept the final eviction receipt"
        )
    return {
        "schema_version": INTENT_SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "state": "DURABLE_INTENT_ALL_SEVEN_ARCHIVES_VERIFIED_BEFORE_FIRST_UNLINK",
        "plan_id": plan_id,
        "plan": stable_identity(root, plan_path(plan_id)),
        "implementation": plan["implementation"],
        "pre_unlink_archive_family": reference["live_archives"],
        "pre_unlink_extracted_family": reference["extracted_payloads"],
        "archive_family_sha256": reference["archive_family_sha256"],
        "extracted_family_sha256": reference["extracted_family_sha256"],
        "terminal_binding": _terminal_binding(root, baseline),
        "removal_allowlist": [item["path"] for item in reference["archives"]],
        "recovery": reference["recovery"],
        "legacy_verifier_breakage_explicitly_acknowledged": True,
    }


def _validate_intent(
    root: Path, plan_id: str, plan: dict[str, Any], intent: dict[str, Any],
) -> None:
    immutable = plan.get("immutable_reference", {})
    terminal = intent.get("terminal_binding")
    if (
        set(intent) != {
            "schema_version", "analysis_id", "state", "plan_id", "plan", "implementation",
            "pre_unlink_archive_family", "pre_unlink_extracted_family",
            "archive_family_sha256", "extracted_family_sha256", "terminal_binding",
            "removal_allowlist", "recovery",
            "legacy_verifier_breakage_explicitly_acknowledged",
        }
        or plan.get("downstream_dependency_audit", {}).get("eviction_compatibility", {}).get("state")
        != "READY_BOTH_REFERENCE_STATES_ACCEPTED_BY_ACTIVE_CONSUMER_GATES"
        or intent.get("schema_version") != INTENT_SCHEMA
        or intent.get("analysis_id") != "track-b-v1.0-local"
        or intent.get("state")
        != "DURABLE_INTENT_ALL_SEVEN_ARCHIVES_VERIFIED_BEFORE_FIRST_UNLINK"
        or intent.get("plan_id") != plan_id
        or intent.get("plan") != stable_identity(root, plan_path(plan_id))
        or intent.get("implementation") != plan.get("implementation")
        or intent.get("pre_unlink_archive_family") != immutable.get("archives")
        or intent.get("pre_unlink_extracted_family") != immutable.get("extracted_payloads")
        or intent.get("archive_family_sha256") != immutable.get("archive_family_sha256")
        or intent.get("extracted_family_sha256") != immutable.get("extracted_family_sha256")
        or intent.get("removal_allowlist") != plan.get("scope", {}).get("removal_allowlist")
        or intent.get("recovery") != immutable.get("recovery")
        or intent.get("legacy_verifier_breakage_explicitly_acknowledged") is not True
        or not isinstance(terminal, dict)
        or set(terminal) != {
            "source_discovery_fingerprint", "continuation_execution_fingerprint",
            "source_family_lock", "active_input_baseline", "current_lineage",
            "terminal_receipt", "terminal_ready", "placo_terminal_gate",
            "placo_terminal_gate_validation",
        }
        or terminal.get("source_discovery_fingerprint") != PINNED_SOURCE_FINGERPRINT
        or terminal.get("continuation_execution_fingerprint")
        != plan.get("lava_dependency_audit", {}).get("current_continuation_fingerprint")
        or terminal.get("source_family_lock")
        != plan.get("lava_dependency_audit", {}).get("source_family_lock")
        or terminal.get("active_input_baseline")
        != plan.get("lava_dependency_audit", {}).get("active_input_baseline")
        or terminal.get("current_lineage")
        != plan.get("lava_dependency_audit", {}).get("current_lineage")
        or terminal.get("placo_terminal_gate_validation") not in {
            "TEST_FIXTURE_STRUCTURAL_VALIDATION",
            "PRODUCTION_NATIVE_GATE_VERIFY_PASS_BEFORE_DURABLE_INTENT",
        }
    ):
        raise ArchiveEvictionError("durable archive-eviction intent drifted")
    live_bindings = {
        "source_family_lock": SOURCE_LOCK_ROOT / f"{PINNED_SOURCE_FINGERPRINT}.json",
        "active_input_baseline": Path(
            str(plan["lava_dependency_audit"]["active_input_baseline"]["path"])
        ),
        "current_lineage": Path(str(plan["lava_dependency_audit"]["current_lineage"]["path"])),
        "terminal_receipt": Path(str(terminal.get("terminal_receipt", {}).get("path", ""))),
        "terminal_ready": Path(str(terminal.get("terminal_ready", {}).get("path", ""))),
        "placo_terminal_gate": PLACO_TERMINAL_GATE,
    }
    for key, relative in live_bindings.items():
        recorded = terminal.get(key)
        if not isinstance(recorded, dict) or stable_identity(root, relative) != recorded:
            raise ArchiveEvictionError(f"durable intent terminal binding drifted: {key}")


def _identity_shape(value: object, *, expected_path: Path | None = None) -> bool:
    if (
        not isinstance(value, dict)
        or set(value) != {"path", "bytes", "sha256"}
        or not isinstance(value.get("path"), str)
        or not isinstance(value.get("bytes"), int)
        or isinstance(value.get("bytes"), bool)
        or int(value["bytes"]) < 0
        or not isinstance(value.get("sha256"), str)
        or not HEX64.fullmatch(str(value["sha256"]))
    ):
        return False
    relative = Path(str(value["path"]))
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        return False
    return expected_path is None or relative == expected_path


def _plan_additional_consumers(plan: dict[str, Any]) -> list[Path]:
    downstream = plan.get("downstream_dependency_audit")
    audit = downstream.get("additional_consumer_audit") if isinstance(downstream, dict) else None
    paths = audit.get("paths") if isinstance(audit, dict) else None
    if not isinstance(paths, list) or any(not isinstance(item, str) for item in paths):
        raise ArchiveEvictionError("archive-eviction plan lacks its additional-consumer inventory")
    normalized = _normalized_additional_consumers(paths)
    if [str(path) for path in normalized] != paths:
        raise ArchiveEvictionError("archive-eviction plan additional-consumer inventory drifted")
    return normalized


def _validate_plan_static(
    root: Path, plan_id: str, plan: dict[str, Any], *, require_live_environment: bool,
) -> None:
    plan_id = _require_plan_id(plan_id)
    scope = plan.get("scope")
    immutable = plan.get("immutable_reference")
    implementation = plan.get("implementation")
    if (
        digest_json(plan) != plan_id
        or set(plan) != {
            "schema_version", "analysis_id", "state", "scope", "implementation",
            "immutable_reference", "lava_dependency_audit", "downstream_dependency_audit",
            "eviction_prerequisites", "scientific_assessment", "legacy_verifier_warning",
            "recovery_protocol",
        }
        or plan.get("schema_version") != PLAN_SCHEMA
        or plan.get("analysis_id") != "track-b-v1.0-local"
        or plan.get("state") != "AUDITED_NO_FILES_REMOVED"
        or not isinstance(scope, dict)
        or set(scope) != {
            "archive_count", "archive_bytes", "scientific_payload_count",
            "scientific_payload_bytes", "removal_allowlist", "all_other_files_forbidden",
        }
        or scope.get("archive_count") != EXPECTED_ARCHIVE_COUNT
        or scope.get("archive_bytes") != EXPECTED_ARCHIVE_BYTES
        or scope.get("scientific_payload_count") != EXPECTED_EXTRACTED_COUNT
        or scope.get("all_other_files_forbidden") is not True
        or not isinstance(immutable, dict)
        or set(immutable) != {
            "policy", "source_registry", "download_manifest", "extracted_manifest",
            "reference_provenance_identity", "archive_family_sha256",
            "extracted_family_sha256", "archives", "extracted_payloads", "recovery",
        }
        or not _identity_shape(implementation, expected_path=IMPLEMENTATION)
    ):
        raise ArchiveEvictionError("archive-eviction plan schema, ID, or scope drifted")

    archives = immutable.get("archives")
    payloads = immutable.get("extracted_payloads")
    recovery = immutable.get("recovery")
    if (
        not isinstance(archives, list)
        or len(archives) != EXPECTED_ARCHIVE_COUNT
        or not isinstance(payloads, list)
        or len(payloads) != EXPECTED_EXTRACTED_COUNT
        or not isinstance(recovery, list)
        or len(recovery) != EXPECTED_ARCHIVE_COUNT
        or any(not _identity_shape(item) for item in archives + payloads)
        or immutable.get("archive_family_sha256") != digest_json(archives)
        or immutable.get("extracted_family_sha256") != digest_json(payloads)
        or scope.get("removal_allowlist") != [item["path"] for item in archives]
        or scope.get("scientific_payload_bytes")
        != sum(int(item["bytes"]) for item in payloads)
        or len({item["path"] for item in archives}) != EXPECTED_ARCHIVE_COUNT
        or len({item["path"] for item in payloads}) != EXPECTED_EXTRACTED_COUNT
    ):
        raise ArchiveEvictionError("archive-eviction plan reference family drifted")
    for archive, recovery_item in zip(archives, recovery, strict=True):
        archive_path = Path(str(archive["path"]))
        if (
            archive_path.parent != REFERENCE_ROOT / "archives"
            or archive_path.suffix != ".zip"
            or not isinstance(recovery_item, dict)
            or set(recovery_item) != {
                "archive_id", "chromosomes", "url", "path", "bytes", "sha256",
            }
            or not isinstance(recovery_item.get("url"), str)
            or not str(recovery_item["url"]).startswith("https://")
            or not identity_matches(recovery_item, archive)
        ):
            raise ArchiveEvictionError("archive-eviction plan recovery family drifted")
    if any(
        not _identity_shape(immutable.get(key))
        for key in (
            "policy", "source_registry", "download_manifest", "extracted_manifest",
            "reference_provenance_identity",
        )
    ):
        raise ArchiveEvictionError("archive-eviction plan provenance identity drifted")

    additional_consumers = _plan_additional_consumers(plan)

    if require_live_environment:
        if not identity_matches(stable_identity(root, IMPLEMENTATION), implementation):
            raise ArchiveEvictionError("archive-eviction implementation differs from plan")
        if audit_downstream_code(root, additional_consumers) != plan.get("downstream_dependency_audit"):
            raise ArchiveEvictionError("downstream consumer code differs from archive-eviction plan")


def _load_plan(
    root: Path, plan_id: str, *, require_live_environment: bool,
) -> dict[str, Any]:
    plan_id = _require_plan_id(plan_id)
    plan = read_json(root, plan_path(plan_id))
    _validate_plan_static(
        root, plan_id, plan, require_live_environment=require_live_environment,
    )
    return plan


def _lifecycle_ids(root: Path, directory: Path, label: str) -> list[str]:
    path = safe_path(root, directory, label)
    if not path.exists():
        return []
    if not path.is_dir():
        raise ArchiveEvictionError(f"{label} is not a directory")
    identifiers: list[str] = []
    for child in sorted(path.iterdir(), key=lambda item: item.name):
        # A process death during no-replace publication can leave only its hidden
        # temporary file.  Such a file is not a published lifecycle artifact.
        if child.name.startswith("."):
            continue
        match = re.fullmatch(r"([0-9a-f]{64})\.json", child.name)
        try:
            observed = os.lstat(child)
        except OSError as error:
            raise ArchiveEvictionError(f"could not inspect {label} entry: {child}") from error
        if match is None or not stat.S_ISREG(observed.st_mode):
            raise ArchiveEvictionError(f"unexpected {label} entry: {child.name}")
        identifiers.append(match.group(1))
    return identifiers


def _preflight_remaining_archives(
    root: Path, expected_archives: list[dict[str, object]], *, allow_absent: bool,
) -> list[dict[str, object]]:
    """Fully hash and hard-link-check the entire remaining family before any unlink."""

    verified: list[dict[str, object]] = []
    for expected in expected_archives:
        relative = Path(str(expected["path"]))
        path = safe_path(root, relative, "archive pre-unlink target")
        present = path.exists() or path.is_symlink()
        if not present:
            if allow_absent:
                continue
            raise ArchiveEvictionError(f"archive is absent before durable intent: {relative}")
        observed = stable_identity(root, relative)
        try:
            current = os.lstat(path)
        except OSError as error:
            raise ArchiveEvictionError(f"could not inspect archive before intent: {relative}") from error
        if (
            not identity_matches(observed, expected)
            or not stat.S_ISREG(current.st_mode)
            or current.st_nlink != 1
        ):
            raise ArchiveEvictionError(f"archive is not one exact private regular file: {relative}")
        verified.append(observed)
    archive_directory = safe_path(root, REFERENCE_ROOT / "archives", "archive directory")
    fsync_directory(archive_directory)
    return verified


def _unlink_exact_archive(root: Path, expected: dict[str, object]) -> None:
    relative = Path(str(expected["path"]))
    if relative.parent != REFERENCE_ROOT / "archives" or relative.suffix != ".zip":
        raise ArchiveEvictionError(f"refusing non-allowlisted unlink target: {relative}")
    archive_dir = safe_path(root, relative.parent, "archive directory")
    name = relative.name
    dir_flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        dir_flags |= os.O_DIRECTORY
    directory_fd = os.open(archive_dir, dir_flags)
    file_fd: int | None = None
    try:
        file_flags = os.O_RDONLY
        if hasattr(os, "O_CLOEXEC"):
            file_flags |= os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            file_flags |= os.O_NOFOLLOW
        file_fd = os.open(name, file_flags, dir_fd=directory_fd)
        opened = os.fstat(file_fd)
        named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if (
            not stat.S_ISREG(opened.st_mode)
            or not stat.S_ISREG(named.st_mode)
            or (opened.st_dev, opened.st_ino) != (named.st_dev, named.st_ino)
            or opened.st_nlink != 1
            or named.st_nlink != 1
            or opened.st_size != int(expected["bytes"])
        ):
            raise ArchiveEvictionError(f"archive is not one exact private regular file: {relative}")
        digest = hashlib.sha256()
        while True:
            block = os.read(file_fd, 4 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
        after = os.fstat(file_fd)
        named_after = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if (
            (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns, after.st_ctime_ns)
            != (opened.st_dev, opened.st_ino, opened.st_size, opened.st_mtime_ns, opened.st_ctime_ns)
            or (
                named_after.st_dev, named_after.st_ino, named_after.st_size,
                named_after.st_mtime_ns, named_after.st_ctime_ns, named_after.st_nlink,
            ) != (
                opened.st_dev, opened.st_ino, opened.st_size,
                opened.st_mtime_ns, opened.st_ctime_ns, opened.st_nlink,
            )
            or digest.hexdigest() != expected["sha256"]
        ):
            raise ArchiveEvictionError(f"archive changed before unlink: {relative}")
        os.unlink(name, dir_fd=directory_fd)
        os.fsync(directory_fd)
        unlinked = os.fstat(file_fd)
        if unlinked.st_nlink != 0:
            raise ArchiveEvictionError(f"opened archive inode was not unlinked: {relative}")
        try:
            os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise ArchiveEvictionError(f"archive path remains after unlink: {relative}")
    except OSError as error:
        raise ArchiveEvictionError(f"could not safely unlink archive {relative}: {error}") from error
    finally:
        if file_fd is not None:
            os.close(file_fd)
        os.close(directory_fd)


def execute_eviction(
    *, root: Path = ROOT, plan_id: str, confirmation: str,
    acknowledge_operational_breakage: bool,
) -> dict[str, Any]:
    plan_id = _require_plan_id(plan_id)
    if confirmation != CONFIRMATION or not acknowledge_operational_breakage:
        raise ArchiveEvictionError("explicit exact eviction confirmation and verifier-breakage acknowledgement required")
    with eviction_mutex(root):
        receipt_ids = _lifecycle_ids(root, RECEIPT_ROOT, "archive-eviction receipt root")
        if receipt_ids:
            if receipt_ids != [plan_id]:
                raise ArchiveEvictionError("a different or ambiguous final eviction receipt exists")
            return verify_final_receipt(root, plan_id)

        plan = _load_plan(root, plan_id, require_live_environment=True)
        intent_ids = _lifecycle_ids(root, INTENT_ROOT, "archive-eviction intent root")
        if any(identifier != plan_id for identifier in intent_ids) or len(intent_ids) > 1:
            raise ArchiveEvictionError("a different or ambiguous durable eviction intent exists")
        target_intent = safe_path(root, intent_path(plan_id), "archive-eviction intent")
        if not (target_intent.exists() or target_intent.is_symlink()):
            # Nothing may be absent before the durable intent.  All seven archive
            # bytes and all 44 payload bytes are rehashed immediately beforehand.
            reference = validate_reference_metadata(root, require_archives=True)
            baseline = validate_lava_dependency_baseline(root, reference)
            current_plan = plan_payload(root, _plan_additional_consumers(plan))
            if current_plan != plan:
                raise ArchiveEvictionError("live eviction audit differs from the frozen plan")
            intent = _intent_payload(root, plan_id, plan, reference, baseline)
            # Check every link count and rehash every archive once more after the
            # terminal-gate verification, immediately before publishing intent.
            _preflight_remaining_archives(
                root, plan["immutable_reference"]["archives"], allow_absent=False,
            )
            publish_no_replace(root, intent_path(plan_id), intent)
        else:
            intent = read_json(root, intent_path(plan_id))

        _validate_intent(root, plan_id, plan, intent)

        expected_archives = plan["immutable_reference"]["archives"]
        # This family-wide pass happens before this invocation removes anything,
        # so a later hard-link or identity anomaly cannot cause a partial batch.
        verified_before_any_new_unlink = _preflight_remaining_archives(
            root, expected_archives, allow_absent=True,
        )
        for expected in expected_archives:
            path = safe_path(root, str(expected["path"]), "archive")
            if path.exists() or path.is_symlink():
                _unlink_exact_archive(root, expected)

        reference_after = validate_reference_metadata(root, require_archives=False)
        if any(reference_after["archive_presence"]):
            raise ArchiveEvictionError("one or more acquisition archives remain after eviction")
        if (
            reference_after["extracted_family_sha256"]
            != plan["immutable_reference"]["extracted_family_sha256"]
            or reference_after["extracted_payloads"]
            != plan["immutable_reference"]["extracted_payloads"]
        ):
            raise ArchiveEvictionError("extracted scientific reference changed during archive eviction")
        receipt = {
            "schema_version": RECEIPT_SCHEMA,
            "analysis_id": "track-b-v1.0-local",
            "state": ARCHIVES_EVICTED,
            "plan_id": plan_id,
            "plan": stable_identity(root, plan_path(plan_id)),
            "intent": stable_identity(root, intent_path(plan_id)),
            "implementation": plan["implementation"],
            "terminal_binding": intent["terminal_binding"],
            "archive_absence": [
                {"path": item["path"], "state": "ABSENT_AFTER_FSYNCED_EXACT_UNLINK"}
                for item in expected_archives
            ],
            "archive_count_removed": EXPECTED_ARCHIVE_COUNT,
            "logical_bytes_removed": EXPECTED_ARCHIVE_BYTES,
            "remaining_archives_rehashed_before_this_invocation": verified_before_any_new_unlink,
            "extracted_payloads_rehashed_after_eviction": reference_after["extracted_payloads"],
            "extracted_family_sha256": reference_after["extracted_family_sha256"],
            "reference_provenance": reference_after["reference_provenance_identity"],
            "recovery": plan["immutable_reference"]["recovery"],
            "legacy_verifier_warning": plan["legacy_verifier_warning"],
            "scientific_reference_modified": False,
        }
        publish_no_replace(root, receipt_path(plan_id), receipt)
        return verify_final_receipt(root, plan_id)


def verify_final_receipt(root: Path = ROOT, plan_id: str = "") -> dict[str, Any]:
    plan_id = _require_plan_id(plan_id)
    plan = _load_plan(root, plan_id, require_live_environment=False)
    intent = read_json(root, intent_path(plan_id))
    _validate_intent(root, plan_id, plan, intent)
    receipt = read_json(root, receipt_path(plan_id))
    if (
        set(receipt) != {
            "schema_version", "analysis_id", "state", "plan_id", "plan", "intent",
            "implementation", "terminal_binding", "archive_absence", "archive_count_removed",
            "logical_bytes_removed", "remaining_archives_rehashed_before_this_invocation",
            "extracted_payloads_rehashed_after_eviction", "extracted_family_sha256",
            "reference_provenance", "recovery", "legacy_verifier_warning",
            "scientific_reference_modified",
        }
        or receipt.get("schema_version") != RECEIPT_SCHEMA
        or receipt.get("analysis_id") != "track-b-v1.0-local"
        or receipt.get("state") != ARCHIVES_EVICTED
        or receipt.get("plan_id") != plan_id
        or receipt.get("plan") != stable_identity(root, plan_path(plan_id))
        or receipt.get("intent") != stable_identity(root, intent_path(plan_id))
        or receipt.get("implementation") != plan.get("implementation")
        or receipt.get("terminal_binding") != intent.get("terminal_binding")
        or receipt.get("archive_count_removed") != EXPECTED_ARCHIVE_COUNT
        or receipt.get("logical_bytes_removed") != EXPECTED_ARCHIVE_BYTES
        or receipt.get("scientific_reference_modified") is not False
        or receipt.get("recovery") != plan.get("immutable_reference", {}).get("recovery")
        or receipt.get("legacy_verifier_warning") != plan.get("legacy_verifier_warning")
    ):
        raise ArchiveEvictionError("final archive-eviction receipt schema or lineage drifted")
    reference = validate_reference_metadata(root, require_archives=False)
    if any(reference["archive_presence"]):
        raise ArchiveEvictionError("an archive is present despite the final eviction receipt")
    expected_immutable = {
        key: reference[key]
        for key in (
            "policy", "source_registry", "download_manifest", "extracted_manifest",
            "reference_provenance_identity", "archive_family_sha256",
            "extracted_family_sha256", "archives", "extracted_payloads", "recovery",
        )
    }
    if (
        plan.get("immutable_reference") != expected_immutable
        or receipt.get("extracted_payloads_rehashed_after_eviction")
        != reference["extracted_payloads"]
        or receipt.get("extracted_family_sha256") != reference["extracted_family_sha256"]
        or reference["extracted_family_sha256"]
        != plan["immutable_reference"]["extracted_family_sha256"]
        or receipt.get("reference_provenance") != reference["reference_provenance_identity"]
    ):
        raise ArchiveEvictionError("post-eviction extracted reference differs from final receipt")
    expected_absence = [
        {"path": item["path"], "state": "ABSENT_AFTER_FSYNCED_EXACT_UNLINK"}
        for item in reference["archives"]
    ]
    if receipt.get("archive_absence") != expected_absence:
        raise ArchiveEvictionError("final receipt archive-absence family drifted")
    remaining = receipt.get("remaining_archives_rehashed_before_this_invocation")
    expected_archives = reference["archives"]
    if (
        not isinstance(remaining, list)
        or len({item.get("path") for item in remaining if isinstance(item, dict)}) != len(remaining)
        or any(
            not isinstance(item, dict)
            or item not in expected_archives
            for item in remaining
        )
    ):
        raise ArchiveEvictionError("final receipt pre-unlink archive subset drifted")
    return receipt


def verify_reference_state(root: Path = ROOT) -> dict[str, Any]:
    """Return one of two exact usable reference states, rehashing all payloads.

    This is the narrow importable API for downstream LAVA/PLACO/fine-mapping
    gates.  It intentionally raises on mixed presence, missing receipt, stale
    provenance, a changed payload, or any state other than the two constants.
    """

    reference = validate_reference_metadata(root, require_archives=False)
    present = sum(reference["archive_presence"])
    receipt_ids = _lifecycle_ids(root, RECEIPT_ROOT, "archive-eviction receipt root")
    intent_ids = _lifecycle_ids(root, INTENT_ROOT, "archive-eviction intent root")

    if receipt_ids:
        if len(receipt_ids) != 1 or intent_ids != receipt_ids:
            raise ArchiveEvictionError("ambiguous or incomplete final archive-eviction lineage")
        plan_id = receipt_ids[0]
        verify_final_receipt(root, plan_id)
        return {
            "state": ARCHIVES_EVICTED,
            "reference_provenance": reference["reference_provenance"],
            "reference_provenance_identity": reference["reference_provenance_identity"],
            "archive_count_present": 0,
            "archive_family_sha256": reference["archive_family_sha256"],
            "extracted_payload_count": EXPECTED_EXTRACTED_COUNT,
            "extracted_family_sha256": reference["extracted_family_sha256"],
            "eviction_plan_id": plan_id,
            "eviction_receipt": stable_identity(root, receipt_path(plan_id)),
        }

    if intent_ids:
        if len(intent_ids) != 1:
            raise ArchiveEvictionError("ambiguous durable eviction intents require manual recovery")
        raise ArchiveEvictionError(
            f"{PARTIAL_STATE}: {present}/{EXPECTED_ARCHIVE_COUNT} archives remain; "
            f"resume --evict --plan-id {intent_ids[0]}"
        )

    if present == EXPECTED_ARCHIVE_COUNT:
        # validate_reference_metadata already rehashed all seven plus all 44.
        return {
            "state": ARCHIVES_PRESENT,
            "reference_provenance": reference["reference_provenance"],
            "reference_provenance_identity": reference["reference_provenance_identity"],
            "archive_count_present": EXPECTED_ARCHIVE_COUNT,
            "archive_family_sha256": reference["archive_family_sha256"],
            "extracted_payload_count": EXPECTED_EXTRACTED_COUNT,
            "extracted_family_sha256": reference["extracted_family_sha256"],
            "eviction_plan_id": None,
            "eviction_receipt": None,
        }
    if present == 0:
        raise ArchiveEvictionError(
            "all acquisition archives are absent without a final eviction receipt; manual recovery required"
        )
    raise ArchiveEvictionError(
        "partial archive family exists without a durable eviction intent; manual recovery required"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--audit", action="store_true")
    action.add_argument("--plan", action="store_true")
    action.add_argument("--verify-state", action="store_true")
    action.add_argument("--verify-receipt", action="store_true")
    action.add_argument("--evict", action="store_true")
    parser.add_argument("--plan-id", default="")
    parser.add_argument(
        "--additional-consumer", action="append", default=[], metavar="SCRIPTS/FILE",
        help="bind an additional post-LAVA consumer into --audit/--plan (repeatable)",
    )
    parser.add_argument("--confirm", default="")
    parser.add_argument("--acknowledge-operational-verifier-breakage", action="store_true")
    args = parser.parse_args()
    if args.audit:
        payload = plan_payload(ROOT, args.additional_consumer)
        print(json.dumps({
            "status": "AUDIT_PASS_NO_FILES_REMOVED",
            "archive_count": payload["scope"]["archive_count"],
            "archive_bytes": payload["scope"]["archive_bytes"],
            "payload_count": payload["scope"]["scientific_payload_count"],
            "continuation": payload["lava_dependency_audit"]["current_continuation_fingerprint"],
            "scientific_assessment": payload["scientific_assessment"],
        }, indent=2, sort_keys=True))
    elif args.plan:
        payload = create_or_verify_plan(ROOT, args.additional_consumer)
        plan_id = digest_json(payload)
        print(json.dumps({
            "status": "PLAN_READY_NO_FILES_REMOVED",
            "plan_id": plan_id,
            "path": str(plan_path(plan_id)),
            "sha256": stable_identity(ROOT, plan_path(plan_id))["sha256"],
            "archive_bytes": payload["scope"]["archive_bytes"],
        }, indent=2, sort_keys=True))
    elif args.verify_state:
        print(json.dumps(verify_reference_state(ROOT), indent=2, sort_keys=True))
    elif args.verify_receipt:
        print(json.dumps(verify_final_receipt(ROOT, args.plan_id), indent=2, sort_keys=True))
    else:
        receipt = execute_eviction(
            root=ROOT,
            plan_id=args.plan_id,
            confirmation=args.confirm,
            acknowledge_operational_breakage=args.acknowledge_operational_verifier_breakage,
        )
        print(json.dumps({
            "status": receipt["state"],
            "plan_id": args.plan_id,
            "receipt": str(receipt_path(args.plan_id)),
            "logical_bytes_removed": receipt["logical_bytes_removed"],
        }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ArchiveEvictionError as error:
        raise SystemExit(f"ERROR: {error}") from error
