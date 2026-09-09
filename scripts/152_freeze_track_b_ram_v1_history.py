#!/usr/bin/env python3
"""Freeze the provisional V1 Track B RAM artifacts as immutable history.

The final RAM-aware report must federate several independently validated
execution families.  The current top-level files are the completed V1 LAVA
measurement source, so preserve their exact bytes before a later finalizer
atomically replaces the top-level view.  This utility does not reinterpret or
repair the historical measurements and never writes a scientific result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import stat
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
SCRIPT_RELATIVE = Path("scripts/152_freeze_track_b_ram_v1_history.py")
SOURCE_FILES = {
    "benchmark": Path("results/track_b/RAM_BENCHMARK.tsv"),
    "namespace": Path("results/track_b/RAM_BENCHMARK.namespace.json"),
    "provenance": Path("results/track_b/RAM_BENCHMARK.provenance.json"),
    "report": Path("results/track_b/RAM_AWARE_EXECUTION_REPORT.md"),
}
HISTORY_ROOT = Path("results/track_b/ram_history/lava_v1")
RECEIPT_NAME = "history.receipt.json"
SCHEMA = "sleep-atlas-track-b-ram-v1-history.1"
BENCHMARK_FIELDS = [
    "analysis", "pair", "locus", "chromosome", "n_snps", "peak_ram_gb",
    "runtime_sec", "exit_status", "output_hash",
]
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class HistoryError(RuntimeError):
    """A historical artifact is absent, unsafe, malformed, or drifted."""


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _safe_relative(root: Path, relative: Path, label: str) -> Path:
    if relative.is_absolute() or ".." in relative.parts:
        raise HistoryError(f"{label} is not a safe repository-relative path: {relative}")
    root_real = root.resolve(strict=True)
    current = root
    for part in relative.parts:
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise HistoryError(f"could not inspect {label}: {current}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise HistoryError(f"{label} contains a symbolic link: {current}")
    try:
        (root / relative).resolve(strict=False).relative_to(root_real)
    except (OSError, RuntimeError, ValueError) as error:
        raise HistoryError(f"{label} escapes the repository: {relative}") from error
    return root / relative


def stable_bytes(path: Path, label: str) -> tuple[bytes, dict[str, Any]]:
    try:
        first = os.lstat(path)
        if not stat.S_ISREG(first.st_mode) or first.st_size <= 0:
            raise HistoryError(f"{label} is not a nonempty regular file: {path}")
        if first.st_nlink != 1:
            raise HistoryError(f"{label} has an unexpected hard-link count: {path}")
        with path.open("rb") as handle:
            opened = os.fstat(handle.fileno())
            content = handle.read()
            after_read = os.fstat(handle.fileno())
        final = os.lstat(path)
    except OSError as error:
        raise HistoryError(f"could not read {label}: {path}: {error}") from error
    identities = {
        (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns)
        for item in (first, opened, after_read, final)
    }
    if len(identities) != 1 or len(content) != final.st_size:
        raise HistoryError(f"{label} changed while it was read: {path}")
    return content, {
        "bytes": len(content),
        "device": int(final.st_dev),
        "inode": int(final.st_ino),
        "mtime_ns": int(final.st_mtime_ns),
        "ctime_ns": int(final.st_ctime_ns),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def portable_identity(relative: Path, content: bytes) -> dict[str, Any]:
    return {
        "path": str(relative),
        "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }


def _json_object(content: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(content.decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise HistoryError(f"{label} is not valid UTF-8 JSON: {error}") from error
    if not isinstance(value, dict):
        raise HistoryError(f"{label} must contain one JSON object")
    return value


def validate_benchmark(content: bytes) -> int:
    try:
        text = content.decode("utf-8")
    except UnicodeError as error:
        raise HistoryError(f"historical benchmark is not UTF-8: {error}") from error
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if list(reader.fieldnames or []) != BENCHMARK_FIELDS:
        raise HistoryError("historical benchmark does not have the exact nine-column schema")
    rows = 0
    for row_number, row in enumerate(reader, start=2):
        rows += 1
        if None in row or any(row[field] == "" for field in BENCHMARK_FIELDS):
            raise HistoryError(f"malformed historical benchmark row {row_number}")
        if not row["analysis"].startswith("LAVA_"):
            raise HistoryError(f"non-LAVA row in V1 historical benchmark at row {row_number}")
        if not SHA256_RE.fullmatch(row["output_hash"]):
            raise HistoryError(f"invalid output hash at historical benchmark row {row_number}")
        try:
            peak = float(row["peak_ram_gb"])
            runtime = float(row["runtime_sec"])
            status = int(row["exit_status"].split(":", 1)[0])
        except (ValueError, OverflowError) as error:
            raise HistoryError(f"invalid numeric field at benchmark row {row_number}") from error
        if not math.isfinite(peak) or peak < 0 or not math.isfinite(runtime) or runtime < 0:
            raise HistoryError(f"invalid RAM/runtime at historical benchmark row {row_number}")
        if status < 0 or status > 255:
            raise HistoryError(f"invalid exit status at historical benchmark row {row_number}")
        if ":" in row["exit_status"]:
            label = row["exit_status"].split(":", 1)[1]
            if not label or not re.fullmatch(r"[A-Z][A-Z0-9_]*", label):
                raise HistoryError(f"invalid exit-status label at benchmark row {row_number}")
        if row["locus"] != "ALL":
            try:
                locus = int(row["locus"])
                chromosome = int(row["chromosome"])
            except ValueError as error:
                raise HistoryError(f"invalid locus identity at benchmark row {row_number}") from error
            if row["n_snps"] == "NA":
                n_snps = None
            else:
                try:
                    n_snps = int(row["n_snps"])
                except ValueError as error:
                    raise HistoryError(
                        f"invalid SNP count at benchmark row {row_number}"
                    ) from error
            if (
                locus < 1
                or chromosome not in range(1, 23)
                or (n_snps is not None and n_snps < 0)
            ):
                raise HistoryError(f"out-of-range locus identity at benchmark row {row_number}")
    if rows < 2495:
        raise HistoryError(f"historical benchmark is incomplete: observed {rows} rows")
    return rows


def provenance_audit(
    benchmark: bytes, row_count: int, provenance: Mapping[str, Any],
) -> dict[str, Any]:
    attempts = provenance.get("all_lava_attempt_rows")
    if not isinstance(attempts, list):
        raise HistoryError("historical provenance lacks its attempted-row family")
    declared_hash = provenance.get("benchmark_sha256")
    declared_bytes = provenance.get("benchmark_bytes")
    if not isinstance(declared_hash, str) or not SHA256_RE.fullmatch(declared_hash):
        raise HistoryError("historical provenance has an invalid benchmark hash")
    if not isinstance(declared_bytes, int) or declared_bytes <= 0:
        raise HistoryError("historical provenance has invalid benchmark byte accounting")
    observed_hash = hashlib.sha256(benchmark).hexdigest()
    observed_bytes = len(benchmark)
    declared_rows = len(attempts)
    return {
        "fully_consistent": (
            declared_hash == observed_hash
            and declared_bytes == observed_bytes
            and declared_rows == row_count
        ),
        "observed_benchmark_bytes": observed_bytes,
        "declared_benchmark_bytes": declared_bytes,
        "benchmark_bytes_match": declared_bytes == observed_bytes,
        "observed_benchmark_sha256": observed_hash,
        "declared_benchmark_sha256": declared_hash,
        "benchmark_sha256_match": declared_hash == observed_hash,
        "observed_benchmark_rows": row_count,
        "declared_attempt_rows": declared_rows,
        "row_count_match": declared_rows == row_count,
        "handling": (
            "PRESERVED_AS_PROVISIONAL_HISTORY_ONLY;FINAL_FEDERATION_MUST_REVALIDATE_RAW_RECEIPTS"
        ),
    }


def load_source(root: Path) -> dict[str, Any]:
    contents: dict[str, bytes] = {}
    live: dict[str, dict[str, Any]] = {}
    for role, relative in SOURCE_FILES.items():
        path = _safe_relative(root, relative, f"historical {role} source")
        content, identity = stable_bytes(path, f"historical {role} source")
        contents[role] = content
        live[role] = {"path": str(relative), **identity}
    row_count = validate_benchmark(contents["benchmark"])
    namespace = _json_object(contents["namespace"], "historical namespace")
    provenance = _json_object(contents["provenance"], "historical provenance")
    fingerprint = namespace.get("execution_fingerprint")
    if (
        namespace.get("schema_version") != 1
        or namespace.get("analysis_id") != "track-b-v1.0-local"
        or not isinstance(fingerprint, str)
        or not SHA256_RE.fullmatch(fingerprint)
    ):
        raise HistoryError("historical namespace identity drifted")
    if (
        provenance.get("schema_version") != 2
        or provenance.get("analysis_id") != "track-b-v1.0-local"
        or provenance.get("execution_fingerprint") != fingerprint
        or provenance.get("benchmark_path") != str(SOURCE_FILES["benchmark"])
    ):
        raise HistoryError("historical benchmark provenance identity drifted")
    audit = provenance_audit(contents["benchmark"], row_count, provenance)
    if b"Track B RAM-aware execution report" not in contents["report"]:
        raise HistoryError("historical RAM report identity drifted")
    return {
        "fingerprint": fingerprint,
        "row_count": row_count,
        "contents": contents,
        "live": live,
        "provenance_audit": audit,
    }


def history_directory(root: Path, fingerprint: str) -> Path:
    return _safe_relative(root, HISTORY_ROOT / fingerprint, "historical destination")


def expected_receipt(root: Path, source: Mapping[str, Any]) -> dict[str, Any]:
    fingerprint = str(source["fingerprint"])
    directory_relative = HISTORY_ROOT / fingerprint
    snapshots = {
        role: portable_identity(directory_relative / relative.name, source["contents"][role])
        for role, relative in SOURCE_FILES.items()
    }
    script_content, _ = stable_bytes(SCRIPT, "history-freezer script")
    return {
        "schema_version": SCHEMA,
        "analysis_id": "track-b-v1.0-local",
        "artifact_role": "IMMUTABLE_PROVISIONAL_V1_RAM_HISTORY_NOT_NEW_SCIENCE",
        "execution_fingerprint": fingerprint,
        "benchmark_row_count": int(source["row_count"]),
        "exact_benchmark_schema": BENCHMARK_FIELDS,
        "source_identities_at_freeze": {
            role: {
                key: value for key, value in identity.items()
                if key in {"path", "bytes", "sha256"}
            }
            for role, identity in source["live"].items()
        },
        "historical_provenance_audit": source["provenance_audit"],
        "snapshots": snapshots,
        "generator": {
            "path": str(SCRIPT_RELATIVE),
            "bytes": len(script_content),
            "sha256": hashlib.sha256(script_content).hexdigest(),
        },
        "future_top_level_replacement_allowed_only_after_snapshot_verifies": True,
    }


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def publish_no_replace(path: Path, content: bytes) -> None:
    if path.exists() or path.is_symlink():
        observed, _ = stable_bytes(path, "existing historical snapshot")
        if observed != content:
            raise HistoryError(f"existing historical snapshot differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    try:
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
    _fsync_directory(path.parent)


def verify_snapshot(root: Path, fingerprint: str) -> dict[str, Any]:
    directory = history_directory(root, fingerprint)
    receipt_path = directory / RECEIPT_NAME
    receipt_content, _ = stable_bytes(receipt_path, "historical receipt")
    receipt = _json_object(receipt_content, "historical receipt")
    if (
        receipt.get("schema_version") != SCHEMA
        or receipt.get("analysis_id") != "track-b-v1.0-local"
        or receipt.get("artifact_role")
        != "IMMUTABLE_PROVISIONAL_V1_RAM_HISTORY_NOT_NEW_SCIENCE"
        or receipt.get("execution_fingerprint") != fingerprint
        or receipt.get("exact_benchmark_schema") != BENCHMARK_FIELDS
        or receipt.get("future_top_level_replacement_allowed_only_after_snapshot_verifies") is not True
    ):
        raise HistoryError("historical receipt identity drifted")
    snapshots = receipt.get("snapshots")
    sources = receipt.get("source_identities_at_freeze")
    if not isinstance(snapshots, dict) or set(snapshots) != set(SOURCE_FILES):
        raise HistoryError("historical receipt snapshot family is incomplete")
    if not isinstance(sources, dict) or set(sources) != set(SOURCE_FILES):
        raise HistoryError("historical receipt source family is incomplete")
    observed_contents: dict[str, bytes] = {}
    for role, record in snapshots.items():
        if not isinstance(record, dict):
            raise HistoryError(f"malformed historical snapshot record: {role}")
        relative = Path(str(record.get("path", "")))
        expected_relative = HISTORY_ROOT / fingerprint / SOURCE_FILES[role].name
        if relative != expected_relative:
            raise HistoryError(f"historical snapshot path drifted: {role}")
        content, identity = stable_bytes(
            _safe_relative(root, relative, f"historical {role} snapshot"),
            f"historical {role} snapshot",
        )
        if identity["bytes"] != record.get("bytes") or identity["sha256"] != record.get("sha256"):
            raise HistoryError(f"historical snapshot identity drifted: {role}")
        source_record = sources[role]
        if (
            not isinstance(source_record, dict)
            or source_record.get("path") != str(SOURCE_FILES[role])
            or source_record.get("bytes") != record.get("bytes")
            or source_record.get("sha256") != record.get("sha256")
        ):
            raise HistoryError(f"historical source/snapshot binding drifted: {role}")
        observed_contents[role] = content
    generator = receipt.get("generator")
    script_content, _ = stable_bytes(SCRIPT, "history-freezer script")
    if (
        not isinstance(generator, dict)
        or generator.get("path") != str(SCRIPT_RELATIVE)
        or generator.get("bytes") != len(script_content)
        or generator.get("sha256") != hashlib.sha256(script_content).hexdigest()
    ):
        raise HistoryError("historical receipt generator identity drifted")
    row_count = validate_benchmark(observed_contents["benchmark"])
    provenance = _json_object(observed_contents["provenance"], "snapshot provenance")
    namespace = _json_object(observed_contents["namespace"], "snapshot namespace")
    audit = provenance_audit(observed_contents["benchmark"], row_count, provenance)
    if (
        row_count != receipt.get("benchmark_row_count")
        or namespace.get("execution_fingerprint") != fingerprint
        or provenance.get("execution_fingerprint") != fingerprint
        or receipt.get("historical_provenance_audit") != audit
    ):
        raise HistoryError("historical snapshot family is internally inconsistent")
    return receipt


def freeze(root: Path = ROOT) -> dict[str, Any]:
    source = load_source(root)
    fingerprint = str(source["fingerprint"])
    directory = history_directory(root, fingerprint)
    directory.mkdir(parents=True, exist_ok=True)
    if directory.is_symlink():
        raise HistoryError("historical destination is a symbolic link")
    receipt = expected_receipt(root, source)
    for role, relative in SOURCE_FILES.items():
        publish_no_replace(directory / relative.name, source["contents"][role])
    # Ensure the live sources did not drift between initial hashing and receipt publication.
    for role, relative in SOURCE_FILES.items():
        _, identity = stable_bytes(
            _safe_relative(root, relative, f"historical {role} source"),
            f"historical {role} source",
        )
        expected = receipt["source_identities_at_freeze"][role]
        if identity["bytes"] != expected["bytes"] or identity["sha256"] != expected["sha256"]:
            raise HistoryError(f"historical live source drifted before receipt publication: {role}")
    receipt_bytes = json.dumps(receipt, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    publish_no_replace(directory / RECEIPT_NAME, receipt_bytes)
    _fsync_directory(directory)
    return verify_snapshot(root, fingerprint)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--freeze", action="store_true")
    mode.add_argument("--verify", action="store_true")
    parser.add_argument("--fingerprint")
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        source = load_source(ROOT)
        print(json.dumps({
            "status": "PREFLIGHT_PASS",
            "execution_fingerprint": source["fingerprint"],
            "benchmark_row_count": source["row_count"],
            "historical_provenance_fully_consistent": source["provenance_audit"][
                "fully_consistent"
            ],
            "writes_performed": False,
        }, indent=2, sort_keys=True))
        return
    if args.freeze:
        if not args.execute:
            raise HistoryError("explicit --execute is required to freeze V1 RAM history")
        receipt = freeze(ROOT)
        print(
            "TRACK_B_RAM_V1_HISTORY_FROZEN "
            f"fingerprint={receipt['execution_fingerprint']} "
            f"rows={receipt['benchmark_row_count']}"
        )
        return
    if not args.fingerprint or not SHA256_RE.fullmatch(args.fingerprint):
        raise HistoryError("--verify requires the exact 64-character fingerprint")
    receipt = verify_snapshot(ROOT, args.fingerprint)
    print(
        "TRACK_B_RAM_V1_HISTORY_VERIFIED "
        f"fingerprint={receipt['execution_fingerprint']} "
        f"rows={receipt['benchmark_row_count']}"
    )


if __name__ == "__main__":
    try:
        main()
    except (HistoryError, OSError) as error:
        raise SystemExit(f"ERROR: {error}") from error
