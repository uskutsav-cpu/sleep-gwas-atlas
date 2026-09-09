#!/usr/bin/env python3
"""Materialize and verify exact chromosome-scoped Track B LAVA inputs.

The source summary statistics remain the frozen genome-wide files.  A shard is
an order-preserving byte-for-byte projection of one source file onto the
normalized SNP identifiers present in one sealed LAVA chromosome reference.
No association-statistic, allele, sample-size, or result-dependent filter is
applied here.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import stat
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO


DEFAULT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = Path(__file__).resolve()
CHROMOSOMES = tuple(range(1, 23))
INPUT_INFO_FIELDS = ["phenotype", "cases", "controls", "prevalence", "filename"]
INPUT_PROVENANCE_FIELDS = [
    "phenotype", "filename", "bytes", "sha256", "observed_columns", "validation_status",
]
EXPECTED_SUMSTAT_FIELDS = ["SNP", "A1", "A2", "Z", "N"]
SCHEMA_VERSION = "track-b-lava-chromosome-inputs.1"
PHENOTYPE_RE = re.compile(r"^[a-z][a-z0-9_]*$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class ShardError(RuntimeError):
    """Fail-closed materialization or verification error."""


@dataclass(frozen=True)
class Snapshot:
    path: Path
    relative: str
    bytes: int
    sha256: str
    stat_key: tuple[int, int, int, int, int]

    def record(self) -> dict[str, object]:
        return {"path": self.relative, "bytes": self.bytes, "sha256": self.sha256}


@dataclass(frozen=True)
class IntersectionMetrics:
    header_sha256: str
    data_row_count: int
    unique_snp_count: int
    duplicate_snp_row_count: int
    decompressed_sha256: str
    data_rows_sha256: str
    normalized_snp_sequence_sha256: str

    def record(self) -> dict[str, object]:
        return {
            "header_sha256": self.header_sha256,
            "data_row_count": self.data_row_count,
            "unique_snp_count": self.unique_snp_count,
            "duplicate_snp_row_count": self.duplicate_snp_row_count,
            "decompressed_sha256": self.decompressed_sha256,
            "data_rows_sha256": self.data_rows_sha256,
            "normalized_snp_sequence_sha256": self.normalized_snp_sequence_sha256,
        }


@dataclass
class Context:
    root: Path
    policy: dict[str, object]
    input_lock: dict[str, object]
    input_rows: list[dict[str, str]]
    source_snapshots: dict[str, Snapshot]
    source_columns: dict[str, list[str]]
    policy_snapshot: Snapshot
    input_lock_snapshot: Snapshot
    input_info_snapshot: Snapshot
    input_provenance_snapshot: Snapshot
    reference_provenance_snapshot: Snapshot
    reference_records: dict[str, dict[str, object]]
    script_snapshot: Snapshot
    output_root: Path
    output_lock: Path


def stat_key(value: os.stat_result) -> tuple[int, int, int, int, int]:
    return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def safe_path(root: Path, value: str, label: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ShardError(f"unsafe {label} path: {value}")
    candidate = root / relative
    try:
        candidate.resolve(strict=False).relative_to(root.resolve())
    except ValueError as error:
        raise ShardError(f"{label} path escapes the repository root: {value}") from error
    return candidate


def relative_path(root: Path, path: Path) -> str:
    try:
        return path.resolve(strict=False).relative_to(root.resolve()).as_posix()
    except ValueError as error:
        raise ShardError(f"artifact is outside the repository root: {path}") from error


def snapshot(path: Path, root: Path, *, allow_outside_root: bool = False) -> Snapshot:
    if path.is_symlink():
        raise ShardError(f"symbolic-link artifacts are forbidden: {path}")
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            if not stat.S_ISREG(before.st_mode) or before.st_size <= 0:
                raise ShardError(f"missing real non-empty artifact: {path}")
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            after = os.fstat(handle.fileno())
        current = path.stat()
    except OSError as error:
        raise ShardError(f"could not hash artifact {path}: {error}") from error
    if stat_key(before) != stat_key(after) or stat_key(after) != stat_key(current):
        raise ShardError(f"artifact changed while hashing: {path}")
    relative = str(path) if allow_outside_root else relative_path(root, path)
    return Snapshot(path, relative, after.st_size, digest.hexdigest(), stat_key(after))


def stable_bytes(item: Snapshot) -> bytes:
    try:
        with item.path.open("rb") as handle:
            before = os.fstat(handle.fileno())
            value = handle.read()
            after = os.fstat(handle.fileno())
        current = item.path.stat()
    except OSError as error:
        raise ShardError(f"could not read artifact {item.path}: {error}") from error
    if stat_key(before) != item.stat_key or stat_key(after) != item.stat_key or stat_key(current) != item.stat_key:
        raise ShardError(f"artifact changed after it was frozen for this run: {item.path}")
    if hashlib.sha256(value).hexdigest() != item.sha256:
        raise ShardError(f"artifact content differs from its frozen identity: {item.path}")
    return value


def load_json(item: Snapshot, label: str) -> dict[str, object]:
    try:
        value = json.loads(stable_bytes(item).decode("utf-8"))
    except (UnicodeError, json.JSONDecodeError) as error:
        raise ShardError(f"{label} is not valid UTF-8 JSON: {error}") from error
    if not isinstance(value, dict):
        raise ShardError(f"{label} must be a JSON object")
    return value


def read_tsv(item: Snapshot, expected_fields: list[str], label: str) -> list[dict[str, str]]:
    try:
        text = stable_bytes(item).decode("utf-8")
    except UnicodeError as error:
        raise ShardError(f"{label} is not valid UTF-8") from error
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    fields = list(reader.fieldnames or [])
    if fields != expected_fields:
        raise ShardError(f"{label} schema drifted: observed={fields} expected={expected_fields}")
    rows = list(reader)
    if any(None in row or None in row.values() for row in rows):
        raise ShardError(f"{label} contains a malformed row")
    return rows


def validate_sha256(value: object, label: str) -> str:
    observed = str(value)
    if not SHA256_RE.fullmatch(observed):
        raise ShardError(f"invalid SHA-256 for {label}: {observed}")
    return observed


def validate_int(value: object, label: str, *, positive: bool = False) -> int:
    try:
        observed = int(value)
    except (TypeError, ValueError) as error:
        raise ShardError(f"invalid integer for {label}: {value}") from error
    if observed < (1 if positive else 0):
        raise ShardError(f"out-of-range integer for {label}: {observed}")
    return observed


def normalize_snp(raw: bytes, label: str) -> str:
    try:
        value = raw.decode("utf-8").strip().lower()
    except UnicodeError as error:
        raise ShardError(f"non-UTF-8 SNP identifier in {label}") from error
    if not value or any(character.isspace() for character in value):
        raise ShardError(f"empty or whitespace-containing SNP identifier in {label}: {value!r}")
    return value


def decode_fields(fields: list[bytes], label: str) -> list[str]:
    try:
        return [field.decode("utf-8") for field in fields]
    except UnicodeError as error:
        raise ShardError(f"non-UTF-8 header in {label}") from error


def split_line(line: bytes, expected_columns: int, label: str) -> list[bytes]:
    if not line:
        raise ShardError(f"unexpected empty line in {label}")
    body = line[:-2] if line.endswith(b"\r\n") else line[:-1] if line.endswith(b"\n") else line
    if not body:
        raise ShardError(f"blank row in {label}")
    fields = body.split(b"\t")
    if len(fields) != expected_columns:
        raise ShardError(
            f"wrong column count in {label}: observed={len(fields)} expected={expected_columns}"
        )
    return fields


def assert_frozen_stat(item: Snapshot, observed: os.stat_result, label: str) -> None:
    if stat_key(observed) != item.stat_key:
        raise ShardError(f"{label} changed after it was frozen for this run: {item.path}")


def inspect_gzip_header(item: Snapshot, expected_columns: list[str]) -> bytes:
    try:
        with item.path.open("rb") as raw:
            assert_frozen_stat(item, os.fstat(raw.fileno()), "source summary statistics")
            with gzip.GzipFile(fileobj=raw, mode="rb") as stream:
                header = stream.readline()
            after = os.fstat(raw.fileno())
        current = item.path.stat()
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise ShardError(f"could not read gzip header {item.path}: {error}") from error
    assert_frozen_stat(item, after, "source summary statistics")
    assert_frozen_stat(item, current, "source summary statistics")
    fields = decode_fields(
        split_line(header, len(expected_columns), str(item.path)), str(item.path),
    )
    if fields != expected_columns:
        raise ShardError(f"summary-statistics header drifted for {item.path}: {fields}")
    return header


def load_context(root: Path) -> Context:
    root = root.resolve()
    policy_path = root / "config/track_b_local_analysis_policy.json"
    input_lock_path = root / "results/track_b/local_analysis_input.lock.json"
    input_info_path = root / "results/track_b/lava_input_info.tsv"
    input_provenance_path = root / "results/track_b/lava_input_provenance.tsv"
    policy_snapshot = snapshot(policy_path, root)
    input_lock_snapshot = snapshot(input_lock_path, root)
    input_info_snapshot = snapshot(input_info_path, root)
    input_provenance_snapshot = snapshot(input_provenance_path, root)
    policy = load_json(policy_snapshot, "Track B local-analysis policy")
    input_lock = load_json(input_lock_snapshot, "Track B local-analysis input lock")

    trait_order = policy.get("trait_order")
    if (
        policy.get("analysis_id") != "track-b-v1.0-local"
        or policy.get("expected_analysis_traits") != 8
        or not isinstance(trait_order, list)
        or len(trait_order) != 8
        or len(set(map(str, trait_order))) != 8
    ):
        raise ShardError("Track B local-analysis policy is not the exact frozen eight-trait family")
    trait_order = [str(value) for value in trait_order]
    if any(not PHENOTYPE_RE.fullmatch(value) for value in trait_order):
        raise ShardError("unsafe or malformed Track B phenotype identifier")
    if (
        input_lock.get("analysis_id") != policy["analysis_id"]
        or input_lock.get("selection_timing") != "BEFORE_LOCAL_RESULT_ACCESS"
        or input_lock.get("local_results_accessed_before_input_freeze") is not False
        or input_lock.get("policy_sha256") != policy_snapshot.sha256
    ):
        raise ShardError("Track B local input lock is stale or post-result")
    locked_artifacts = input_lock.get("input_artifact_sha256")
    if not isinstance(locked_artifacts, dict):
        raise ShardError("Track B local input lock has no artifact hash map")
    for item in (input_info_snapshot, input_provenance_snapshot):
        if locked_artifacts.get(item.relative) != item.sha256:
            raise ShardError(f"frozen Track B input differs from its lock: {item.relative}")

    input_rows = read_tsv(input_info_snapshot, INPUT_INFO_FIELDS, "Track B LAVA input info")
    provenance_rows = read_tsv(
        input_provenance_snapshot, INPUT_PROVENANCE_FIELDS, "Track B LAVA input provenance",
    )
    if [row["phenotype"] for row in input_rows] != trait_order:
        raise ShardError("Track B LAVA input-info trait order differs from policy")
    if [row["phenotype"] for row in provenance_rows] != trait_order:
        raise ShardError("Track B LAVA input-provenance trait order differs from policy")
    provenance_by_trait = {row["phenotype"]: row for row in provenance_rows}
    source_snapshots: dict[str, Snapshot] = {}
    source_columns: dict[str, list[str]] = {}
    for row in input_rows:
        phenotype = row["phenotype"]
        provenance = provenance_by_trait[phenotype]
        if row["filename"] != provenance["filename"]:
            raise ShardError(f"input-info/provenance filename disagreement for {phenotype}")
        source_path = safe_path(root, row["filename"], f"{phenotype} summary statistics")
        source = snapshot(source_path, root)
        if (
            source.bytes != validate_int(provenance["bytes"], f"{phenotype} source bytes", positive=True)
            or source.sha256 != validate_sha256(provenance["sha256"], f"{phenotype} source")
        ):
            raise ShardError(f"live Track B source differs from frozen provenance: {phenotype}")
        columns = provenance["observed_columns"].split(",")
        if columns != EXPECTED_SUMSTAT_FIELDS:
            raise ShardError(f"unsupported frozen LAVA source schema for {phenotype}: {columns}")
        if provenance["validation_status"] != "VALIDATED_FOR_LAVA_LOCAL_NOT_FINE_MAPPING":
            raise ShardError(f"source is not frozen for LAVA local analysis: {phenotype}")
        inspect_gzip_header(source, columns)
        source_snapshots[phenotype] = source
        source_columns[phenotype] = columns

    reference_relative = input_lock.get("required_reference_provenance")
    if not isinstance(reference_relative, str):
        raise ShardError("Track B input lock does not identify reference provenance")
    reference_provenance_path = safe_path(root, reference_relative, "LAVA reference provenance")
    reference_provenance_snapshot = snapshot(reference_provenance_path, root)
    reference_provenance = load_json(reference_provenance_snapshot, "sealed LAVA reference provenance")
    extracted = reference_provenance.get("extracted_files")
    if (
        reference_provenance.get("schema_version") != "sleep-atlas-lava-reference.1"
        or reference_provenance.get("verification") != "SHA-256 verified after official HTTPS acquisition"
        or not isinstance(extracted, list)
        or len(extracted) < len(CHROMOSOMES)
        or reference_provenance.get("extracted_file_count") != len(extracted)
    ):
        raise ShardError("sealed LAVA reference provenance has no extracted-file family")
    reference_records: dict[str, dict[str, object]] = {}
    for record in extracted:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise ShardError("malformed extracted-file record in LAVA reference provenance")
        relative = str(record["path"])
        if relative in reference_records:
            raise ShardError(f"duplicate LAVA reference provenance path: {relative}")
        reference_records[relative] = record

    reference_prefix = policy.get("reference_prefix")
    if not isinstance(reference_prefix, str):
        raise ShardError("Track B policy does not define a reference prefix")
    for chromosome in CHROMOSOMES:
        info_path = safe_path(
            root, f"{reference_prefix}_chr{chromosome}.info", f"chromosome {chromosome} reference info",
        )
        relative = relative_path(root, info_path)
        record = reference_records.get(relative)
        if record is None:
            raise ShardError(f"sealed reference omits chromosome {chromosome} info: {relative}")
        observed = snapshot(info_path, root)
        if (
            observed.bytes != validate_int(record.get("bytes"), f"reference bytes {relative}", positive=True)
            or observed.sha256 != validate_sha256(record.get("sha256"), f"reference {relative}")
        ):
            raise ShardError(f"live chromosome reference differs from sealed provenance: {relative}")

    script_snapshot = snapshot(SCRIPT, root, allow_outside_root=not SCRIPT.is_relative_to(root))
    return Context(
        root=root,
        policy=policy,
        input_lock=input_lock,
        input_rows=input_rows,
        source_snapshots=source_snapshots,
        source_columns=source_columns,
        policy_snapshot=policy_snapshot,
        input_lock_snapshot=input_lock_snapshot,
        input_info_snapshot=input_info_snapshot,
        input_provenance_snapshot=input_provenance_snapshot,
        reference_provenance_snapshot=reference_provenance_snapshot,
        reference_records=reference_records,
        script_snapshot=script_snapshot,
        output_root=root / "results/track_b/lava_chromosome_inputs",
        output_lock=root / "results/track_b/lava_chromosome_inputs.provenance.json",
    )


def reference_snapshot(context: Context, chromosome: int) -> Snapshot:
    prefix = str(context.policy["reference_prefix"])
    observed = snapshot(
        safe_path(context.root, f"{prefix}_chr{chromosome}.info", "LAVA reference info"),
        context.root,
    )
    sealed = context.reference_records.get(observed.relative)
    if sealed is None or (
        observed.bytes
        != validate_int(sealed.get("bytes"), f"reference bytes {observed.relative}", positive=True)
        or observed.sha256
        != validate_sha256(sealed.get("sha256"), f"reference {observed.relative}")
    ):
        raise ShardError(
            f"live chromosome reference differs from sealed provenance: {observed.relative}"
        )
    return observed


def load_reference_snps(item: Snapshot, chromosome: int) -> tuple[set[str], dict[str, object]]:
    normalized_digest = hashlib.sha256()
    identifiers: set[str] = set()
    rows = 0
    try:
        with item.path.open("rb") as handle:
            assert_frozen_stat(item, os.fstat(handle.fileno()), "reference info")
            header = handle.readline()
            header_fields = split_line(header, len(header.rstrip(b"\r\n").split(b"\t")), str(item.path))
            names = decode_fields(header_fields, str(item.path))
            if "SNP" not in names or "CHR" not in names:
                raise ShardError(f"reference info lacks SNP/CHR columns: {item.path}")
            snp_index, chromosome_index = names.index("SNP"), names.index("CHR")
            for line_number, line in enumerate(handle, start=2):
                fields = split_line(line, len(names), f"{item.path}:{line_number}")
                try:
                    observed_chromosome = int(fields[chromosome_index])
                except ValueError as error:
                    raise ShardError(f"invalid reference chromosome at {item.path}:{line_number}") from error
                if observed_chromosome != chromosome:
                    raise ShardError(
                        f"wrong chromosome in reference info {item.path}:{line_number}: {observed_chromosome}"
                    )
                identifier = normalize_snp(fields[snp_index], f"{item.path}:{line_number}")
                if identifier in identifiers:
                    raise ShardError(f"duplicate normalized SNP in chromosome reference: chr{chromosome}/{identifier}")
                identifiers.add(identifier)
                normalized_digest.update(identifier.encode("utf-8") + b"\n")
                rows += 1
            after = os.fstat(handle.fileno())
        current = item.path.stat()
    except OSError as error:
        raise ShardError(f"could not read reference info {item.path}: {error}") from error
    assert_frozen_stat(item, after, "reference info")
    assert_frozen_stat(item, current, "reference info")
    if not identifiers:
        raise ShardError(f"chromosome {chromosome} reference SNP set is empty")
    record = item.record() | {
        "chromosome": chromosome,
        "row_count": rows,
        "unique_snp_count": len(identifiers),
        "normalized_snp_sequence_sha256": normalized_digest.hexdigest(),
    }
    return identifiers, record


def scan_source_intersection(
    source: Snapshot,
    expected_columns: list[str],
    reference_snps: set[str],
    *,
    output: BinaryIO | None = None,
) -> tuple[IntersectionMetrics, bytes]:
    decompressed = hashlib.sha256()
    data_rows = hashlib.sha256()
    snp_sequence = hashlib.sha256()
    identifiers: set[str] = set()
    count = 0
    try:
        with source.path.open("rb") as raw:
            assert_frozen_stat(source, os.fstat(raw.fileno()), "source summary statistics")
            with gzip.GzipFile(fileobj=raw, mode="rb") as stream:
                header = stream.readline()
                fields = decode_fields(
                    split_line(header, len(expected_columns), str(source.path)), str(source.path),
                )
                if fields != expected_columns:
                    raise ShardError(f"source header drifted while scanning: {source.path}")
                decompressed.update(header)
                if output is not None:
                    output.write(header)
                for line_number, line in enumerate(stream, start=2):
                    row = split_line(line, len(expected_columns), f"{source.path}:{line_number}")
                    identifier = normalize_snp(row[0], f"{source.path}:{line_number}")
                    if identifier not in reference_snps:
                        continue
                    if output is not None:
                        output.write(line)
                    decompressed.update(line)
                    data_rows.update(line)
                    snp_sequence.update(identifier.encode("utf-8") + b"\n")
                    identifiers.add(identifier)
                    count += 1
            after = os.fstat(raw.fileno())
        current = source.path.stat()
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise ShardError(f"could not scan source summary statistics {source.path}: {error}") from error
    assert_frozen_stat(source, after, "source summary statistics")
    assert_frozen_stat(source, current, "source summary statistics")
    metrics = IntersectionMetrics(
        header_sha256=hashlib.sha256(header).hexdigest(),
        data_row_count=count,
        unique_snp_count=len(identifiers),
        duplicate_snp_row_count=count - len(identifiers),
        decompressed_sha256=decompressed.hexdigest(),
        data_rows_sha256=data_rows.hexdigest(),
        normalized_snp_sequence_sha256=snp_sequence.hexdigest(),
    )
    return metrics, header


def scan_shard(
    item: Snapshot,
    expected_columns: list[str],
    expected_header: bytes,
    reference_snps: set[str],
) -> IntersectionMetrics:
    decompressed = hashlib.sha256()
    data_rows = hashlib.sha256()
    snp_sequence = hashlib.sha256()
    identifiers: set[str] = set()
    count = 0
    try:
        with item.path.open("rb") as raw:
            assert_frozen_stat(item, os.fstat(raw.fileno()), "chromosome shard")
            with gzip.GzipFile(fileobj=raw, mode="rb") as stream:
                header = stream.readline()
                if header != expected_header:
                    raise ShardError(f"shard header is not byte-identical to its source: {item.relative}")
                fields = decode_fields(
                    split_line(header, len(expected_columns), item.relative), item.relative,
                )
                if fields != expected_columns:
                    raise ShardError(f"shard schema drifted: {item.relative}")
                decompressed.update(header)
                for line_number, line in enumerate(stream, start=2):
                    row = split_line(line, len(expected_columns), f"{item.relative}:{line_number}")
                    identifier = normalize_snp(row[0], f"{item.relative}:{line_number}")
                    if identifier not in reference_snps:
                        raise ShardError(f"shard contains a SNP outside its chromosome reference: {item.relative}")
                    decompressed.update(line)
                    data_rows.update(line)
                    snp_sequence.update(identifier.encode("utf-8") + b"\n")
                    identifiers.add(identifier)
                    count += 1
            after = os.fstat(raw.fileno())
        current = item.path.stat()
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise ShardError(f"could not scan chromosome shard {item.path}: {error}") from error
    assert_frozen_stat(item, after, "chromosome shard")
    assert_frozen_stat(item, current, "chromosome shard")
    return IntersectionMetrics(
        header_sha256=hashlib.sha256(header).hexdigest(),
        data_row_count=count,
        unique_snp_count=len(identifiers),
        duplicate_snp_row_count=count - len(identifiers),
        decompressed_sha256=decompressed.hexdigest(),
        data_rows_sha256=data_rows.hexdigest(),
        normalized_snp_sequence_sha256=snp_sequence.hexdigest(),
    )


def fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def link_no_replace(temporary: Path, destination: Path) -> None:
    try:
        os.link(temporary, destination)
    except FileExistsError:
        return
    except OSError as error:
        raise ShardError(f"could not publish artifact without replacement: {destination}: {error}") from error
    fsync_directory(destination.parent)


def publish_bytes_no_replace(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        try:
            observed = path.read_bytes()
        except OSError as error:
            raise ShardError(f"could not inspect existing artifact {path}: {error}") from error
        if observed != value:
            raise ShardError(f"refusing to replace differing existing artifact: {path}")
        return
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        link_no_replace(temporary, path)
        if path.read_bytes() != value:
            raise ShardError(f"concurrent no-replace publication differs: {path}")
    finally:
        temporary.unlink(missing_ok=True)


def materialize_shard(
    destination: Path,
    source: Snapshot,
    columns: list[str],
    reference_snps: set[str],
) -> tuple[IntersectionMetrics, bytes]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        return scan_source_intersection(source, columns, reference_snps)
    temporary = destination.with_name(f".{destination.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", compresslevel=6, fileobj=raw, mtime=0) as output:
                metrics, header = scan_source_intersection(
                    source, columns, reference_snps, output=output,
                )
            raw.flush()
            os.fsync(raw.fileno())
        link_no_replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return metrics, header


def input_info_bytes(context: Context, chromosome: int) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=INPUT_INFO_FIELDS, delimiter="\t", lineterminator="\n",
    )
    writer.writeheader()
    for original in context.input_rows:
        row = dict(original)
        row["filename"] = relative_path(
            context.root,
            context.output_root / f"chr{chromosome:02d}" / f"{row['phenotype']}.sumstats.gz",
        )
        writer.writerow(row)
    return buffer.getvalue().encode("utf-8")


def shard_record(
    context: Context,
    chromosome: int,
    phenotype: str,
    reference: Snapshot,
    reference_summary: dict[str, object],
    shard: Snapshot,
    metrics: IntersectionMetrics,
) -> dict[str, object]:
    source = context.source_snapshots[phenotype]
    return {
        "chromosome": chromosome,
        "phenotype": phenotype,
        "path": shard.relative,
        "bytes": shard.bytes,
        "sha256": shard.sha256,
        **metrics.record(),
        "source_path": source.relative,
        "source_bytes": source.bytes,
        "source_sha256": source.sha256,
        "reference_info_path": reference.relative,
        "reference_info_bytes": reference.bytes,
        "reference_info_sha256": reference.sha256,
        "reference_unique_snp_count": reference_summary["unique_snp_count"],
        "selection_rule": "NORMALIZED_LOWERCASE_SNP_ID_IN_COMPLETE_SEALED_CHROMOSOME_REFERENCE",
        "row_transform": "NONE_EXACT_SOURCE_ROW_BYTES_IN_ORIGINAL_ORDER",
    }


def build_or_validate_family(context: Context, *, create_missing: bool) -> dict[str, object]:
    reference_summaries: list[dict[str, object]] = []
    shard_records: list[dict[str, object]] = []
    input_info_records: list[dict[str, object]] = []
    for chromosome in CHROMOSOMES:
        reference = reference_snapshot(context, chromosome)
        reference_snps, reference_summary = load_reference_snps(reference, chromosome)
        reference_summaries.append(reference_summary)
        chromosome_dir = context.output_root / f"chr{chromosome:02d}"
        for row in context.input_rows:
            phenotype = row["phenotype"]
            destination = chromosome_dir / f"{phenotype}.sumstats.gz"
            if not destination.exists() and not create_missing:
                raise ShardError(f"chromosome shard is missing: {relative_path(context.root, destination)}")
            if create_missing:
                expected, header = materialize_shard(
                    destination,
                    context.source_snapshots[phenotype],
                    context.source_columns[phenotype],
                    reference_snps,
                )
            else:
                expected, header = scan_source_intersection(
                    context.source_snapshots[phenotype],
                    context.source_columns[phenotype],
                    reference_snps,
                )
            live_shard = snapshot(destination, context.root)
            observed = scan_shard(
                live_shard, context.source_columns[phenotype], header, reference_snps,
            )
            if observed != expected:
                raise ShardError(
                    f"shard is not the exact ordered source/reference intersection: {live_shard.relative}"
                )
            shard_records.append(
                shard_record(
                    context, chromosome, phenotype, reference, reference_summary,
                    live_shard, observed,
                )
            )
        info_path = chromosome_dir / "lava_input_info.tsv"
        expected_info = input_info_bytes(context, chromosome)
        if create_missing:
            publish_bytes_no_replace(info_path, expected_info)
        elif not info_path.is_file():
            raise ShardError(f"chromosome input-info is missing: {relative_path(context.root, info_path)}")
        if info_path.read_bytes() != expected_info:
            raise ShardError(f"chromosome input-info differs from frozen values: {info_path}")
        info_snapshot = snapshot(info_path, context.root)
        input_info_records.append(
            info_snapshot.record() | {"chromosome": chromosome, "trait_count": len(context.input_rows)}
        )

    source_records = [
        context.source_snapshots[row["phenotype"]].record()
        | {
            "phenotype": row["phenotype"],
            "observed_columns": context.source_columns[row["phenotype"]],
        }
        for row in context.input_rows
    ]
    total_rows = sum(int(item["data_row_count"]) for item in shard_records)
    return {
        "schema_version": SCHEMA_VERSION,
        "analysis_id": context.policy["analysis_id"],
        "artifact_role": "EXACT_CHROMOSOME_SCOPED_LAVA_INPUT_MATERIALIZATION",
        "selection_timing": "BEFORE_LOCAL_RESULT_ACCESS",
        "scientific_filtering": "NONE_ONLY_EXACT_SEALED_REFERENCE_SNP_INTERSECTION",
        "chromosomes": list(CHROMOSOMES),
        "trait_order": [row["phenotype"] for row in context.input_rows],
        "bindings": {
            "policy_path": context.policy_snapshot.relative,
            "policy_sha256": context.policy_snapshot.sha256,
            "local_input_lock_path": context.input_lock_snapshot.relative,
            "local_input_lock_sha256": context.input_lock_snapshot.sha256,
            "source_input_info_path": context.input_info_snapshot.relative,
            "source_input_info_sha256": context.input_info_snapshot.sha256,
            "source_input_provenance_path": context.input_provenance_snapshot.relative,
            "source_input_provenance_sha256": context.input_provenance_snapshot.sha256,
            "reference_provenance_path": context.reference_provenance_snapshot.relative,
            "reference_provenance_sha256": context.reference_provenance_snapshot.sha256,
            "script_path": context.script_snapshot.relative,
            "script_sha256": context.script_snapshot.sha256,
        },
        "sources": source_records,
        "reference_info": reference_summaries,
        "shards": shard_records,
        "chromosome_input_info": input_info_records,
        "totals": {
            "chromosome_count": len(CHROMOSOMES),
            "trait_count": len(context.input_rows),
            "shard_count": len(shard_records),
            "chromosome_input_info_count": len(input_info_records),
            "selected_source_row_count_across_all_shards": total_rows,
        },
    }


def lock_bytes(payload: dict[str, object]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")


def validate_output_identities(context: Context, payload: dict[str, object]) -> None:
    records = list(payload["shards"]) + list(payload["chromosome_input_info"])
    for record in records:
        path = safe_path(context.root, str(record["path"]), "chromosome-input artifact")
        observed = snapshot(path, context.root)
        if observed.bytes != int(record["bytes"]) or observed.sha256 != str(record["sha256"]):
            raise ShardError(f"chromosome-input artifact changed before sealing: {observed.relative}")


def validate_frozen_build_inputs(context: Context) -> None:
    for original in (
        context.policy_snapshot, context.input_lock_snapshot, context.input_info_snapshot,
        context.input_provenance_snapshot, context.reference_provenance_snapshot,
        context.script_snapshot,
    ):
        observed = snapshot(
            original.path, context.root, allow_outside_root=not original.path.is_relative_to(context.root),
        )
        if observed != original:
            raise ShardError(f"build input changed during chromosome materialization: {original.path}")
    for phenotype, original in context.source_snapshots.items():
        if snapshot(original.path, context.root) != original:
            raise ShardError(f"source summary statistics changed during chromosome materialization: {phenotype}")
    for chromosome in CHROMOSOMES:
        reference_snapshot(context, chromosome)


def verify(context: Context) -> dict[str, object]:
    if not context.output_lock.is_file():
        raise ShardError(f"chromosome-input provenance lock is missing: {context.output_lock}")
    observed_snapshot = snapshot(context.output_lock, context.root)
    observed = load_json(observed_snapshot, "Track B LAVA chromosome-input provenance")
    expected = build_or_validate_family(context, create_missing=False)
    if observed != expected:
        raise ShardError("chromosome-input provenance differs from the fully revalidated shard family")
    validate_output_identities(context, expected)
    validate_frozen_build_inputs(context)
    return expected


def build(context: Context) -> dict[str, object]:
    if context.output_lock.exists():
        return verify(context)
    payload = build_or_validate_family(context, create_missing=True)
    validate_output_identities(context, payload)
    validate_frozen_build_inputs(context)
    publish_bytes_no_replace(context.output_lock, lock_bytes(payload))
    observed = load_json(snapshot(context.output_lock, context.root), "Track B LAVA chromosome-input provenance")
    if observed != payload:
        raise ShardError("published chromosome-input provenance differs from materialized family")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="fully verify existing shards without writing")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT, help=argparse.SUPPRESS)
    args = parser.parse_args()
    context = load_context(args.root)
    payload = verify(context) if args.verify else build(context)
    action = "VERIFIED" if args.verify else "READY"
    print(
        f"TRACK_B_LAVA_CHROMOSOME_INPUTS_{action} "
        f"chromosomes={payload['totals']['chromosome_count']} "
        f"traits={payload['totals']['trait_count']} shards={payload['totals']['shard_count']} "
        f"lock={relative_path(context.root, context.output_lock)}"
    )


if __name__ == "__main__":
    try:
        main()
    except ShardError as error:
        raise SystemExit(f"ERROR: {error}") from error
