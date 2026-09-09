#!/usr/bin/env python3
"""Shared fail-closed I/O for Track B phases 24 and 25.

The helpers deliberately avoid implicit directories, symbolic links, in-place
replacement, and Python bytecode.  Read-only preflight callers therefore have
no filesystem side effects, while publication is idempotent after a crash only
when every already-present byte is exactly the byte that would be published.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import stat
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


HEX64 = re.compile(r"^[0-9a-f]{64}$")
MAX_STRUCTURED_BYTES = 256 * 1024**2


class ContractError(RuntimeError):
    """An immutable input, schema, or publication contract was violated."""


class UpstreamBlocked(ContractError):
    """Required real upstream results are absent or explicitly non-runnable."""


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def digest_json(value: Any) -> str:
    return digest_bytes(canonical_json_bytes(value))


def _stat_identity(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def safe_path(
    root: Path, relative: str | Path, label: str, *, must_exist: bool = False,
) -> Path:
    """Resolve a repository-relative path without following any symlink."""

    relative = Path(relative)
    if relative.is_absolute() or not relative.parts or ".." in relative.parts:
        raise ContractError(f"unsafe {label} path: {relative}")
    try:
        root_real = root.resolve(strict=True)
    except OSError as error:
        raise ContractError(f"repository root is unavailable: {error}") from error
    current = root
    for offset, part in enumerate(relative.parts):
        current /= part
        try:
            observed = os.lstat(current)
        except FileNotFoundError:
            break
        except OSError as error:
            raise ContractError(f"could not inspect {label} {relative}: {error}") from error
        if stat.S_ISLNK(observed.st_mode):
            raise ContractError(f"{label} contains a symbolic link: {relative}")
        if offset < len(relative.parts) - 1 and not stat.S_ISDIR(observed.st_mode):
            raise ContractError(f"{label} has a non-directory ancestor: {relative}")
    path = root / relative
    try:
        path.resolve(strict=False).relative_to(root_real)
    except (OSError, RuntimeError, ValueError) as error:
        raise ContractError(f"{label} escapes repository: {relative}") from error
    if must_exist and not path.exists():
        raise ContractError(f"required {label} is absent: {relative}")
    return path


def stable_bytes(
    root: Path,
    relative: str | Path,
    *,
    allow_empty: bool = False,
    maximum_bytes: int = MAX_STRUCTURED_BYTES,
) -> tuple[bytes, dict[str, object]]:
    """Read one regular file and reject named-inode mutation during the read."""

    path = safe_path(root, relative, "artifact", must_exist=True)
    relative_string = str(Path(relative))
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        named_before = os.lstat(path)
        if not stat.S_ISREG(named_before.st_mode):
            raise ContractError(f"artifact is not a regular file: {relative_string}")
        if named_before.st_size == 0 and not allow_empty:
            raise ContractError(f"artifact is empty: {relative_string}")
        if named_before.st_size > maximum_bytes:
            raise ContractError(f"structured artifact is too large: {relative_string}")
        descriptor = os.open(path, flags)
        try:
            opened = os.fstat(descriptor)
            chunks: list[bytes] = []
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                chunks.append(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        named_after = os.lstat(path)
    except OSError as error:
        raise ContractError(f"could not read stable artifact {relative_string}: {error}") from error
    identities = {_stat_identity(item) for item in (named_before, opened, after, named_after)}
    if len(identities) != 1:
        raise ContractError(f"artifact changed while reading: {relative_string}")
    payload = b"".join(chunks)
    return payload, {
        "path": relative_string,
        "bytes": len(payload),
        "sha256": digest_bytes(payload),
    }


def stable_identity(root: Path, relative: str | Path) -> dict[str, object]:
    """Hash a possibly large regular file without retaining it in memory."""

    path = safe_path(root, relative, "artifact", must_exist=True)
    relative_string = str(Path(relative))
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    digest = hashlib.sha256()
    try:
        named_before = os.lstat(path)
        if not stat.S_ISREG(named_before.st_mode) or named_before.st_size <= 0:
            raise ContractError(f"artifact is not a non-empty regular file: {relative_string}")
        descriptor = os.open(path, flags)
        try:
            opened = os.fstat(descriptor)
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        named_after = os.lstat(path)
    except OSError as error:
        raise ContractError(f"could not hash stable artifact {relative_string}: {error}") from error
    if len({_stat_identity(item) for item in (named_before, opened, after, named_after)}) != 1:
        raise ContractError(f"artifact changed while hashing: {relative_string}")
    return {
        "path": relative_string,
        "bytes": named_after.st_size,
        "sha256": digest.hexdigest(),
    }


def stable_json(root: Path, relative: str | Path) -> tuple[dict[str, Any], dict[str, object]]:
    payload, identity = stable_bytes(root, relative)
    try:
        value = json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError(f"invalid JSON artifact {relative}: {error}") from error
    if not isinstance(value, dict):
        raise ContractError(f"JSON artifact must contain an object: {relative}")
    return value, identity


def stable_tsv(
    root: Path, relative: str | Path, expected_fields: Sequence[str] | None = None,
) -> tuple[list[str], list[dict[str, str]], dict[str, object]]:
    payload, identity = stable_bytes(root, relative)
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ContractError(f"invalid UTF-8 table {relative}: {error}") from error
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    fields = list(reader.fieldnames or [])
    if expected_fields is not None and fields != list(expected_fields):
        raise ContractError(f"table header differs from contract: {relative}")
    try:
        rows = list(reader)
    except csv.Error as error:
        raise ContractError(f"invalid TSV artifact {relative}: {error}") from error
    if any(None in row for row in rows):
        raise ContractError(f"table contains an over-wide row: {relative}")
    return fields, rows, identity


def table_bytes(fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output,
        fieldnames=list(fields),
        delimiter="\t",
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    for row in rows:
        if list(row) != list(fields):
            raise ContractError("row field order differs from the output contract")
        writer.writerow(row)
    return output.getvalue().encode("utf-8")


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        raise ContractError(f"invalid SHA-256 for {label}")
    return value


def _ensure_parent(root: Path, relative: Path) -> Path:
    parent_relative = relative.parent
    if str(parent_relative) == ".":
        return root
    current = root
    for part in parent_relative.parts:
        current /= part
        try:
            info = os.lstat(current)
        except FileNotFoundError:
            try:
                os.mkdir(current, 0o755)
            except FileExistsError:
                pass
            info = os.lstat(current)
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise ContractError(f"output parent is not a real directory: {parent_relative}")
    safe_path(root, parent_relative, "output parent", must_exist=True)
    return root / parent_relative


def no_replace_publish(
    root: Path,
    payloads: Sequence[tuple[str | Path, bytes]],
) -> dict[str, str]:
    """Publish exact payloads without replacement; matching partials resume."""

    if not payloads:
        raise ContractError("no output payloads were supplied")
    observed_paths: set[Path] = set()
    for relative, _payload in payloads:
        key = Path(relative)
        if key in observed_paths:
            raise ContractError(f"duplicate output path: {key}")
        observed_paths.add(key)

    state: dict[str, str] = {}
    for relative_value, payload in payloads:
        relative = Path(relative_value)
        path = safe_path(root, relative, "output")
        try:
            info = os.lstat(path)
        except FileNotFoundError:
            continue
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
            raise ContractError(f"existing output is not a real regular file: {relative}")
        existing, _ = stable_bytes(root, relative, allow_empty=True)
        if existing != payload:
            raise ContractError(f"refusing to replace divergent existing output: {relative}")
        state[str(relative)] = "ALREADY_PRESENT_BYTE_IDENTICAL"

    for relative_value, payload in payloads:
        relative = Path(relative_value)
        if str(relative) in state:
            continue
        parent = _ensure_parent(root, relative)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{relative.name}.", suffix=".tmp", dir=parent,
        )
        temporary = Path(temporary_name)
        try:
            with os.fdopen(descriptor, "wb") as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            target = root / relative
            try:
                os.link(temporary, target, follow_symlinks=False)
            except FileExistsError:
                existing, _ = stable_bytes(root, relative, allow_empty=True)
                if existing != payload:
                    raise ContractError(
                        f"concurrent divergent output won publication: {relative}"
                    )
                state[str(relative)] = "CONCURRENT_BYTE_IDENTICAL"
            else:
                state[str(relative)] = "PUBLISHED_NO_REPLACE"
                directory_fd = os.open(parent, os.O_RDONLY)
                try:
                    os.fsync(directory_fd)
                finally:
                    os.close(directory_fd)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
    return state
