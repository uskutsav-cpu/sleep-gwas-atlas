"""Streaming IO, immutable artifacts, checksums and sandboxed write paths."""
from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import tempfile
import time
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, TextIO


class ContractError(ValueError):
    """An input or scientific contract failed; do not silently continue."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def json_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def read_json(path: str | Path) -> Any:
    with Path(path).open(encoding="utf-8") as f:
        def pairs_hook(items):
            result = {}
            for key, value in items:
                require(key not in result, f"Duplicate JSON key: {key}")
                result[key] = value
            return result
        def invalid_constant(value):
            raise ContractError(f"Nonfinite JSON constant: {value}")
        return json.load(f, object_pairs_hook=pairs_hook, parse_constant=invalid_constant)


def check_hash(path: str | Path, expected: str) -> None:
    require(bool(re.fullmatch(r"[a-fA-F0-9]{64}", expected or "")),
            f"Missing/invalid SHA-256 for {path}")
    actual = sha256(path)
    require(actual == expected.lower(), f"Checksum mismatch for {path}: {actual}")


def safe_id(value: str) -> str:
    require(bool(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,119}", value or "")),
            f"Unsafe identifier: {value!r}")
    require(value not in {".", ".."}, "Invalid identifier")
    return value


def safe_write_path(root: str | Path, relative: str | Path) -> Path:
    root = Path(root).resolve()
    relative = Path(relative)
    require(not relative.is_absolute(), f"Absolute output path forbidden: {relative}")
    target = root / relative
    # Check lexical and resolved ancestry. In particular, reject symlink escape.
    require(".." not in relative.parts, f"Parent traversal forbidden: {relative}")
    require(target.resolve().is_relative_to(root), f"Output escapes root: {target}")
    cur = root
    for part in relative.parts:
        cur /= part
        require(not cur.is_symlink(), f"Symlink output forbidden: {cur}")
    return target


def ensure_free(path: str | Path, needed_bytes: int) -> None:
    path = Path(path)
    while not path.exists():
        path = path.parent
    free = shutil.disk_usage(path).free
    require(free >= needed_bytes,
            f"BLOCKED_BY_STORAGE: {needed_bytes / 2**30:.2f} GiB required; "
            f"{free / 2**30:.2f} GiB available")


@contextlib.contextmanager
def open_text(path: str | Path) -> Iterator[TextIO]:
    """Read gzip/bgzip by magic, not by unreliable filename extension."""
    with Path(path).open("rb") as raw:
        magic = raw.read(2)
        raw.seek(0)
        stream = gzip.GzipFile(fileobj=raw) if magic == b"\x1f\x8b" else raw
        with io.TextIOWrapper(stream, encoding="utf-8-sig", newline="") as text:
            yield text


@contextlib.contextmanager
def atomic_text(path: str | Path, *, immutable: bool = True,
                compress: bool | None = None) -> Iterator[TextIO]:
    """Never promote an interrupted file. Gzip is deterministic (mtime=0)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    require(not path.is_symlink(), f"Refusing symlink output: {path}")
    if immutable:
        require(not path.exists(), f"Refusing to overwrite existing artifact: {path}")
    fd, tmp = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".partial", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as raw:
            gz = (path.suffix == ".gz") if compress is None else compress
            stream = gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) if gz else raw
            wrapper = io.TextIOWrapper(stream, encoding="utf-8", newline="")
            try:
                yield wrapper
                wrapper.flush()
            finally:
                if gz:
                    wrapper.close()
                else:
                    wrapper.detach()
            raw.flush()
            os.fsync(raw.fileno())
        if immutable:
            # link() is an atomic no-clobber publication; same filesystem.
            os.link(tmp, path)
            os.unlink(tmp)
        else:
            os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def write_json(path: str | Path, value: Any, *, immutable: bool = True) -> None:
    with atomic_text(path, immutable=immutable, compress=False) as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write("\n")


def read_tsv(path: str | Path, required: Iterable[str] = ()) -> Iterator[dict[str, str]]:
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        require(reader.fieldnames is not None, f"Empty table: {path}")
        require(len(set(reader.fieldnames)) == len(reader.fieldnames), f"Duplicate columns: {path}")
        require(set(required) <= set(reader.fieldnames),
                f"{path}: missing columns {set(required) - set(reader.fieldnames)}")
        for line, row in enumerate(reader, 2):
            require(None not in row and all(v is not None for v in row.values()),
                    f"{path}:{line}: malformed field count")
            yield row


def write_tsv(path: str | Path, fields: Iterable[str], rows: Iterable[Mapping[str, Any]],
              *, immutable: bool = True) -> int:
    count = 0
    with atomic_text(path, immutable=immutable) as f:
        writer = csv.DictWriter(f, fieldnames=list(fields), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
            count += 1
    return count


def file_record(path: str | Path) -> dict[str, Any]:
    p = Path(path).resolve()
    require(p.is_file(), f"Missing input: {p}")
    return {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size}


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


@contextlib.contextmanager
def exclusive_lock(path: str | Path) -> Iterator[None]:
    """Single writer. Stale locks require explicit human inspection/removal."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise ContractError(f"LOCKED: {path}; inspect owner before removing a stale lock") from exc
    try:
        with os.fdopen(fd, "w") as f:
            json.dump({"pid": os.getpid(), "created": utc_now()}, f)
        yield
    finally:
        path.unlink(missing_ok=True)
