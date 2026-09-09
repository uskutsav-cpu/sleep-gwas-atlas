#!/usr/bin/env python3
"""Checksum-verifying readers for large local or version-pinned remote gzip objects."""

from __future__ import annotations

import gzip
import hashlib
import io
import ssl
import urllib.parse
import urllib.request
from contextlib import contextmanager
from pathlib import Path
from typing import BinaryIO, Iterator, TextIO


class DigestingReader(io.RawIOBase):
    """Minimal binary wrapper that hashes/counts bytes consumed by gzip.GzipFile."""

    def __init__(self, raw: BinaryIO):
        self.raw = raw
        self.md5 = hashlib.md5(usedforsecurity=False)
        self.sha256 = hashlib.sha256()
        self.bytes_read = 0

    def read(self, size: int = -1) -> bytes:
        block = self.raw.read(size)
        if block:
            self.md5.update(block)
            self.sha256.update(block)
            self.bytes_read += len(block)
        return block

    def readinto(self, buffer: bytearray | memoryview) -> int:
        block = self.raw.read(len(buffer))
        if not block:
            return 0
        self.md5.update(block)
        self.sha256.update(block)
        self.bytes_read += len(block)
        buffer[: len(block)] = block
        return len(block)

    def readable(self) -> bool:
        return True

    def close(self) -> None:
        if not self.closed:
            self.raw.close()
        super().close()


def _parse_expected_md5(checksum: str | None) -> str | None:
    if checksum is None or checksum in {"", "NA", "UNKNOWN"}:
        return None
    if checksum.startswith("md5:"):
        value = checksum.removeprefix("md5:").lower()
    elif checksum.startswith("etag:") and "-" not in checksum:
        value = checksum.removeprefix("etag:").lower()
    else:
        raise ValueError(
            "streaming verification requires md5:<hex> or a single-part etag:<hex>; "
            f"found {checksum!r}"
        )
    if len(value) != 32 or any(character not in "0123456789abcdef" for character in value):
        raise ValueError(f"invalid MD5 checksum: {checksum!r}")
    return value


def _open_remote(url: str, timeout_seconds: int) -> BinaryIO:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https":
        raise ValueError(f"remote production sources must use https, found {parsed.scheme!r}")
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "sleep-gwas-atlas-extension/1.0"},
        method="GET",
    )
    try:
        import certifi

        context = ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        context = ssl.create_default_context()
    return urllib.request.urlopen(  # nosec B310: HTTPS enforced above
        request, timeout=timeout_seconds, context=context
    )


@contextmanager
def open_verified_gzip_text(
    *,
    local_path: Path | None = None,
    url: str | None = None,
    expected_md5: str | None = None,
    expected_size_bytes: int | None = None,
    timeout_seconds: int = 120,
) -> Iterator[tuple[TextIO, dict[str, object]]]:
    """Yield decompressed text and verify the entire compressed object at EOF.

    The caller must consume the text iterator to EOF. Verification happens only after
    the context exits successfully; a mismatch raises and downstream temporary output
    must therefore not be promoted before this context has closed.
    """

    if (local_path is None) == (url is None):
        raise ValueError("provide exactly one of local_path or url")
    locked_md5 = _parse_expected_md5(expected_md5)
    raw: BinaryIO
    source_label: str
    response_headers: dict[str, str] = {}
    if local_path is not None:
        raw = local_path.open("rb")
        source_label = str(local_path)
    else:
        assert url is not None
        raw = _open_remote(url, timeout_seconds)
        source_label = url
        response_headers = {key.lower(): value for key, value in raw.headers.items()}

    digesting = DigestingReader(raw)
    buffered = io.BufferedReader(digesting, buffer_size=1024 * 1024)
    receipt: dict[str, object] = {
        "source": source_label,
        "expected_md5": locked_md5,
        "expected_size_bytes": expected_size_bytes,
        "http_etag": response_headers.get("etag", "NA"),
        "http_last_modified": response_headers.get("last-modified", "NA"),
        "http_version_id": response_headers.get("x-amz-version-id", "NA"),
    }
    try:
        with gzip.GzipFile(fileobj=buffered, mode="rb") as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                yield text, receipt
        observed_md5 = digesting.md5.hexdigest()
        observed_sha256 = digesting.sha256.hexdigest()
        receipt.update(
            {
                "observed_md5": observed_md5,
                "observed_sha256": observed_sha256,
                "observed_size_bytes": digesting.bytes_read,
            }
        )
        if expected_size_bytes is not None and digesting.bytes_read != expected_size_bytes:
            raise ValueError(
                f"compressed size mismatch for {source_label}: "
                f"observed {digesting.bytes_read}, expected {expected_size_bytes}"
            )
        if locked_md5 is not None and observed_md5 != locked_md5:
            raise ValueError(
                f"MD5 mismatch for {source_label}: observed {observed_md5}, expected {locked_md5}"
            )
        receipt["verification_status"] = "PASS"
    finally:
        buffered.close()
