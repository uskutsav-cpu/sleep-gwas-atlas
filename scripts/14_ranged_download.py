#!/usr/bin/env python3
"""Download a large public file by validated HTTP byte ranges.

Some public Google Drive releases are too large to safely restart in a single
interactive transfer.  This downloader writes only to ``OUT.partial`` until
every requested range has the expected HTTP ``206 Content-Range`` response,
then atomically promotes the completed file.  It deliberately refuses a
server that returns a whole-file ``200`` response to a range request, avoiding
the duplicated/corrupted archives that a blind ``curl --continue-at`` can
create on such endpoints.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess


CONTENT_RANGE = re.compile(r"bytes (\d+)-(\d+)/(\d+)$")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fetch_range(url: str, start: int, end: int, total: int, destination: Path) -> Path:
    """Fetch one exact byte range to an atomic part file."""
    expected_size = end - start + 1
    if destination.is_file() and destination.stat().st_size == expected_size:
        return destination
    temporary = destination.with_name(f".{destination.name}.partial")
    headers = destination.with_name(f".{destination.name}.headers")
    temporary.unlink(missing_ok=True)
    headers.unlink(missing_ok=True)
    curl = shutil.which("curl")
    if curl is None:
        fail("curl is required for certificate-aware ranged downloads")
    try:
        result = subprocess.run(
            [
                curl,
                "--fail",
                "--silent",
                "--show-error",
                "--location",
                "--retry",
                "3",
                "--range",
                f"{start}-{end}",
                "--dump-header",
                str(headers),
                "--output",
                str(temporary),
                url,
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode:
            fail(f"curl range request failed: {result.stderr.strip()}")
        header_text = headers.read_text(encoding="iso-8859-1")
        statuses = re.findall(r"^HTTP/\S+\s+(\d{3})", header_text, re.MULTILINE)
        content_ranges = re.findall(
            r"^Content-Range:\s*([^\r\n]+)", header_text, re.MULTILINE | re.IGNORECASE
        )
        status = int(statuses[-1]) if statuses else 0
        content_range = content_ranges[-1].strip() if content_ranges else ""
        match = CONTENT_RANGE.fullmatch(content_range)
        if status != 206 or match is None:
            fail(
                "range request was not honored: expected HTTP 206 with a valid "
                f"Content-Range, got status={status}, Content-Range={content_range!r}"
            )
        actual = tuple(int(value) for value in match.groups())
        if actual != (start, end, total):
            fail(
                f"unexpected Content-Range: expected bytes {start}-{end}/{total}, "
                f"got {content_range!r}"
            )
        if temporary.stat().st_size != expected_size:
            fail(
                f"range body has {temporary.stat().st_size} bytes, expected {expected_size}"
            )
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        headers.unlink(missing_ok=True)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected-bytes", type=int, required=True)
    parser.add_argument("--expected-sha256")
    parser.add_argument("--chunk-bytes", type=int, default=128 * 1024 * 1024)
    parser.add_argument(
        "--max-chunks", type=int, default=1,
        help="Number of validated chunks to retrieve in this invocation (default: 1).",
    )
    parser.add_argument(
        "--workers", type=int, default=1,
        help="Concurrent validated range requests (default: 1).",
    )
    args = parser.parse_args()
    if (
        args.expected_bytes <= 0
        or args.chunk_bytes <= 0
        or args.max_chunks <= 0
        or args.workers <= 0
    ):
        fail("expected bytes, chunk bytes, max chunks, and workers must be positive")

    output = args.out
    partial = output.with_name(f"{output.name}.partial")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        if output.stat().st_size != args.expected_bytes:
            fail(f"existing output has {output.stat().st_size} bytes, expected {args.expected_bytes}")
        actual_sha256 = sha256(output)
        if args.expected_sha256 and actual_sha256 != args.expected_sha256:
            fail("existing output SHA-256 does not match the registered value")
        print(f"Already complete: {output} ({actual_sha256})")
        return

    offset = partial.stat().st_size if partial.exists() else 0
    if offset > args.expected_bytes:
        fail(f"partial file has {offset} bytes, expected at most {args.expected_bytes}")

    ranges = []
    cursor = offset
    while cursor < args.expected_bytes and len(ranges) < args.max_chunks:
        end = min(cursor + args.chunk_bytes, args.expected_bytes) - 1
        ranges.append((cursor, end))
        cursor = end + 1

    if ranges:
        parts_dir = output.with_name(f".{output.name}.range-parts")
        parts_dir.mkdir(exist_ok=True)
        parts = {
            (start, end): parts_dir / f"{start:020d}-{end:020d}.part"
            for start, end in ranges
        }
        with ThreadPoolExecutor(max_workers=min(args.workers, len(ranges))) as executor:
            futures = {
                executor.submit(
                    fetch_range, args.url, start, end, args.expected_bytes, parts[(start, end)]
                ): (start, end)
                for start, end in ranges
            }
            for future in as_completed(futures):
                start, end = futures[future]
                future.result()
                print(f"Validated range {start:,}-{end:,}", flush=True)

        with partial.open("ab") as output_handle:
            for start, end in ranges:
                part = parts[(start, end)]
                with part.open("rb") as input_handle:
                    for block in iter(lambda: input_handle.read(1024 * 1024), b""):
                        output_handle.write(block)
                offset = end + 1
                part.unlink()
            output_handle.flush()
            os.fsync(output_handle.fileno())
        try:
            parts_dir.rmdir()
        except OSError:
            pass
    print(f"Validated {offset:,} / {args.expected_bytes:,} bytes", flush=True)

    if offset != args.expected_bytes:
        print(f"Partial download retained at {partial}; run again to resume.")
        return

    actual_sha256 = sha256(partial)
    if args.expected_sha256 and actual_sha256 != args.expected_sha256:
        fail("completed file SHA-256 does not match the registered value")
    os.replace(partial, output)
    print(f"Completed {output}")
    print(f"SHA-256: {actual_sha256}")


if __name__ == "__main__":
    main()
