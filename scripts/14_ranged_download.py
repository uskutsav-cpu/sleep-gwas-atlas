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
import hashlib
import os
from pathlib import Path
import re
import sys
from urllib.request import Request, urlopen


CONTENT_RANGE = re.compile(r"bytes (\d+)-(\d+)/(\d+)$")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    args = parser.parse_args()
    if args.expected_bytes <= 0 or args.chunk_bytes <= 0 or args.max_chunks <= 0:
        fail("expected bytes, chunk bytes, and max chunks must be positive")

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

    retrieved = 0
    while offset < args.expected_bytes and retrieved < args.max_chunks:
        end = min(offset + args.chunk_bytes, args.expected_bytes) - 1
        request = Request(args.url, headers={"Range": f"bytes={offset}-{end}"})
        with urlopen(request, timeout=120) as response:
            content_range = response.headers.get("Content-Range", "")
            match = CONTENT_RANGE.fullmatch(content_range.strip())
            if response.status != 206 or match is None:
                fail(
                    "range request was not honored: expected HTTP 206 with a valid "
                    f"Content-Range, got status={response.status}, Content-Range={content_range!r}"
                )
            start, actual_end, total = (int(value) for value in match.groups())
            if (start, actual_end, total) != (offset, end, args.expected_bytes):
                fail(
                    "unexpected Content-Range: expected "
                    f"bytes {offset}-{end}/{args.expected_bytes}, got {content_range!r}"
                )
            remaining = end - offset + 1
            with partial.open("ab") as handle:
                while remaining:
                    block = response.read(min(1024 * 1024, remaining))
                    if not block:
                        fail(f"response ended early with {remaining} bytes still expected")
                    handle.write(block)
                    remaining -= len(block)
                if response.read(1):
                    fail("range response contained more bytes than declared")
                handle.flush()
                os.fsync(handle.fileno())
        offset = end + 1
        retrieved += 1
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
