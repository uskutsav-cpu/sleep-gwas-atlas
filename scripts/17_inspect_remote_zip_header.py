#!/usr/bin/env python3
"""Inspect the first records of a registered remote ZIP member by byte range.

This curation helper deliberately avoids downloading a full GWAS archive.  It
validates every HTTP range response, parses the classic ZIP central directory,
and inflates only enough of one named member to print its first text records.
It is inspection-only: the output is not an accepted raw GWAS artifact.
"""

from __future__ import annotations

import argparse
import re
import ssl
import struct
import sys
import zlib
from urllib.request import Request, urlopen

try:
    import certifi
except ImportError:  # use the platform trust store outside the pinned environment
    certifi = None


EOCD_SIGNATURE = b"PK\x05\x06"
CENTRAL_SIGNATURE = b"PK\x01\x02"
LOCAL_SIGNATURE = b"PK\x03\x04"
TLS_CONTEXT = ssl.create_default_context(cafile=certifi.where() if certifi else None)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--expected-bytes", required=True, type=int)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--member")
    action.add_argument("--list-members", action="store_true")
    parser.add_argument("--lines", type=int, default=2)
    parser.add_argument("--match", help="print the first decoded line containing this text")
    parser.add_argument("--max-compressed-bytes", type=int, default=1_048_576)
    return parser.parse_args()


def fetch_range(url: str, start: int, end: int, total: int) -> bytes:
    if not 0 <= start <= end < total:
        raise ValueError(f"invalid byte range {start}-{end} for {total}-byte file")
    request = Request(
        url,
        headers={
            "Range": f"bytes={start}-{end}",
            "User-Agent": "sleep-gwas-atlas/remote-zip-header-inspector",
        },
    )
    with urlopen(request, timeout=60, context=TLS_CONTEXT) as response:
        content_range = response.headers.get("Content-Range", "")
        match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", content_range)
        if response.status != 206 or not match:
            raise RuntimeError(
                "server did not honor byte range: "
                f"status={response.status}, Content-Range={content_range!r}"
            )
        observed = tuple(map(int, match.groups()))
        if observed != (start, end, total):
            raise RuntimeError(
                f"unexpected Content-Range {observed}; expected {(start, end, total)}"
            )
        payload = response.read()
    expected = end - start + 1
    if len(payload) != expected:
        raise RuntimeError(f"short range response: expected {expected}, got {len(payload)}")
    return payload


def central_directory(url: str, total: int) -> bytes:
    tail_size = min(total, 65_557)
    tail_start = total - tail_size
    tail = fetch_range(url, tail_start, total - 1, total)
    offset = tail.rfind(EOCD_SIGNATURE)
    if offset < 0 or len(tail) - offset < 22:
        raise RuntimeError("classic ZIP end-of-central-directory record not found")
    fields = struct.unpack_from("<4s4H2LH", tail, offset)
    _, disk, central_disk, disk_entries, entries, size, start, comment_length = fields
    if disk != 0 or central_disk != 0 or disk_entries != entries:
        raise RuntimeError("multi-disk ZIP archives are unsupported")
    if entries == 0xFFFF or size == 0xFFFFFFFF or start == 0xFFFFFFFF:
        raise RuntimeError("ZIP64 central directories are unsupported")
    if offset + 22 + comment_length != len(tail):
        raise RuntimeError("ZIP end-of-central-directory length mismatch")
    if start + size > total:
        raise RuntimeError("central directory lies outside the registered byte count")
    return fetch_range(url, start, start + size - 1, total)


def member_entries(central: bytes) -> list[tuple[str, tuple[int, int, int, int]]]:
    cursor = 0
    entries = []
    while cursor < len(central):
        if central[cursor : cursor + 4] != CENTRAL_SIGNATURE or cursor + 46 > len(central):
            raise RuntimeError(f"malformed central-directory entry at byte {cursor}")
        flags, method = struct.unpack_from("<HH", central, cursor + 8)
        compressed_size, uncompressed_size = struct.unpack_from("<LL", central, cursor + 20)
        name_len, extra_len, comment_len = struct.unpack_from("<HHH", central, cursor + 28)
        local_offset = struct.unpack_from("<L", central, cursor + 42)[0]
        end = cursor + 46 + name_len + extra_len + comment_len
        if end > len(central):
            raise RuntimeError("truncated central-directory entry")
        encoding = "utf-8" if flags & 0x800 else "cp437"
        name = central[cursor + 46 : cursor + 46 + name_len].decode(encoding)
        entries.append(
            (name, (local_offset, method, compressed_size, uncompressed_size))
        )
        cursor = end
    if cursor != len(central):
        raise RuntimeError("central-directory length mismatch")
    return entries


def find_member(
    entries: list[tuple[str, tuple[int, int, int, int]]], wanted: str
) -> tuple[int, int, int, int]:
    matches = [entry for name, entry in entries if name == wanted]
    if len(matches) != 1:
        raise RuntimeError(f"expected exactly one member named {wanted!r}; found {len(matches)}")
    return matches[0]


def member_data_start(
    url: str, total: int, member: str, entry: tuple[int, int, int, int]
) -> int:
    local_offset, _, _, _ = entry
    local = fetch_range(url, local_offset, local_offset + 29, total)
    if local[:4] != LOCAL_SIGNATURE:
        raise RuntimeError(f"bad local-file signature for {member}")
    name_len, extra_len = struct.unpack_from("<HH", local, 26)
    return local_offset + 30 + name_len + extra_len


def first_lines(
    url: str,
    total: int,
    member: str,
    entry: tuple[int, int, int, int],
    line_count: int,
    max_compressed_bytes: int,
) -> list[str]:
    _, method, compressed_size, _ = entry
    data_start = member_data_start(url, total, member, entry)
    read_size = min(compressed_size, max_compressed_bytes)
    compressed = fetch_range(url, data_start, data_start + read_size - 1, total)
    if method == 0:
        text = compressed
    elif method == 8:
        text = zlib.decompressobj(-zlib.MAX_WBITS).decompress(compressed)
    else:
        raise RuntimeError(f"unsupported ZIP compression method {method}")
    lines = text.decode("utf-8-sig", errors="strict").splitlines()
    if len(lines) < line_count and read_size < compressed_size:
        raise RuntimeError(
            f"only decoded {len(lines)} line(s); raise --max-compressed-bytes"
        )
    return lines[:line_count]


def find_text_line(
    url: str,
    total: int,
    member: str,
    entry: tuple[int, int, int, int],
    wanted: str,
    max_compressed_bytes: int,
) -> str:
    _, method, compressed_size, _ = entry
    if method not in {0, 8}:
        raise RuntimeError(f"unsupported ZIP compression method {method}")
    data_start = member_data_start(url, total, member, entry)
    limit = min(compressed_size, max_compressed_bytes)
    decompressor = zlib.decompressobj(-zlib.MAX_WBITS) if method == 8 else None
    needle = wanted.encode("utf-8")
    carry = b""
    offset = 0
    chunk_size = 4 * 1024 * 1024
    while offset < limit:
        length = min(chunk_size, limit - offset)
        compressed = fetch_range(
            url, data_start + offset, data_start + offset + length - 1, total
        )
        decoded = decompressor.decompress(compressed) if decompressor else compressed
        lines = (carry + decoded).split(b"\n")
        carry = lines.pop()
        for line in lines:
            if needle in line:
                return line.rstrip(b"\r").decode("utf-8-sig", errors="strict")
        offset += length
    if offset == compressed_size and needle in carry:
        return carry.rstrip(b"\r").decode("utf-8-sig", errors="strict")
    raise RuntimeError(
        f"text {wanted!r} not found in the first {limit} compressed member bytes"
    )


def main() -> int:
    args = parse_args()
    if args.lines < 1 or args.max_compressed_bytes < 1:
        raise SystemExit("ERROR: --lines and --max-compressed-bytes must be positive")
    try:
        central = central_directory(args.url, args.expected_bytes)
        entries = member_entries(central)
        if args.list_members:
            for name, (_, method, compressed_size, uncompressed_size) in entries:
                print(f"{name}\t{method}\t{compressed_size}\t{uncompressed_size}")
            return 0
        entry = find_member(entries, args.member)
        if args.match:
            matched_line = find_text_line(
                args.url,
                args.expected_bytes,
                args.member,
                entry,
                args.match,
                args.max_compressed_bytes,
            )
            lines = []
        else:
            matched_line = None
            lines = first_lines(
                args.url,
                args.expected_bytes,
                args.member,
                entry,
                args.lines,
                args.max_compressed_bytes,
            )
    except Exception as error:  # keep this CLI's failure mode concise
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    local_offset, method, compressed_size, uncompressed_size = entry
    print(f"member\t{args.member}")
    print(f"compression_method\t{method}")
    print(f"compressed_bytes\t{compressed_size}")
    print(f"uncompressed_bytes\t{uncompressed_size}")
    print(f"local_header_offset\t{local_offset}")
    if matched_line is not None:
        print(f"matched_line\t{matched_line}")
    for index, line in enumerate(lines, start=1):
        print(f"line_{index}\t{line}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
