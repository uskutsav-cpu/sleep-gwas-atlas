"""Resumable public downloads with exact byte count + independently supplied hash."""
from __future__ import annotations
import os
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from .io import check_hash, ensure_free, exclusive_lock, read_json, require, sha256, write_json


def download(url: str, destination: str | Path, expected_hash: str, expected_bytes: int,
             *, allow_local_http=False, timeout=60) -> Path:
    parsed = urllib.parse.urlsplit(url)
    local = parsed.hostname in {"127.0.0.1", "localhost", "::1"}
    require(parsed.scheme == "https" or (allow_local_http and local and parsed.scheme == "http"),
            "Use HTTPS; HTTP is restricted to loopback tests")
    require(not parsed.username and not parsed.password, "Do not put credentials in source URLs")
    require(expected_bytes > 0, "Exact positive byte count required")
    require(bool(re.fullmatch(r"[a-fA-F0-9]{64}", expected_hash)), "Expected SHA-256 required")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    require(not destination.is_symlink(), "Refusing destination symlink")
    if destination.exists():
        require(destination.stat().st_size == expected_bytes, "Existing destination has different length")
        check_hash(destination, expected_hash)
        return destination
    partial = destination.with_name(destination.name + ".partial")
    stamp = destination.with_name(destination.name + ".download.json")
    identity = {"url": url, "sha256": expected_hash.lower(), "bytes": expected_bytes}
    with exclusive_lock(destination.with_name(destination.name + ".download.lock")):
        require(not partial.is_symlink() and not stamp.is_symlink(), "Refusing partial-file symlink")
        if stamp.exists():
            require(read_json(stamp) == identity, "Partial download belongs to different source")
        else:
            require(not partial.exists(), "Unattributed partial download; inspect it before reuse")
            write_json(stamp, identity)
        offset = partial.stat().st_size if partial.exists() else 0
        require(offset <= expected_bytes, "Partial download longer than expected")
        ensure_free(destination.parent, expected_bytes - offset + 1024**2)
        if offset < expected_bytes:
            headers = {"Accept-Encoding": "identity", "User-Agent": "sleep-gwas-brain6/0.1"}
            if offset:
                headers["Range"] = f"bytes={offset}-"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as response:
                final_url = urllib.parse.urlsplit(response.geturl())
                require(final_url.scheme == "https" or (allow_local_http and final_url.hostname in {"127.0.0.1", "localhost", "::1"}),
                        "HTTPS redirect downgrade forbidden")
                if offset:
                    require(response.status == 206, "Server ignored resume Range; partial left intact")
                    match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", response.headers.get("Content-Range", ""))
                    require(match is not None and int(match[1]) == offset and int(match[3]) == expected_bytes
                            and int(match[2]) == expected_bytes - 1, "Invalid Content-Range")
                else:
                    require(response.status == 200, "Unexpected download HTTP status")
                length = response.headers.get("Content-Length")
                if length is not None:
                    require(int(length) == expected_bytes - offset, "Content-Length differs from locked source")
                mode = "ab" if offset else "wb"
                with partial.open(mode) as output:
                    for block in iter(lambda: response.read(4 * 1024 * 1024), b""):
                        require(output.tell() + len(block) <= expected_bytes, "Download exceeds registered length")
                        output.write(block)
                    output.flush(); os.fsync(output.fileno())
        require(partial.stat().st_size == expected_bytes, "Incomplete download; resumable partial preserved")
        check_hash(partial, expected_hash)
        os.link(partial, destination)
        partial.unlink()
    return destination
