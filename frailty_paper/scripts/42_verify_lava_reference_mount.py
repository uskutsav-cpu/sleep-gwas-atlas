#!/usr/bin/env python3
"""Verify the pinned LAVA reference when its repository path is an external mount.

The repository-wide LAVA verifier intentionally rejects symlinks.  Frailty
inputs are stored on the designated external volume, so this read-only helper
accepts the physical reference directory explicitly and verifies it against
the already-pinned atlas provenance before any frailty LAVA use.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
import sys
from pathlib import Path


PINNED_REFERENCE_SHA256 = "35ab371935b9c260dab21be248f4907b693565787ef9f62334e7d8030204270e"
PINNED_EXTRACTED_MANIFEST_SHA256 = "5e250575a5afba843d2bed12cbdc4cc87071267d3096ca1f4b25e50da9dba579"
REFERENCE_ID = "LAVA UK Biobank LD v1.1"
EXPECTED_FILES = 44


class VerificationError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    before_path = os.lstat(path)
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
    after_path = os.lstat(path)
    identities = {
        (item.st_dev, item.st_ino, item.st_size, item.st_mtime_ns, item.st_ctime_ns)
        for item in (before_path, before, after, after_path)
    }
    if len(identities) != 1:
        raise VerificationError(f"file changed while hashing: {path}")
    return digest.hexdigest()


def require_regular_file(path: Path, label: str) -> None:
    try:
        info = os.lstat(path)
    except OSError as error:
        raise VerificationError(f"cannot inspect {label} {path}: {error}") from error
    if not stat.S_ISREG(info.st_mode):
        raise VerificationError(f"{label} is not a regular non-symlink file: {path}")
    if info.st_size <= 0:
        raise VerificationError(f"{label} is empty: {path}")


def verify(root: Path) -> tuple[int, int]:
    if root.is_symlink() or not root.is_dir():
        raise VerificationError(f"reference root must be a physical directory: {root}")

    provenance_path = root / "reference.provenance.json"
    require_regular_file(provenance_path, "reference provenance")
    provenance_sha = sha256_file(provenance_path)
    if provenance_sha != PINNED_REFERENCE_SHA256:
        raise VerificationError(
            f"reference provenance SHA-256 mismatch: {provenance_sha}"
        )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    if provenance.get("reference") != REFERENCE_ID:
        raise VerificationError("unexpected LAVA reference identity")
    if provenance.get("extracted_file_count") != EXPECTED_FILES:
        raise VerificationError("unexpected extracted reference file count")

    manifest_path = root / "extracted_manifest.tsv"
    require_regular_file(manifest_path, "extracted manifest")
    manifest_sha = sha256_file(manifest_path)
    if manifest_sha != PINNED_EXTRACTED_MANIFEST_SHA256:
        raise VerificationError(f"extracted manifest SHA-256 mismatch: {manifest_sha}")
    if provenance.get("extracted_manifest_sha256") != manifest_sha:
        raise VerificationError("provenance and extracted manifest hashes disagree")

    entries = provenance.get("extracted_files")
    if not isinstance(entries, list) or len(entries) != EXPECTED_FILES:
        raise VerificationError("provenance does not list exactly 44 payloads")

    observed_bytes = 0
    seen: set[str] = set()
    observed_types = {".bcor": 0, ".info": 0}
    for entry in entries:
        rel = Path(str(entry.get("path", "")))
        prefix = Path("ref/lava/ukb_v1.1")
        if rel.is_absolute() or ".." in rel.parts or rel.parts[:3] != prefix.parts:
            raise VerificationError(f"unsafe reference payload path: {rel}")
        local_rel = Path(*rel.parts[3:])
        if not local_rel.parts or local_rel.as_posix() in seen:
            raise VerificationError(f"empty or duplicate reference path: {rel}")
        seen.add(local_rel.as_posix())
        path = root / local_rel
        require_regular_file(path, "reference payload")
        info = os.lstat(path)
        expected_bytes = int(entry["bytes"])
        if info.st_size != expected_bytes:
            raise VerificationError(f"byte-size mismatch for {path}")
        observed_sha = sha256_file(path)
        if observed_sha != entry["sha256"]:
            raise VerificationError(f"SHA-256 mismatch for {path}")
        if path.suffix not in observed_types:
            raise VerificationError(f"unexpected reference payload type: {path}")
        observed_types[path.suffix] += 1
        observed_bytes += info.st_size

    if observed_types != {".bcor": 22, ".info": 22}:
        raise VerificationError(f"unexpected reference payload types: {observed_types}")
    if observed_bytes != provenance.get("extracted_bytes"):
        raise VerificationError("extracted payload byte total disagrees with provenance")
    return len(entries), observed_bytes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    args = parser.parse_args()
    try:
        count, total_bytes = verify(args.reference_root)
    except (KeyError, OSError, ValueError, json.JSONDecodeError, VerificationError) as error:
        print(f"FAIL: {error}", file=sys.stderr)
        return 1
    print(
        f"PASS: {REFERENCE_ID}; provenance and manifest pinned; "
        f"{count} payloads (22 .bcor, 22 .info), {total_bytes} bytes; all sizes and SHA-256 values match."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
