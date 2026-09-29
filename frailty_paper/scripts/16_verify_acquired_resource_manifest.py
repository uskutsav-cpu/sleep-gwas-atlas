#!/usr/bin/env python3
"""Read-only verification of the acquired-resource size and SHA-256 manifest."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

REQUIRED = {"resource_id", "file", "bytes", "sha256"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest(root: Path, manifest: Path, external_root: Path | None = None) -> tuple[int, list[str]]:
    errors: list[str] = []
    try:
        with manifest.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = set(reader.fieldnames or [])
            missing_fields = REQUIRED - fields
            if missing_fields:
                return 0, [f"manifest missing required columns: {', '.join(sorted(missing_fields))}"]
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        return 0, [f"cannot read manifest {manifest}: {exc}"]

    seen: set[str] = set()
    for line, row in enumerate(rows, start=2):
        resource_id = row.get("resource_id", "").strip()
        relpath = row.get("file", "").strip()
        if not resource_id:
            errors.append(f"line {line}: empty resource_id")
        elif resource_id in seen:
            errors.append(f"line {line}: duplicate resource_id {resource_id}")
        seen.add(resource_id)
        relative = Path(relpath)
        if not relpath or relative.is_absolute() or ".." in relative.parts:
            errors.append(f"line {line}: file path must be a nonempty repo-relative path: {relpath!r}")
            continue
        path = (root / relative).resolve()
        allowed_roots = [root.resolve()]
        if external_root is not None:
            allowed_roots.append(external_root.resolve())
        if not any(path == allowed or allowed in path.parents for allowed in allowed_roots):
            errors.append(f"line {line}: file path escapes repository root: {relpath}")
            continue
        if not path.is_file():
            errors.append(f"line {line}: missing file: {relpath}")
            continue
        try:
            expected_bytes = int(row.get("bytes", ""))
        except ValueError:
            errors.append(f"line {line}: invalid byte count for {relpath}")
            continue
        if expected_bytes < 0:
            errors.append(f"line {line}: negative byte count for {relpath}")
            continue
        if path.stat().st_size != expected_bytes:
            errors.append(f"line {line}: byte count mismatch for {relpath}")
            continue
        expected_hash = row.get("sha256", "").strip().lower()
        if len(expected_hash) != 64 or any(c not in "0123456789abcdef" for c in expected_hash):
            errors.append(f"line {line}: invalid SHA-256 for {relpath}")
            continue
        actual_hash = sha256(path)
        if actual_hash != expected_hash:
            errors.append(f"line {line}: SHA-256 mismatch for {relpath}")
    return len(rows), errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", default=".")
    parser.add_argument("--manifest", default="frailty_paper/manifests/all_acquired_resources.tsv")
    parser.add_argument("--external-storage-root", help="Explicitly allow manifest paths resolving under this external data root")
    args = parser.parse_args()
    root = Path(args.repo).resolve()
    manifest = Path(args.manifest)
    if not manifest.is_absolute():
        manifest = root / manifest
    external_root = Path(args.external_storage_root) if args.external_storage_root else None
    count, errors = verify_manifest(root, manifest, external_root)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        print(f"RESOURCE_MANIFEST_FAILED rows={count} errors={len(errors)}")
        return 1
    print(f"RESOURCE_MANIFEST_OK rows={count} files_verified={count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
