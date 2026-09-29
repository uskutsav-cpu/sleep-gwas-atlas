#!/usr/bin/env python3
"""Fail-closed structural audit of the acquired-resource metadata manifest."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

REQUIRED_COLUMNS = (
    "resource_id", "trait", "resource_type", "study", "publication",
    "accession", "source", "download_date", "file", "bytes", "sha256",
    "genome_build", "ancestry", "sample_size", "cases", "controls",
    "effect_type", "cohorts", "ukb_overlap", "finngen_overlap", "license",
    "notes",
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def audit_manifest(path: Path) -> tuple[int, list[str]]:
    errors: list[str] = []
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames or []
            if len(columns) != len(set(columns)):
                return 0, ["duplicate column names"]
            missing_columns = [name for name in REQUIRED_COLUMNS if name not in columns]
            if missing_columns:
                return 0, [f"missing required columns: {', '.join(missing_columns)}"]
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        return 0, [f"cannot read manifest {path}: {exc}"]

    if not rows:
        return 0, ["manifest has no resource rows"]

    seen: set[str] = set()
    for line, row in enumerate(rows, start=2):
        resource_id = (row.get("resource_id") or "").strip()
        if not resource_id:
            errors.append(f"line {line}: empty resource_id")
        elif resource_id in seen:
            errors.append(f"line {line}: duplicate resource_id {resource_id}")
        seen.add(resource_id)

        for column in REQUIRED_COLUMNS:
            if not (row.get(column) or "").strip():
                errors.append(f"line {line}: empty required field {column}")

        try:
            size = int((row.get("bytes") or "").strip())
        except ValueError:
            errors.append(f"line {line}: invalid byte count")
        else:
            if size <= 0:
                errors.append(f"line {line}: byte count must be positive")

        checksum = (row.get("sha256") or "").strip().lower()
        if not SHA256_RE.fullmatch(checksum):
            errors.append(f"line {line}: invalid SHA-256")

    return len(rows), errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("frailty_paper/manifests/all_acquired_resources.tsv"),
    )
    args = parser.parse_args()
    count, errors = audit_manifest(args.manifest)
    if errors:
        for error in errors:
            print(f"ERROR: {error}")
        print(f"MANIFEST_METADATA_FAILED rows={count} errors={len(errors)}")
        return 1
    print(f"MANIFEST_METADATA_OK rows={count} required_columns={len(REQUIRED_COLUMNS)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
