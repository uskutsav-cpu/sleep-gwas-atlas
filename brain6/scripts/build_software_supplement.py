#!/usr/bin/env python3
"""Publish the source-of-truth software ledger as supplementary Table S17."""
from __future__ import annotations

import hashlib
import json
import os
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "brain6/manifests/software_versions.tsv"
OUTPUT = ROOT / "brain6/results/supplement/table_S17_software_resources.tsv"
PROVENANCE = ROOT / "brain6/results/supplement/table_S17_software_resources.provenance.json"


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def replace_if_changed(path: Path, data: bytes) -> None:
    if path.is_file() and path.read_bytes() == data:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_name, path)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)


def main() -> None:
    source_bytes = SOURCE.read_bytes()
    digest = sha256(source_bytes)
    lines = source_bytes.decode("utf-8").splitlines()
    if not lines or lines[0] != "component\tversion\tstatus_or_evidence":
        raise ValueError("Software ledger has an unexpected TSV schema")
    rows = [line for line in lines[1:] if line]
    if not rows or any(len(line.split("\t")) != 3 for line in rows):
        raise ValueError("Software ledger contains empty or malformed rows")

    replace_if_changed(OUTPUT, source_bytes)
    record = {
        "schema_version": 1,
        "status": "PASS",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": digest,
        "output": str(OUTPUT.relative_to(ROOT)),
        "output_sha256": sha256(OUTPUT.read_bytes()),
        "data_rows": len(rows),
    }
    serialized = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode("utf-8")
    replace_if_changed(PROVENANCE, serialized)
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    main()
