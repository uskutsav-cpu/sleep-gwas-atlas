#!/usr/bin/env python3
"""Materialize only h2-eligible rg jobs from the frozen replication manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=Path("discovery_extension/config/replication_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/replication_manifest.lock.json"))
    parser.add_argument("--h2", type=Path, default=Path("discovery_extension/results/replication/replication_source_h2.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/replication/replication_rg_jobs.tsv"))
    args = parser.parse_args()
    manifest = read_tsv(args.manifest)
    lock = json.loads(args.lock.read_text(encoding="utf-8"))
    if sha256(args.manifest) != lock.get("manifest_sha256"):
        raise SystemExit("ERROR: manifest differs from its pre-result lock")
    h2_rows = read_tsv(args.h2)
    h2_by_source = {row["replication_source_id"]: row for row in h2_rows}
    if len(h2_rows) != 13 or len(h2_by_source) != 13:
        raise SystemExit("ERROR: h2 table is not the frozen 13-source family")
    passing = {source_id for source_id, row in h2_by_source.items() if row["primary_status"] == "PASS"}
    testable = set(lock["testable_pair_ids_in_locked_order"])
    jobs: dict[str, list[dict[str, str]]] = {}
    for row in manifest:
        if row["pair_id"] in testable and row["replication_source_id"] in passing:
            jobs.setdefault(row["sleep_trait"], []).append(row)
    output = [
        {
            "sleep_trait": sleep, "source_ids": ",".join(row["replication_source_id"] for row in rows),
            "pair_ids": ",".join(row["pair_id"] for row in rows), "source_count": str(len(rows)),
        }
        for sleep, rows in jobs.items()
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sleep_trait", "source_ids", "pair_ids", "source_count"]
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    print(f"REPLICATION_RG_JOBS_OK jobs={len(output)} pairs={sum(int(row['source_count']) for row in output)} passing_sources={len(passing)} sha256={sha256(args.out)}")


if __name__ == "__main__":
    main()
