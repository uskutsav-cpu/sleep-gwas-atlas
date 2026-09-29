#!/usr/bin/env python3
"""Bind inherited dense-input MAF filtering to its upstream QC receipts."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def parse_qc(path: Path) -> tuple[dict[str, str], dict[str, tuple[int, int]]]:
    metadata: dict[str, str] = {}
    filters: dict[str, tuple[int, int]] = {}
    in_steps = False
    with path.open(encoding="utf-8") as stream:
        for raw in stream:
            line = raw.rstrip("\n")
            if line == "step\tdropped\tremaining":
                in_steps = True
                continue
            fields = line.split("\t")
            if in_steps and len(fields) == 3:
                try:
                    filters[fields[0]] = (int(fields[1]), int(fields[2]))
                    continue
                except ValueError:
                    pass
            if len(fields) == 2:
                metadata[fields[0]] = fields[1]
    return metadata, filters


def run(locked_path: Path, qc_master_path: Path, output_path: Path) -> int:
    locked = read_tsv(locked_path)
    qc_master = {r["trait_id"]: r for r in read_tsv(qc_master_path)}
    if len(locked) != 7 or len({r["trait_id"] for r in locked}) != 7:
        raise ValueError("Expected exactly seven unique locked dense inputs")
    rows = []
    for item in sorted(locked, key=lambda r: r["trait_id"]):
        trait = item["trait_id"]
        if item["status"] != "VERIFIED" or item["expected_sha256"] != item["observed_sha256"]:
            raise ValueError(f"Locked input is not checksum-verified: {trait}")
        dense = Path(item["dense_path"])
        qc_path = dense.with_name(f"{trait}.qc.txt")
        if not dense.is_file() or sha256(dense) != item["expected_sha256"]:
            raise ValueError(f"Dense input changed since lock: {trait}")
        if not qc_path.is_file():
            raise FileNotFoundError(qc_path)
        metadata, filters = parse_qc(qc_path)
        policy = "invalid FRQ or MAF <= 0.01 [CDG3]"
        if metadata.get("trait") != trait or policy not in filters:
            raise ValueError(f"Upstream QC lacks expected trait/policy: {trait}")
        if metadata.get("output_build") != "hg19":
            raise ValueError(f"Unexpected upstream output build for {trait}")
        upstream_rows = int(metadata["rows_out"])
        local_rows = int(qc_master[trait]["retained_rows"])
        if upstream_rows != local_rows:
            raise ValueError(f"Row-count lineage mismatch for {trait}: {upstream_rows} != {local_rows}")
        dropped, remaining = filters[policy]
        rows.append({
            "trait_id": trait,
            "dense_sha256": item["expected_sha256"],
            "upstream_qc_path": str(qc_path),
            "upstream_qc_sha256": sha256(qc_path),
            "upstream_input_sha256": metadata.get("infile_sha256", ""),
            "upstream_filter_record": policy,
            "upstream_frequency_rows_dropped": str(dropped),
            "upstream_rows_after_frequency_filter": str(remaining),
            "upstream_rows_out": str(upstream_rows),
            "current_dense_rows": str(local_rows),
            "row_count_match": "YES",
            "frequency_orientation": "NOT_ESTABLISHED_BY_QC_LOG",
            "info_policy": "NOT_INFERRED_FROM_FREQUENCY_AUDIT",
            "status": "PASS_INHERITED_MAF_FILTER_LINEAGE",
        })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with output_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locked-inputs", type=Path, default=Path("brain6/manifests/locked_dense_input_audit.tsv"))
    parser.add_argument("--qc-master", type=Path, default=Path("brain6/qc/dense_source_qc_master.tsv"))
    parser.add_argument("--output", type=Path, default=Path("brain6/qc/upstream_frequency_filter_audit.tsv"))
    args = parser.parse_args()
    print(f"audited_traits={run(args.locked_inputs, args.qc_master, args.output)} output={args.output}")


if __name__ == "__main__":
    main()
