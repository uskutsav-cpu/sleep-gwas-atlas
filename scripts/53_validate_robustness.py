#!/usr/bin/env python3
"""Validate the complete nine-family robustness matrix for major conclusions."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


SHA256 = re.compile(r"^[0-9a-f]{64}$")
FORBIDDEN = ("SYNTHETIC", "PLACEHOLDER", "FAKE_RESULT", "SMOKE_TEST")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--input", default="results/tables/robustness_summary.tsv")
    parser.add_argument("--policy", default="config/downstream_analysis_policy.json")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    policy = json.loads((root / args.policy).read_text(encoding="utf-8"))["robustness"]
    path = root / args.input
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing real non-empty artifact: {args.input}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        rows = list(reader)
    if fields != policy["table_fields"]:
        fail("robustness table header differs from the frozen policy")
    if not rows:
        fail("robustness table has no major-conclusion rows")
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for line, row in enumerate(rows, start=2):
        identity = row["conclusion_id"]
        if not identity or not row["conclusion"] or not row["provenance_id"]:
            fail(f"missing conclusion/provenance identity at line {line}")
        if any(marker in "\t".join(row.values()).upper() for marker in FORBIDDEN):
            fail(f"forbidden synthetic/placeholder marker at line {line}")
        if row["robustness_family"] not in policy["required_families"]:
            fail(f"unknown robustness family for {identity}")
        if row["applicable"] not in {"TRUE", "FALSE"}:
            fail(f"invalid applicability for {identity}/{row['robustness_family']}")
        if row["material_contradiction"] not in {"TRUE", "FALSE"}:
            fail(f"invalid contradiction flag for {identity}/{row['robustness_family']}")
        for field in ("direction_concordant", "significance_concordant"):
            if row[field] not in {"TRUE", "FALSE", "NA"}:
                fail(f"invalid {field} for {identity}/{row['robustness_family']}")
        if row["material_contradiction"] == "TRUE" and row["resolution"] in {"", "NA", "UNRESOLVED"}:
            fail(f"unresolved material contradiction for {identity}/{row['robustness_family']}")
        if row["applicable"] == "TRUE":
            relative = Path(row["evidence_path"])
            if relative.is_absolute() or ".." in relative.parts or not relative.parts or relative.parts[0] != "results":
                fail(f"unsafe/non-release evidence path for {identity}/{row['robustness_family']}")
            evidence = root / relative
            if not evidence.is_file() or evidence.stat().st_size == 0:
                fail(f"missing robustness evidence for {identity}/{row['robustness_family']}")
            if not SHA256.fullmatch(row["evidence_sha256"]) or sha256(evidence) != row["evidence_sha256"]:
                fail(f"robustness evidence checksum mismatch for {identity}/{row['robustness_family']}")
            if row["primary_result"] in {"", "NA"} or row["sensitivity_result"] in {"", "NA"}:
                fail(f"missing applicable result for {identity}/{row['robustness_family']}")
            if row["direction_concordant"] == "NA" or row["significance_concordant"] == "NA":
                fail(f"missing concordance assessment for {identity}/{row['robustness_family']}")
        elif row["resolution"] in {"", "NA"}:
            fail(f"non-applicability is unexplained for {identity}/{row['robustness_family']}")
        grouped[identity].append(row)
    required = policy["required_families"]
    for identity, values in grouped.items():
        observed = [row["robustness_family"] for row in values]
        if len(values) != len(required) or set(observed) != set(required) or len(observed) != len(set(observed)):
            fail(f"{identity} does not cover each locked robustness family exactly once")
        conclusions = {row["conclusion"] for row in values}
        if len(conclusions) != 1:
            fail(f"conclusion text drifted within {identity}")
    if not args.quiet:
        print(f"ROBUSTNESS_OK conclusions={len(grouped)} rows={len(rows)} families={len(required)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
