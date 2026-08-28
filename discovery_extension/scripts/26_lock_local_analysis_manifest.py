#!/usr/bin/env python3
"""Freeze the exact pair-by-locus local-analysis family before result access."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


SELECTED = {"SELECTED_PRIORITY_DISCOVERY", "SELECTED_GLOBAL_NULL_SECONDARY"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_locked_file(row: dict[str, str], path_field: str, checksum_field: str) -> None:
    path = Path(row[path_field])
    expected = row[checksum_field]
    if not path.is_file():
        raise SystemExit(f"ERROR: selected pair {row['pair_id']} lacks {path_field}: {path}")
    if len(expected) != 64 or sha256(path) != expected:
        raise SystemExit(f"ERROR: selected pair {row['pair_id']} has invalid {checksum_field}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/local/local_analysis_queue.tsv"))
    parser.add_argument("--readiness", type=Path, default=Path("discovery_extension/results/local/local_architecture_readiness.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/config/local_analysis_manifest.tsv"))
    parser.add_argument("--lock", type=Path, default=Path("discovery_extension/config/local_analysis_manifest.lock.json"))
    args = parser.parse_args()

    queue = read_tsv(args.queue)
    readiness = {row["method_id"]: row for row in read_tsv(args.readiness)}
    selected = [row for row in queue if row["selection_status"] in SELECTED]
    if not selected:
        raise SystemExit("ERROR: no local-analysis pair has a completed pre-result selection")
    if len({row["pair_id"] for row in selected}) != len(selected):
        raise SystemExit("ERROR: duplicate selected pair_id")
    for row in selected:
        if row["result_access_status"] != "NOT_ACCESSED":
            raise SystemExit(f"ERROR: local results were accessed before lock: {row['pair_id']}")
        if row["global_analysis_status"] != "PRIMARY_EXTENSION_RG_COMPLETE":
            raise SystemExit(f"ERROR: selected pair lacks clean global rg: {row['pair_id']}")
        if row["selection_status"] == "SELECTED_GLOBAL_NULL_SECONDARY":
            if float(row["global_rg_fdr"]) < 0.05 or row["selection_stratum"] != "GLOBAL_NULL_SECONDARY":
                raise SystemExit(f"ERROR: secondary pair is not globally null: {row['pair_id']}")
            relevance = row["biological_relevance"] in {"HIGH", "MODERATE"}
            overlap_evidence = row["polygenic_overlap_evidence"].strip() not in {"", "PENDING", "NONE"}
            if not relevance and not overlap_evidence:
                raise SystemExit(f"ERROR: globally-null selection lacks relevance/overlap rationale: {row['pair_id']}")
        if row["dense_summary_statistics_available"] != "YES":
            raise SystemExit(f"ERROR: selected pair lacks dense summary statistics: {row['pair_id']}")
        try:
            lava_seed = int(row["lava_simulation_random_seed"])
        except ValueError as exc:
            raise SystemExit(f"ERROR: invalid LAVA simulation seed for {row['pair_id']}") from exc
        if not 1 <= lava_seed <= 2**31 - 1:
            raise SystemExit(f"ERROR: LAVA simulation seed out of range for {row['pair_id']}")
        for field in ("curator", "curation_date", "secondary_selection_rationale", "hdl_feasibility_rationale"):
            if not row[field].strip() or row[field] == "PENDING":
                raise SystemExit(f"ERROR: selected pair lacks {field}: {row['pair_id']}")
        try:
            locus_count = int(row["locus_count"])
        except ValueError as exc:
            raise SystemExit(f"ERROR: invalid locus_count for {row['pair_id']}") from exc
        if locus_count < 1:
            raise SystemExit(f"ERROR: locus_count must be positive: {row['pair_id']}")
        for path_field, checksum_field in (
            ("trait_a_dense_path", "trait_a_dense_sha256"),
            ("trait_b_dense_path", "trait_b_dense_sha256"),
            ("locus_definition_path", "locus_definition_sha256"),
        ):
            require_locked_file(row, path_field, checksum_field)
        methods = row["planned_methods"].split(";")
        if not methods or methods[0] != "LAVA_PRIMARY" or not set(methods) <= {"LAVA_PRIMARY", "HDL_L_ROBUSTNESS"}:
            raise SystemExit(f"ERROR: invalid planned_methods for {row['pair_id']}")
        for method in methods:
            if method not in readiness or readiness[method]["readiness"] != "READY":
                raise SystemExit(f"ERROR: {method} dependencies are not READY for {row['pair_id']}")

    fields = list(queue[0])
    selected.sort(key=lambda row: (row["selection_stratum"], row["pair_id"]))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected)
    pair_locus_tests = sum(int(row["locus_count"]) for row in selected)
    local_h2_test_count = 2 * pair_locus_tests
    local_h2_alpha = 0.05 / local_h2_test_count
    lock = {
        "schema_version": "1.0.0", "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "results_accessed_before_lock": False, "manifest": str(args.out),
        "manifest_sha256": sha256(args.out), "pair_ids": [row["pair_id"] for row in selected],
        "pair_count": len(selected), "lava_pair_locus_family_size": pair_locus_tests,
        "lava_local_h2_test_count": local_h2_test_count,
        "lava_local_h2_bonferroni_alpha": local_h2_alpha,
        "primary_method": "LAVA_PRIMARY", "primary_local_rg_fdr_method": "Benjamini-Hochberg",
        "primary_local_rg_fdr_threshold": 0.05,
        "primary_family_definition": "all clean estimable LAVA_PRIMARY bivariate tests in this exact locked pair-by-locus family",
        "hdl_l_role": "robustness only; excluded from the LAVA primary FDR family",
        "globally_null_pair_count": sum(row["selection_stratum"] == "GLOBAL_NULL_SECONDARY" for row in selected),
        "readiness_sha256": sha256(args.readiness), "queue_sha256": sha256(args.queue),
        "warning": "Locking a family authorizes no claim; failed, inestimable, and QC-failed loci must remain in normalized results.",
    }
    args.lock.parent.mkdir(parents=True, exist_ok=True)
    args.lock.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"LOCAL_ANALYSIS_MANIFEST_LOCKED pairs={len(selected)} lava_pair_locus_family={pair_locus_tests} result_free=true")


if __name__ == "__main__":
    main()
