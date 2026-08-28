#!/usr/bin/env python3
"""Lock independently curated replication sources before results are accessed."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_COMMON = {
    "results_accessed_before_lock", "source_curation_status",
    "replication_search_databases", "replication_search_queries", "replication_search_date",
    "replication_search_evidence", "unavailable_reason",
}
REQUIRED_TESTABLE = {
    "replication_source_id", "replication_study_accession", "replication_publication",
    "replication_PMID", "replication_DOI", "replication_source_url", "replication_checksum",
    "replication_local_path",
    "replication_phenotype_definition", "phenotype_match_status", "ancestry", "build",
    "sample_size", "cases", "controls", "discovery_cohort_relation",
    "participant_overlap_status", "participant_overlap_evidence",
    "source_identity_status", "schema_status", "effect_allele_status", "full_resolution_availability",
}


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--queue", type=Path,
        default=Path("discovery_extension/results/replication/replication_source_queue.tsv"),
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("discovery_extension/config/replication_manifest.tsv"),
    )
    parser.add_argument(
        "--candidate-lock", type=Path,
        default=Path("discovery_extension/config/replication_candidate_family.lock.json"),
    )
    parser.add_argument(
        "--lock", type=Path,
        default=Path("discovery_extension/config/replication_manifest.lock.json"),
    )
    args = parser.parse_args()
    fields, rows = read_tsv(args.queue)
    candidate_lock = json.loads(args.candidate_lock.read_text(encoding="utf-8"))
    if not rows:
        raise SystemExit("ERROR: no replication candidates can be locked")
    missing_columns = (REQUIRED_COMMON | REQUIRED_TESTABLE) - set(fields)
    if missing_columns:
        raise SystemExit(f"ERROR: replication queue lacks columns: {sorted(missing_columns)}")
    if len({row["pair_id"] for row in rows}) != len(rows):
        raise SystemExit("ERROR: duplicate replication pair_id")
    if [row["pair_id"] for row in rows] != candidate_lock.get("pair_ids_in_locked_order"):
        raise SystemExit("ERROR: curated replication queue membership/order differs from candidate-family lock")
    if len(rows) != candidate_lock.get("pair_count") or candidate_lock.get("results_accessed_before_lock") is not False:
        raise SystemExit("ERROR: replication candidate-family lock is invalid")
    family_sizes = {int(row["replication_family_size"]) for row in rows}
    alphas = {float(row["replication_alpha"]) for row in rows}
    if family_sizes != {len(rows)} or len(alphas) != 1 or abs(next(iter(alphas)) - 0.05 / len(rows)) > 1e-12:
        raise SystemExit("ERROR: replication family size/alpha is not the frozen Bonferroni family")
    testable_rows: list[dict[str, str]] = []
    unavailable_rows: list[dict[str, str]] = []
    for row in rows:
        if row["results_accessed_before_lock"] != "NO":
            raise SystemExit(f"ERROR: replication results were accessed before source lock: {row['pair_id']}")
        for field in ("replication_search_databases", "replication_search_queries", "replication_search_date", "replication_search_evidence"):
            if not row[field].strip() or row[field].startswith("PENDING"):
                raise SystemExit(f"ERROR: {row['pair_id']} lacks completed independent-source search field {field}")
        if row["source_curation_status"] == "NO_INDEPENDENT_DATASET_COMPLETE_BEFORE_RESULTS":
            if not row["unavailable_reason"].strip() or row["unavailable_reason"].startswith("PENDING"):
                raise SystemExit(f"ERROR: unavailable replication candidate lacks reason: {row['pair_id']}")
            unavailable_rows.append(row)
            continue
        for field in REQUIRED_TESTABLE:
            if not row[field].strip() or row[field].startswith("PENDING"):
                raise SystemExit(f"ERROR: {row['pair_id']} lacks curated {field}")
        if row["phenotype_match_status"] not in {"EXACT", "COMPARABLE_WITH_DOCUMENTED_DIFFERENCES"}:
            raise SystemExit(f"ERROR: unsupported phenotype match: {row['pair_id']}")
        if row["ancestry"] != "EUR":
            raise SystemExit(f"ERROR: primary replication manifest requires EUR ancestry: {row['pair_id']}")
        if row["participant_overlap_status"] not in {"NON_OVERLAPPING_CONFIRMED", "OVERLAP_QUANTIFIED_NEGLIGIBLE"}:
            raise SystemExit(f"ERROR: replication independence not established: {row['pair_id']}")
        if row["discovery_cohort_relation"] not in {"NON_UKB", "INDEPENDENT_COHORT"}:
            raise SystemExit(f"ERROR: same-cohort/internal split is not independent replication: {row['pair_id']}")
        if row["source_identity_status"] != "VERIFIED" or row["schema_status"] != "VERIFIED":
            raise SystemExit(f"ERROR: replication source/schema is not verified: {row['pair_id']}")
        if row["effect_allele_status"] != "UNAMBIGUOUS" or row["full_resolution_availability"] != "YES":
            raise SystemExit(f"ERROR: replication effect/full-resolution contract failed: {row['pair_id']}")
        if row["source_curation_status"] != "COMPLETE_BEFORE_RESULTS":
            raise SystemExit(f"ERROR: replication curation is incomplete: {row['pair_id']}")
        if not row["replication_source_url"].startswith("https://") or not row["replication_checksum"].startswith("sha256:") or len(row["replication_checksum"].removeprefix("sha256:")) != 64:
            raise SystemExit(f"ERROR: replication source URL/checksum is invalid: {row['pair_id']}")
        local_path = Path(row["replication_local_path"])
        if not local_path.is_file() or sha256(local_path) != row["replication_checksum"].removeprefix("sha256:"):
            raise SystemExit(f"ERROR: replication source is missing or differs from its SHA256: {row['pair_id']}")
        testable_rows.append(row)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    lock = {
        "schema_version": "1.0.0", "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": "before_replication_result_access", "pair_count": len(rows),
        "pair_ids_in_locked_order": [row["pair_id"] for row in rows],
        "testable_pair_ids_in_locked_order": [row["pair_id"] for row in testable_rows],
        "unavailable_pair_ids_in_locked_order": [row["pair_id"] for row in unavailable_rows],
        "bonferroni_alpha": next(iter(alphas)), "manifest_sha256": sha256(args.out),
        "source_queue_sha256": sha256(args.queue),
        "candidate_family_lock_sha256": sha256(args.candidate_lock),
        "independence_required": True, "same_cohort_internal_split_allowed": False,
    }
    args.lock.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    print(f"REPLICATION_MANIFEST_LOCKED pairs={len(rows)} testable={len(testable_rows)} unavailable={len(unavailable_rows)} alpha={lock['bonferroni_alpha']:.12g} sha256={lock['manifest_sha256']}")


if __name__ == "__main__":
    main()
