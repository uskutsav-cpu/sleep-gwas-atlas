#!/usr/bin/env python3
"""Freeze the replicated high-priority pair family for PLACO+ input curation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


FIELDS = [
    "pair_id", "sleep_trait", "extension_trait_id", "external_phenotype_name",
    "discovery_rg", "discovery_fdr", "replication_rg", "replication_p",
    "replication_class", "priority_tier", "selection_status",
    "merged_genomewide_path", "merged_genomewide_sha256", "input_variant_count",
    "eligible_variant_count", "z2_excluded_count",
    "input_scope", "build", "ancestry", "schema_status", "effect_allele_alignment_status",
    "maf_info_filter_status", "z_squared_filter_status", "ld_reference_id",
    "ld_reference_path", "ld_reference_sha256", "clumping_parameters",
    "placo_source_path", "placo_source_sha256", "results_accessed_before_lock",
    "curator", "curation_date",
]


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
    parser.add_argument("--priority", type=Path, default=Path("discovery_extension/results/prioritization/novel_hit_priority.tsv"))
    parser.add_argument("--replication", type=Path, default=Path("discovery_extension/results/replication/replication_results.tsv"))
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/pleiotropy_sources.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/pleiotropy/pleiotropy_input_queue.tsv"))
    parser.add_argument("--candidate-lock-out", type=Path, default=Path("discovery_extension/config/pleiotropy_candidate_family.lock.json"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/pleiotropy_input_queue.json"))
    args = parser.parse_args()

    priority = read_tsv(args.priority)
    replication = {row["pair_id"]: row for row in read_tsv(args.replication)}
    source_rows = read_tsv(args.sources)
    if len(source_rows) != 1:
        raise SystemExit("ERROR: exact PLACO+ source row missing")
    source = source_rows[0]
    candidates = [
        row for row in priority
        if row["priority_tier"] == "B"
        and row["pair_id"] in replication
        and replication[row["pair_id"]]["replication_class"] == "REPLICATED"
    ]
    if len({row["pair_id"] for row in candidates}) != len(candidates):
        raise SystemExit("ERROR: duplicate PLACO+ candidate pair")
    output: list[dict[str, str]] = []
    for row in candidates:
        repl = replication[row["pair_id"]]
        record = {field: "PENDING" for field in FIELDS}
        record.update({
            "pair_id": row["pair_id"], "sleep_trait": row["sleep_trait"],
            "extension_trait_id": row["extension_trait_id"],
            "external_phenotype_name": row["phenotype_name"],
            "discovery_rg": row["rg"], "discovery_fdr": row["extension_fdr"],
            "replication_rg": repl["replication_rg"], "replication_p": repl["replication_p"],
            "replication_class": repl["replication_class"], "priority_tier": row["priority_tier"],
            "selection_status": "SELECTED_REPLICATED_PRIORITY_PAIR",
            "input_scope": "PENDING_FULL_GENOME", "build": "GRCh37", "ancestry": "EUR",
            "placo_source_path": source["installed_source_path"],
            "placo_source_sha256": source["source_sha256"],
            "results_accessed_before_lock": "NO",
        })
        output.append(record)
    output.sort(key=lambda row: row["pair_id"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    family_lock = {
        "schema_version": "1.0.0", "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": "before_PLACO_input_curation_or_result_access",
        "results_accessed_before_lock": False, "pair_count": len(output),
        "pair_ids_in_locked_order": [row["pair_id"] for row in output],
        "priority_sha256": sha256(args.priority), "replication_sha256": sha256(args.replication),
        "queue_template_sha256": sha256(args.out),
        "selection_rule": "priority Tier B and independent replication class REPLICATED",
    }
    args.candidate_lock_out.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_lock_out.write_text(json.dumps(family_lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance = {
        "schema_version": "1.0.0", "result_free_template": True,
        "candidate_pair_count": len(output), "priority_sha256": sha256(args.priority),
        "replication_sha256": sha256(args.replication), "sources_sha256": sha256(args.sources),
        "output": str(args.out), "output_sha256": sha256(args.out),
        "warning": "The queue contains no PLACO+ result and cannot be used until all full-genome input fields are curated.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"PLEIOTROPY_QUEUE_OK pairs={len(output)} result_free=true")


if __name__ == "__main__":
    main()
