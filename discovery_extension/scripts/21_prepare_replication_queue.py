#!/usr/bin/env python3
"""Create a result-free source-curation queue for Tier A replication pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


SOURCE_FIELDS = [
    "pair_id", "sleep_trait", "extension_trait_id", "external_phenotype_name",
    "discovery_rg", "discovery_se", "discovery_fdr", "pair_novelty_class", "pair_novelty_strength",
    "discovery_sleep_source_id", "discovery_sleep_PMID", "discovery_sleep_DOI",
    "discovery_external_source", "discovery_external_study_accession",
    "discovery_external_PMID", "discovery_external_DOI",
    "replication_source_id", "replication_study_accession", "replication_publication",
    "replication_PMID", "replication_DOI", "replication_source_url", "replication_checksum",
    "replication_local_path",
    "replication_phenotype_definition", "phenotype_match_status", "ancestry", "build",
    "sample_size", "cases", "controls", "discovery_cohort_relation",
    "participant_overlap_status", "participant_overlap_evidence",
    "source_identity_status", "schema_status", "effect_allele_status",
    "full_resolution_availability", "results_accessed_before_lock",
    "replication_search_databases", "replication_search_queries", "replication_search_date",
    "replication_search_evidence", "unavailable_reason",
    "source_curation_status", "replication_family_size", "replication_alpha",
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
    parser.add_argument(
        "--priority", type=Path,
        default=Path("discovery_extension/results/prioritization/novel_hit_priority.tsv"),
    )
    parser.add_argument(
        "--panel", type=Path, default=Path("discovery_extension/config/candidate_traits.tsv")
    )
    parser.add_argument("--core", type=Path, default=Path("config/analysis_panel.tsv"))
    parser.add_argument(
        "--out", type=Path,
        default=Path("discovery_extension/results/replication/replication_source_queue.tsv"),
    )
    parser.add_argument(
        "--provenance-out", type=Path,
        default=Path("discovery_extension/provenance/replication_source_queue.json"),
    )
    parser.add_argument(
        "--candidate-lock-out", type=Path,
        default=Path("discovery_extension/config/replication_candidate_family.lock.json"),
    )
    args = parser.parse_args()
    priority = read_tsv(args.priority)
    panel = {row["extension_trait_id"]: row for row in read_tsv(args.panel)}
    core = {row["trait_id"]: row for row in read_tsv(args.core)}
    candidates = [
        row for row in priority
        if row["priority_tier"] in {"A", "B"} and row["replication_class"] != "REPLICATED"
    ]
    if len({row["pair_id"] for row in candidates}) != len(candidates):
        raise SystemExit("ERROR: duplicate priority pair selected for replication")
    family_size = len(candidates)
    alpha = 0.05 / family_size if family_size else 0.0
    output: list[dict[str, str]] = []
    for row in candidates:
        trait = panel[row["extension_trait_id"]]
        sleep = core[row["sleep_trait"]]
        record = {field: "PENDING" for field in SOURCE_FIELDS}
        record.update({
            "pair_id": row["pair_id"], "sleep_trait": row["sleep_trait"],
            "extension_trait_id": row["extension_trait_id"],
            "external_phenotype_name": row["phenotype_name"],
            "discovery_rg": row["rg"], "discovery_se": row["se"], "discovery_fdr": row["extension_fdr"],
            "pair_novelty_class": row["pair_novelty_class"],
            "pair_novelty_strength": row["pair_novelty_strength"],
            "discovery_sleep_source_id": sleep["source_id"],
            "discovery_sleep_PMID": sleep["pmid"], "discovery_sleep_DOI": sleep["doi"],
            "discovery_external_source": trait["source"],
            "discovery_external_study_accession": trait["study_accession"],
            "discovery_external_PMID": trait["PMID"], "discovery_external_DOI": trait["DOI"],
            "results_accessed_before_lock": "NO",
            "source_curation_status": "PENDING_INDEPENDENT_SOURCE_CURATION",
            "replication_family_size": str(family_size), "replication_alpha": f"{alpha:.12g}",
        })
        output.append(record)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=SOURCE_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    provenance = {
        "schema_version": "1.0.0", "source_priority": str(args.priority),
        "source_priority_sha256": sha256(args.priority), "candidate_pair_count": family_size,
        "replication_alpha": alpha,
        "selection": "Tier A/B pairs without replication class REPLICATED",
        "result_free_template": True,
        "warning": "This queue is not a locked replication manifest and contains no replication result.",
        "output": str(args.out), "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    candidate_lock = {
        "schema_version": "1.0.0",
        "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": "before_independent_source_curation_or_replication_result_access",
        "results_accessed_before_lock": False,
        "pair_count": family_size, "pair_ids_in_locked_order": [row["pair_id"] for row in output],
        "bonferroni_alpha": alpha, "source_priority_sha256": sha256(args.priority),
        "queue_template_sha256": sha256(args.out),
        "family_definition": "all Tier A/B pairs not already classified REPLICATED",
    }
    args.candidate_lock_out.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_lock_out.write_text(json.dumps(candidate_lock, indent=2, sort_keys=True) + "\n")
    print(f"REPLICATION_SOURCE_QUEUE_OK pairs={family_size} alpha={alpha:.12g} pending={family_size}")


if __name__ == "__main__":
    main()
