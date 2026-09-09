#!/usr/bin/env python3
"""Prepare priority and globally-null pair rows for pre-result local curation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


FIELDS = [
    "pair_id", "sleep_trait", "extension_trait_id", "phenotype_name", "phenotype_domain",
    "global_rg", "global_rg_se", "global_rg_p", "global_rg_fdr", "global_analysis_status",
    "discovery_priority_tier", "replication_class", "pair_novelty_class", "pair_novelty_strength",
    "snp_overlap_valid_alleles", "selection_stratum", "selection_status",
    "biological_relevance", "polygenic_overlap_evidence", "secondary_selection_rationale",
    "dense_summary_statistics_available", "trait_a_dense_path", "trait_a_dense_sha256",
    "trait_b_dense_path", "trait_b_dense_sha256", "locus_definition_path",
    "locus_definition_sha256", "locus_count", "planned_methods", "hdl_feasibility_rationale",
    "lava_simulation_random_seed", "result_access_status", "curator", "curation_date",
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
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/local/local_analysis_queue.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/local_analysis_queue.json"))
    args = parser.parse_args()

    rows = read_tsv(args.priority)
    required = {
        "pair_id", "sleep_trait", "extension_trait_id", "phenotype_name", "phenotype_domain",
        "rg", "se", "p", "extension_fdr", "analysis_status", "priority_tier",
        "replication_class", "pair_novelty_class", "pair_novelty_strength",
        "snp_overlap_valid_alleles", "extension_h2_pass", "sleep_h2_pass",
        "dense_summary_statistics_available", "biological_plausibility",
    }
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(f"ERROR: priority table lacks fields: {sorted(required - set(rows[0] if rows else []))}")
    if len({row["pair_id"] for row in rows}) != len(rows):
        raise SystemExit("ERROR: duplicate priority pair_id")

    output: list[dict[str, str]] = []
    for row in rows:
        clean = row["analysis_status"] == "PRIMARY_EXTENSION_RG_COMPLETE"
        both_h2 = row["extension_h2_pass"] == "True" and row["sleep_h2_pass"] == "True"
        fdr = float(row["extension_fdr"])
        if row["priority_tier"] in {"A", "B"}:
            stratum = "PRIORITY_DISCOVERY"
            selected = row["dense_summary_statistics_available"] == "YES"
            selection_status = "SELECTED_PRIORITY_DISCOVERY" if selected else "BLOCKED_DENSE_SUMMARY_STATISTICS"
            rationale = "Tier A/B discovery after novelty and replication prioritization"
            relevance = row["biological_plausibility"]
        elif clean and both_h2 and fdr >= 0.05:
            stratum = "GLOBAL_NULL_SECONDARY"
            selection_status = "PENDING_GLOBAL_NULL_CURATION"
            rationale = "PENDING: document biological relevance or polygenic-overlap evidence without local-result access"
            relevance = "PENDING"
        else:
            continue
        output.append({
            "pair_id": row["pair_id"], "sleep_trait": row["sleep_trait"],
            "extension_trait_id": row["extension_trait_id"], "phenotype_name": row["phenotype_name"],
            "phenotype_domain": row["phenotype_domain"], "global_rg": row["rg"],
            "global_rg_se": row["se"], "global_rg_p": row["p"], "global_rg_fdr": row["extension_fdr"],
            "global_analysis_status": row["analysis_status"], "discovery_priority_tier": row["priority_tier"],
            "replication_class": row["replication_class"], "pair_novelty_class": row["pair_novelty_class"],
            "pair_novelty_strength": row["pair_novelty_strength"],
            "snp_overlap_valid_alleles": row["snp_overlap_valid_alleles"],
            "selection_stratum": stratum, "selection_status": selection_status,
            "biological_relevance": relevance, "polygenic_overlap_evidence": "PENDING",
            "secondary_selection_rationale": rationale,
            "dense_summary_statistics_available": row["dense_summary_statistics_available"],
            "trait_a_dense_path": "PENDING", "trait_a_dense_sha256": "PENDING",
            "trait_b_dense_path": "PENDING", "trait_b_dense_sha256": "PENDING",
            "locus_definition_path": "PENDING", "locus_definition_sha256": "PENDING",
            "locus_count": "PENDING", "planned_methods": "LAVA_PRIMARY;HDL_L_ROBUSTNESS",
            "hdl_feasibility_rationale": "PENDING", "result_access_status": "NOT_ACCESSED",
            "lava_simulation_random_seed": str(int(hashlib.sha256(row["pair_id"].encode()).hexdigest()[:8], 16) % (2**31 - 1) + 1),
            "curator": "PENDING", "curation_date": "PENDING",
        })
    output.sort(key=lambda row: (
        0 if row["selection_stratum"] == "PRIORITY_DISCOVERY" else 1,
        float(row["global_rg_fdr"]), -abs(float(row["global_rg"])), row["pair_id"],
    ))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    counts = Counter(row["selection_stratum"] for row in output)
    provenance = {
        "schema_version": "1.0.0", "source_priority": str(args.priority),
        "source_priority_sha256": sha256(args.priority), "row_count": len(output),
        "stratum_counts": dict(sorted(counts.items())), "contains_local_results": False,
        "global_null_rule": "extension_fdr>=0.05, clean primary global rg, both traits h2-pass",
        "warning": "PENDING_GLOBAL_NULL_CURATION is not selection; no local results may be accessed before lock.",
        "output": str(args.out), "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"LOCAL_ANALYSIS_QUEUE_OK rows={len(output)} strata={dict(sorted(counts.items()))} result_free=true")


if __name__ == "__main__":
    main()
