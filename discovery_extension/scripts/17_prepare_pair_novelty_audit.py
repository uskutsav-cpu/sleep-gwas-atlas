#!/usr/bin/env python3
"""Create, but never pre-answer, the audit for extension-FDR-significant pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


FIELDS = [
    "pair_id", "sleep_trait", "extension_trait_id", "phenotype_name", "phenotype_domain",
    "rg", "se", "p", "extension_fdr", "abs_rg_ge_0.15",
    "preanalysis_novelty_priority", "preanalysis_prior_screen_coverage",
    "audit_status", "direct_prior_same_pair", "same_sleep_trait_context",
    "same_or_equivalent_phenotype", "same_direction", "broad_phenome_screen_overlap",
    "near_neighbor_evidence", "discovery_vs_replication_in_prior_work",
    "exact_prior_rg_found", "closest_prior_result", "prior_method", "prior_effect",
    "prior_publication", "prior_DOI", "prior_PMID", "search_databases",
    "search_queries_used", "search_date", "evidence_PMIDs_DOIs_URLs",
    "biological_plausibility", "connection_obviousness",
    "independent_replication_dataset_availability", "dense_summary_statistics_available",
    "molecular_qtl_data_available", "independent_replication_class",
    "novelty_class", "novelty_strength",
    "decision_rationale", "reviewer_notes", "reviewer",
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
        "--rg", type=Path,
        default=Path("discovery_extension/results/ldsc/extension_rg_matrix.tsv"),
    )
    parser.add_argument(
        "--panel", type=Path,
        default=Path("discovery_extension/config/candidate_traits.tsv"),
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path("discovery_extension/results/novelty/extension_novelty_audit.tsv"),
    )
    parser.add_argument(
        "--provenance-out", type=Path,
        default=Path("discovery_extension/provenance/pair_novelty_audit_template.json"),
    )
    args = parser.parse_args()
    if not args.rg.is_file():
        raise SystemExit(f"ERROR: primary extension rg table is missing: {args.rg}")
    panel = {row["extension_trait_id"]: row for row in read_tsv(args.panel)}
    rg = read_tsv(args.rg)
    hits = [row for row in rg if float(row["extension_fdr"]) < 0.05]
    output: list[dict[str, str]] = []
    for row in hits:
        trait_id = row["extension_trait_id"]
        if trait_id not in panel:
            raise SystemExit(f"ERROR: rg result is outside locked extension panel: {trait_id}")
        trait = panel[trait_id]
        output.append({
            "pair_id": f"{row['sleep_trait']}__{trait_id}",
            "sleep_trait": row["sleep_trait"], "extension_trait_id": trait_id,
            "phenotype_name": trait["phenotype_name"], "phenotype_domain": trait["phenotype_domain"],
            "rg": row["rg"], "se": row["se"], "p": row["p"],
            "extension_fdr": row["extension_fdr"],
            "abs_rg_ge_0.15": str(abs(float(row["rg"])) >= 0.15),
            "preanalysis_novelty_priority": trait["novelty_priority"],
            "preanalysis_prior_screen_coverage": trait["prior_sleep_screen_coverage"],
            "audit_status": "PENDING_MANUAL_LITERATURE_AUDIT",
            "direct_prior_same_pair": "PENDING", "same_sleep_trait_context": "PENDING",
            "same_or_equivalent_phenotype": "PENDING", "same_direction": "PENDING",
            "broad_phenome_screen_overlap": "PENDING", "near_neighbor_evidence": "PENDING",
            "discovery_vs_replication_in_prior_work": "PENDING",
            "exact_prior_rg_found": "PENDING", "closest_prior_result": "PENDING",
            "prior_method": "PENDING", "prior_effect": "PENDING",
            "prior_publication": "PENDING", "prior_DOI": "PENDING", "prior_PMID": "PENDING",
            "search_databases": "PENDING", "search_queries_used": "PENDING", "search_date": "PENDING",
            "evidence_PMIDs_DOIs_URLs": "PENDING",
            "biological_plausibility": "PENDING", "connection_obviousness": "PENDING",
            "independent_replication_dataset_availability": "PENDING",
            "dense_summary_statistics_available": "PENDING", "molecular_qtl_data_available": "PENDING",
            "independent_replication_class": "NOT_YET_ATTEMPTED",
            "novelty_class": "PENDING", "novelty_strength": "PENDING",
            "decision_rationale": "PENDING", "reviewer_notes": "PENDING", "reviewer": "PENDING",
        })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    provenance = {
        "schema_version": "1.0.0", "source_rg": str(args.rg),
        "source_rg_sha256": sha256(args.rg), "extension_fdr_threshold": 0.05,
        "audit_pair_count": len(output), "template_only": True,
        "allowed_novelty_classes": ["KNOWN_REPLICATION", "KNOWN_BUT_NEW_DATASET", "PARTIAL_EXTENSION", "NO_DIRECT_RG_FOUND", "APPARENTLY_NOVEL", "UNCERTAIN"],
        "strong_novelty_rule": "APPARENTLY_NOVEL after several targeted searches, no direct prior rg, and replication class REPLICATED",
        "warning": "PENDING fields are not negative evidence and preanalysis categories are not pair-level novelty conclusions",
        "output": str(args.out), "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"PAIR_NOVELTY_AUDIT_TEMPLATE_OK pairs={len(output)} pending={len(output)}")


if __name__ == "__main__":
    main()
