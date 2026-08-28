#!/usr/bin/env python3
"""Freeze the strongest replicated locus family before fine-mapping/QTL curation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


FIELDS = [
    "queue_row_id", "pair_id", "sleep_trait", "extension_trait_id", "external_phenotype_name",
    "locus_id", "CHR", "start", "end", "lead_snp", "lead_P_PLACO_PLUS",
    "discovery_rg", "discovery_fdr", "replication_class", "priority_tier",
    "comparison_type", "comparison_role", "qtl_modality",
    "dataset1_id", "dataset1_role", "dataset1_type", "dataset1_N", "dataset1_effective_N",
    "dataset1_case_fraction", "dataset1_sdY", "dataset1_prior_method",
    "dataset1_prior_source_path", "dataset1_prior_source_sha256",
    "dataset2_id", "dataset2_role", "dataset2_type", "dataset2_N", "dataset2_effective_N",
    "dataset2_case_fraction", "dataset2_sdY", "dataset2_prior_method",
    "dataset2_prior_source_path", "dataset2_prior_source_sha256",
    "molecular_feature_id", "molecular_feature_name", "tissue_cell_context",
    "source_search_status", "source_search_log_path", "source_search_log_sha256",
    "source_id", "source_url", "source_access_date", "source_license",
    "phenotype_compatibility", "sample_overlap_status",
    "summary1_path", "summary1_sha256", "summary2_path", "summary2_sha256",
    "ld_reference_id", "ld_path", "ld_sha256", "ld_variant_order_path", "ld_variant_order_sha256",
    "ld_source_type", "ld_ancestry", "ld_build", "ld_sample_size",
    "build", "ancestry", "schema_status", "effect_allele_alignment_status", "dense_variant_status",
    "single_signal_fallback_justification", "results_accessed_before_lock", "curator", "curation_date", "notes",
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


def template(locus: dict[str, str], priority: dict[str, str], modality: str) -> dict[str, str]:
    comparison = "TRAIT_TRAIT" if modality == "NONE" else f"{modality}_SEARCH"
    row_id = f"{locus['pair_id']}__{locus['locus_id']}__{comparison}"
    record = {field: "PENDING" for field in FIELDS}
    record.update({
        "queue_row_id": row_id,
        "pair_id": locus["pair_id"],
        "sleep_trait": locus["sleep_trait"],
        "extension_trait_id": locus["extension_trait_id"],
        "external_phenotype_name": locus["external_phenotype_name"],
        "locus_id": locus["locus_id"],
        "CHR": locus["CHR"],
        "start": locus["start"],
        "end": locus["end"],
        "lead_snp": locus["lead_snp"],
        "lead_P_PLACO_PLUS": locus["lead_P_PLACO_PLUS"],
        "discovery_rg": priority["rg"],
        "discovery_fdr": priority["extension_fdr"],
        "replication_class": "REPLICATED",
        "priority_tier": "B",
        "comparison_type": comparison,
        "comparison_role": "PRIMARY" if modality == "NONE" else "SOURCE_SEARCH",
        "qtl_modality": modality,
        "dataset1_id": locus["sleep_trait"],
        "dataset1_role": "SLEEP_TRAIT",
        "dataset1_prior_method": "FLAT",
        "dataset2_id": locus["extension_trait_id"] if modality == "NONE" else "PENDING_MOLECULAR_QTL",
        "dataset2_role": "EXTERNAL_TRAIT" if modality == "NONE" else "MOLECULAR_QTL",
        "dataset2_prior_method": "FLAT",
        "source_search_status": "PENDING",
        "build": "GRCh37",
        "ancestry": "EUR",
        "ld_ancestry": "EUR",
        "ld_build": "GRCh37",
        "single_signal_fallback_justification": "NOT_JUSTIFIED",
        "results_accessed_before_lock": "NO",
    })
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci", type=Path, default=Path("discovery_extension/results/pleiotropy/novel_shared_loci.tsv"))
    parser.add_argument("--priority", type=Path, default=Path("discovery_extension/results/prioritization/novel_hit_priority.tsv"))
    parser.add_argument("--replication", type=Path, default=Path("discovery_extension/results/replication/replication_results.tsv"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/fine_mapping_colocalization_contract.json"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_input_queue.tsv"))
    parser.add_argument("--candidate-lock-out", type=Path, default=Path("discovery_extension/config/fine_mapping_candidate_family.lock.json"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/fine_mapping_input_queue.json"))
    args = parser.parse_args()

    loci = read_tsv(args.loci)
    priority = {row["pair_id"]: row for row in read_tsv(args.priority)}
    replication = {row["pair_id"]: row for row in read_tsv(args.replication)}
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    eligible: list[tuple[dict[str, str], dict[str, str]]] = []
    seen_loci: set[tuple[str, str]] = set()
    for locus in loci:
        key = (locus["pair_id"], locus["locus_id"])
        if key in seen_loci:
            raise SystemExit(f"ERROR: duplicate upstream locus: {key}")
        seen_loci.add(key)
        pair_priority = priority.get(locus["pair_id"])
        pair_replication = replication.get(locus["pair_id"])
        if pair_priority is None or pair_replication is None:
            continue
        if pair_priority["priority_tier"] == "B" and pair_replication["replication_class"] == "REPLICATED":
            eligible.append((locus, pair_priority))
    if not eligible:
        raise SystemExit("ERROR: no independently replicated Tier B pleiotropic locus is available for fine-mapping")
    eligible.sort(key=lambda item: (
        float(item[0]["lead_P_PLACO_PLUS"]),
        float(item[1]["extension_fdr"]),
        -abs(float(item[1]["rg"])),
        item[0]["pair_id"],
        int(item[0]["CHR"]),
        int(item[0]["start"]),
        item[0]["locus_id"],
    ))
    maximum = int(contract["selection"]["maximum_loci"])
    selected = eligible[:maximum]
    output: list[dict[str, str]] = []
    for locus, pair_priority in selected:
        for modality in ("NONE", "EQTL", "SQTL", "PQTL"):
            output.append(template(locus, pair_priority, modality))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)

    locus_keys = [f"{locus['pair_id']}__{locus['locus_id']}" for locus, _ in selected]
    family_lock = {
        "schema_version": "1.0.0",
        "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": "after_upstream_locus_evidence_but_before_fine_mapping_QTL_curation_or_result_access",
        "results_accessed_before_lock": False,
        "eligible_locus_count": len(eligible),
        "selected_locus_count": len(selected),
        "maximum_loci": maximum,
        "locus_keys_in_locked_order": locus_keys,
        "required_template_rows_per_locus": ["TRAIT_TRAIT", "EQTL_SEARCH", "SQTL_SEARCH", "PQTL_SEARCH"],
        "loci_sha256": sha256(args.loci),
        "priority_sha256": sha256(args.priority),
        "replication_sha256": sha256(args.replication),
        "contract_sha256": sha256(args.contract),
        "queue_template_sha256": sha256(args.out),
        "selection_rule": contract["selection"]["ranking"],
    }
    args.candidate_lock_out.parent.mkdir(parents=True, exist_ok=True)
    args.candidate_lock_out.write_text(json.dumps(family_lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance = {
        "schema_version": "1.0.0",
        "result_free_template": True,
        "eligible_locus_count": len(eligible),
        "selected_locus_count": len(selected),
        "queue_row_count": len(output),
        "locus_keys_in_locked_order": locus_keys,
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "warning": "The queue contains no fine-mapping or colocalization result. All trait and molecular-QTL inputs, LD, context, and unavailable-source searches must be curated before lock.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"FINEMAPPING_QUEUE_OK loci={len(selected)} rows={len(output)} result_free=true")


if __name__ == "__main__":
    main()
