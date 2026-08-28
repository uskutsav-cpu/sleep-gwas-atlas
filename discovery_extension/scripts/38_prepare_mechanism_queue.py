#!/usr/bin/env python3
"""Freeze the prior-robust shared-signal family and its mechanistic search tasks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


TASK_SOURCES = {
    "MOLECULAR_QTL": ["GTEX_PORTAL", "EQTL_CATALOGUE", "PSYCHENCODE_SYNAPSE", "GWAS_CATALOG", "OPEN_TARGETS_PLATFORM"],
    "VARIANT_REGULATORY_ELEMENT": ["ENCODE_DCC", "SCREEN_CCRE", "OPEN_TARGETS_PLATFORM"],
    "REGULATORY_ELEMENT_TARGET_GENE": ["SCREEN_CCRE", "ENCODE_DCC", "OPEN_TARGETS_PLATFORM", "PSYCHENCODE_SYNAPSE"],
    "TISSUE_CELL_EXPRESSION": ["GTEX_PORTAL", "SINGLE_CELL_EXPRESSION_ATLAS", "PSYCHENCODE_SYNAPSE"],
    "CELL_STATE_ACCESSIBILITY": ["ENCODE_DCC", "SCREEN_CCRE", "PSYCHENCODE_SYNAPSE"],
    "PATHWAY": ["REACTOME", "OPEN_TARGETS_PLATFORM"],
    "MODEL_SYSTEM_PERTURBATION": ["BIOSTUDIES_GEO", "ENCODE_DCC", "OPEN_TARGETS_PLATFORM"],
}
FIELDS = [
    "search_task_id", "comparison_id", "pair_id", "locus_id", "comparison_type",
    "dataset1_id", "dataset2_id", "molecular_feature_id", "signal1", "signal2",
    "top_shared_variant", "top_shared_variant_PP_H4", "PP_H4", "PP_H4_over_PP_H3",
    "prior_robust", "search_task", "allowed_source_ids", "search_query_scope",
    "source_search_status", "exact_source_id", "exact_dataset_release", "exact_accession",
    "source_url", "source_access_date", "source_license_or_terms", "snapshot_path", "snapshot_sha256",
    "results_accessed_before_plan_lock", "curator", "curation_date", "notes",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fine-mapping", type=Path, default=Path("discovery_extension/results/fine_mapping/fine_mapping_colocalization.tsv"))
    parser.add_argument("--fine-mapping-provenance", type=Path, default=Path("discovery_extension/provenance/fine_mapping_colocalization_results.json"))
    parser.add_argument("--contract", type=Path, default=Path("discovery_extension/config/mechanistic_annotation_contract.json"))
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/mechanistic_sources.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/mechanism/mechanistic_search_queue.tsv"))
    parser.add_argument("--family-lock-out", type=Path, default=Path("discovery_extension/config/mechanistic_signal_family.lock.json"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/mechanistic_search_queue.json"))
    args = parser.parse_args()

    results = read_tsv(args.fine_mapping)
    fine_mapping_provenance = json.loads(args.fine_mapping_provenance.read_text(encoding="utf-8"))
    if fine_mapping_provenance.get("output_sha256") != sha256(args.fine_mapping):
        raise SystemExit("ERROR: fine-mapping/colocalization artifact differs from its provenance")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    sources = {row["source_id"] for row in read_tsv(args.sources)}
    if set(TASK_SOURCES) != set(contract["required_search_tasks_per_signal"]):
        raise SystemExit("ERROR: mechanistic task/source routing differs from contract")
    if any(not set(route).issubset(sources) for route in TASK_SOURCES.values()):
        raise SystemExit("ERROR: mechanistic task route points outside source registry")
    eligible = [
        row for row in results
        if row["prior_role"] == "PRIMARY"
        and row["coloc_method"] == "COLOC_SUSIE"
        and row["primary_shared_signal_rule_pass"] == "True"
        and row["prior_robust"] == "True"
        and row["colocalization_interpretation"] == "SHARED_SIGNAL_MODEL_SUPPORTED"
        and row["analysis_status"] == "COLOC_SUSIE_COMPLETE"
    ]
    if not eligible:
        raise SystemExit("ERROR: no prior-robust primary coloc-SuSiE shared signal is available for mechanistic annotation")
    eligible.sort(key=lambda row: (
        row["pair_id"], row["locus_id"], row["comparison_type"], row["signal1"], row["signal2"], row["comparison_id"],
    ))
    signal_keys: list[str] = []
    output: list[dict[str, str]] = []
    for row in eligible:
        signal_key = f"{row['comparison_id']}__{row['signal1']}__{row['signal2']}"
        if signal_key in signal_keys:
            raise SystemExit(f"ERROR: duplicate eligible mechanistic signal: {signal_key}")
        signal_keys.append(signal_key)
        for task in contract["required_search_tasks_per_signal"]:
            task_id = f"{signal_key}__{task}"
            record = {field: "PENDING" for field in FIELDS}
            record.update({
                "search_task_id": task_id, "comparison_id": row["comparison_id"],
                "pair_id": row["pair_id"], "locus_id": row["locus_id"],
                "comparison_type": row["comparison_type"], "dataset1_id": row["dataset1_id"],
                "dataset2_id": row["dataset2_id"], "molecular_feature_id": row["molecular_feature_id"],
                "signal1": row["signal1"], "signal2": row["signal2"],
                "top_shared_variant": row["top_shared_variant"],
                "top_shared_variant_PP_H4": row["top_shared_variant_PP_H4"],
                "PP_H4": row["PP_H4"], "PP_H4_over_PP_H3": row["PP_H4_over_PP_H3"],
                "prior_robust": "True", "search_task": task,
                "allowed_source_ids": ";".join(TASK_SOURCES[task]),
                "search_query_scope": f"PENDING_EXACT_{task}_QUERY_FOR_{row['locus_id']}_{row['signal1']}_{row['signal2']}",
                "source_search_status": "PENDING", "results_accessed_before_plan_lock": "NO",
            })
            output.append(record)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    lock = {
        "schema_version": "1.0.0",
        "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": contract["selection_timing"],
        "results_accessed_before_plan_lock": False,
        "signal_count": len(signal_keys),
        "signal_keys_in_locked_order": signal_keys,
        "search_task_count": len(output),
        "search_task_ids_in_locked_order": [row["search_task_id"] for row in output],
        "required_search_tasks_per_signal": contract["required_search_tasks_per_signal"],
        "fine_mapping_sha256": sha256(args.fine_mapping),
        "fine_mapping_provenance_sha256": sha256(args.fine_mapping_provenance),
        "contract_sha256": sha256(args.contract),
        "sources_sha256": sha256(args.sources),
        "queue_template_sha256": sha256(args.out),
        "claim_limit": contract["claim_limit"],
    }
    args.family_lock_out.parent.mkdir(parents=True, exist_ok=True)
    args.family_lock_out.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance = {
        "schema_version": "1.0.0",
        "result_free_search_plan": True,
        "signal_count": len(signal_keys),
        "search_task_count": len(output),
        "output": str(args.out),
        "output_sha256": sha256(args.out),
        "warning": "The queue is a result-free search plan. Landing pages and PENDING rows are not mechanistic evidence.",
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"MECHANISM_QUEUE_OK signals={len(signal_keys)} tasks={len(output)} result_free=true")


if __name__ == "__main__":
    main()
