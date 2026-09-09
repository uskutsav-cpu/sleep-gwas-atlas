#!/usr/bin/env python3
"""Classify the exact locked independent-replication family."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


OUTPUT_FIELDS = [
    "pair_id", "sleep_trait", "extension_trait_id", "external_phenotype_name",
    "replication_source_id", "replication_study_accession", "replication_phenotype_definition",
    "participant_overlap_status", "discovery_rg", "discovery_se", "discovery_fdr",
    "replication_rg", "replication_se", "replication_z", "replication_p", "replication_alpha",
    "direction_concordant", "heterogeneity_z", "heterogeneity_Q", "heterogeneity_p",
    "replication_h2", "replication_h2_se", "replication_h2_z", "replication_LDSC_intercept",
    "replication_h2_pass", "cross_trait_LDSC_intercept", "cross_trait_LDSC_intercept_se",
    "snp_overlap_valid_alleles", "ancestry", "analysis_status",
    "replication_search_databases", "replication_search_queries", "replication_search_date",
    "replication_search_evidence", "unavailable_reason", "replication_class",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest", type=Path,
        default=Path("discovery_extension/config/replication_manifest.tsv"),
    )
    parser.add_argument(
        "--lock", type=Path,
        default=Path("discovery_extension/config/replication_manifest.lock.json"),
    )
    parser.add_argument("--h2", type=Path, required=True)
    parser.add_argument("--rg", type=Path, required=True)
    parser.add_argument(
        "--out", type=Path,
        default=Path("discovery_extension/results/replication/replication_results.tsv"),
    )
    parser.add_argument(
        "--provenance-out", type=Path,
        default=Path("discovery_extension/provenance/replication_results.json"),
    )
    args = parser.parse_args()
    manifest = read_tsv(args.manifest)
    lock = json.loads(args.lock.read_text())
    if sha256(args.manifest) != lock["manifest_sha256"]:
        raise SystemExit("ERROR: replication manifest differs from its pre-result lock")
    pair_ids = [row["pair_id"] for row in manifest]
    if pair_ids != lock["pair_ids_in_locked_order"] or len(pair_ids) != lock["pair_count"]:
        raise SystemExit("ERROR: replication manifest membership/order differs from lock")
    manifest_by_pair = {row["pair_id"]: row for row in manifest}
    h2_rows = read_tsv(args.h2)
    h2_by_source = {row["replication_source_id"]: row for row in h2_rows}
    if len(h2_by_source) != len(h2_rows):
        raise SystemExit("ERROR: duplicate replication_source_id in replication h2 table")
    testable_pair_ids = list(lock.get("testable_pair_ids_in_locked_order", pair_ids))
    unavailable_pair_ids = list(lock.get("unavailable_pair_ids_in_locked_order", []))
    if set(testable_pair_ids) | set(unavailable_pair_ids) != set(pair_ids) or set(testable_pair_ids) & set(unavailable_pair_ids):
        raise SystemExit("ERROR: lock does not partition the replication candidate family")
    rg_rows = read_tsv(args.rg)
    expected_sources = {manifest_by_pair[pair_id]["replication_source_id"] for pair_id in testable_pair_ids}
    if set(h2_by_source) != expected_sources:
        raise SystemExit("ERROR: replication h2 results are not the exact locked source family")
    h2_pass_sources = {
        source_id for source_id, row in h2_by_source.items()
        if row["primary_status"] == "PASS" and float(row["h2_z"]) >= 4 and float(row["LDSC_intercept"]) <= 1.2
    }
    expected_rg_pair_ids = {
        pair_id for pair_id in testable_pair_ids
        if manifest_by_pair[pair_id]["replication_source_id"] in h2_pass_sources
    }
    if {row["pair_id"] for row in rg_rows} != expected_rg_pair_ids or len(rg_rows) != len(expected_rg_pair_ids):
        raise SystemExit("ERROR: replication rg results are not the exact h2-pass subset of the locked source-available family")
    rg_by_pair = {row["pair_id"]: row for row in rg_rows}

    alpha = float(lock["bonferroni_alpha"])
    output: list[dict[str, object]] = []
    for pair_id in pair_ids:
        source = manifest_by_pair[pair_id]
        if pair_id in unavailable_pair_ids:
            output.append({
                "pair_id": pair_id, "sleep_trait": source["sleep_trait"],
                "extension_trait_id": source["extension_trait_id"],
                "external_phenotype_name": source["external_phenotype_name"],
                "replication_source_id": "NA", "replication_study_accession": "NA",
                "replication_phenotype_definition": "NA", "participant_overlap_status": "NA",
                "discovery_rg": source["discovery_rg"], "discovery_se": source["discovery_se"],
                "discovery_fdr": source["discovery_fdr"], "replication_rg": "NA",
                "replication_se": "NA", "replication_z": "NA", "replication_p": "NA",
                "replication_alpha": alpha, "direction_concordant": "NA",
                "heterogeneity_z": "NA", "heterogeneity_Q": "NA", "heterogeneity_p": "NA",
                "replication_h2": "NA", "replication_h2_se": "NA", "replication_h2_z": "NA",
                "replication_LDSC_intercept": "NA", "replication_h2_pass": "NA",
                "cross_trait_LDSC_intercept": "NA", "cross_trait_LDSC_intercept_se": "NA",
                "snp_overlap_valid_alleles": "NA", "ancestry": "NA",
                "analysis_status": "NO_INDEPENDENT_DATASET_PRE_RESULT_CLASSIFICATION",
                "replication_search_databases": source["replication_search_databases"],
                "replication_search_queries": source["replication_search_queries"],
                "replication_search_date": source["replication_search_date"],
                "replication_search_evidence": source["replication_search_evidence"],
                "unavailable_reason": source["unavailable_reason"],
                "replication_class": "NO_INDEPENDENT_DATASET",
            })
            continue
        if source["replication_source_id"] not in h2_by_source:
            raise SystemExit(f"ERROR: replication h2 result missing: {pair_id}")
        h2 = h2_by_source[source["replication_source_id"]]
        h2_z, h2_intercept = float(h2["h2_z"]), float(h2["LDSC_intercept"])
        h2_pass = h2["primary_status"] == "PASS" and h2_z >= 4 and h2_intercept <= 1.2
        if not h2_pass:
            output.append({
                "pair_id": pair_id, "sleep_trait": source["sleep_trait"],
                "extension_trait_id": source["extension_trait_id"],
                "external_phenotype_name": source["external_phenotype_name"],
                "replication_source_id": source["replication_source_id"],
                "replication_study_accession": source["replication_study_accession"],
                "replication_phenotype_definition": source["replication_phenotype_definition"],
                "participant_overlap_status": source["participant_overlap_status"],
                "discovery_rg": source["discovery_rg"], "discovery_se": source["discovery_se"],
                "discovery_fdr": source["discovery_fdr"], "replication_rg": "NA",
                "replication_se": "NA", "replication_z": "NA", "replication_p": "NA",
                "replication_alpha": alpha, "direction_concordant": "NA",
                "heterogeneity_z": "NA", "heterogeneity_Q": "NA", "heterogeneity_p": "NA",
                "replication_h2": h2["h2"], "replication_h2_se": h2["h2_se"],
                "replication_h2_z": h2_z, "replication_LDSC_intercept": h2_intercept,
                "replication_h2_pass": "False", "cross_trait_LDSC_intercept": "NA",
                "cross_trait_LDSC_intercept_se": "NA", "snp_overlap_valid_alleles": "NA",
                "ancestry": source["ancestry"], "analysis_status": "REPLICATION_H2_FAILED_NO_PAIR_TEST",
                "replication_search_databases": source["replication_search_databases"],
                "replication_search_queries": source["replication_search_queries"],
                "replication_search_date": source["replication_search_date"],
                "replication_search_evidence": source["replication_search_evidence"],
                "unavailable_reason": "Replication source acquired but failed the prespecified h2 Z/intercept gate; pairwise rg was not run.",
                "replication_class": "UNDERPOWERED",
            })
            continue
        result = rg_by_pair[pair_id]
        if result["sleep_trait"] != source["sleep_trait"] or result["replication_source_id"] != source["replication_source_id"]:
            raise SystemExit(f"ERROR: replication result identity differs from lock: {pair_id}")
        if result["ancestry"] != source["ancestry"] or result["analysis_status"] != "REPLICATION_RG_COMPLETE":
            raise SystemExit(f"ERROR: replication ancestry/status failed: {pair_id}")
        discovery_rg, discovery_se = float(source["discovery_rg"]), float(source["discovery_se"])
        replication_rg, replication_se = float(result["rg"]), float(result["se"])
        p_value = float(result["p"])
        if not all(math.isfinite(value) for value in (discovery_rg, discovery_se, replication_rg, replication_se, p_value)) or discovery_se <= 0 or replication_se <= 0 or not 0 <= p_value <= 1:
            raise SystemExit(f"ERROR: invalid replication estimate or p-value: {pair_id}")
        concordant = discovery_rg * replication_rg > 0
        heterogeneity_z = (replication_rg - discovery_rg) / math.sqrt(discovery_se**2 + replication_se**2)
        heterogeneity_q = heterogeneity_z**2
        heterogeneity_p = math.erfc(math.sqrt(heterogeneity_q / 2))
        if concordant and p_value < alpha:
            replication_class = "REPLICATED"
        elif concordant:
            replication_class = "DIRECTIONALLY_CONCORDANT"
        else:
            replication_class = "FAILED_REPLICATION"
        output.append({
            "pair_id": pair_id, "sleep_trait": source["sleep_trait"],
            "extension_trait_id": source["extension_trait_id"],
            "external_phenotype_name": source["external_phenotype_name"],
            "replication_source_id": source["replication_source_id"],
            "replication_study_accession": source["replication_study_accession"],
            "replication_phenotype_definition": source["replication_phenotype_definition"],
            "participant_overlap_status": source["participant_overlap_status"],
            "discovery_rg": discovery_rg, "discovery_se": discovery_se,
            "discovery_fdr": source["discovery_fdr"],
            "replication_rg": replication_rg, "replication_se": replication_se,
            "replication_z": result["z"], "replication_p": p_value,
            "replication_alpha": alpha, "direction_concordant": str(concordant),
            "heterogeneity_z": heterogeneity_z, "heterogeneity_Q": heterogeneity_q,
            "heterogeneity_p": heterogeneity_p,
            "replication_h2": h2["h2"], "replication_h2_se": h2["h2_se"],
            "replication_h2_z": h2_z, "replication_LDSC_intercept": h2_intercept,
            "replication_h2_pass": str(h2_pass),
            "cross_trait_LDSC_intercept": result["cross_trait_LDSC_intercept"],
            "cross_trait_LDSC_intercept_se": result["cross_trait_LDSC_intercept_se"],
            "snp_overlap_valid_alleles": result["snp_overlap_valid_alleles"],
            "ancestry": result["ancestry"], "analysis_status": result["analysis_status"],
            "replication_search_databases": source["replication_search_databases"],
            "replication_search_queries": source["replication_search_queries"],
            "replication_search_date": source["replication_search_date"],
            "replication_search_evidence": source["replication_search_evidence"],
            "unavailable_reason": "NA",
            "replication_class": replication_class,
        })
    output.sort(key=lambda row: pair_ids.index(str(row["pair_id"])))
    write_tsv(args.out, OUTPUT_FIELDS, output)
    counts: dict[str, int] = {}
    for row in output:
        key = str(row["replication_class"])
        counts[key] = counts.get(key, 0) + 1
    provenance = {
        "schema_version": "1.0.0", "manifest_sha256": sha256(args.manifest),
        "lock_sha256": sha256(args.lock), "h2_sha256": sha256(args.h2),
        "rg_sha256": sha256(args.rg), "pair_count": len(output),
        "bonferroni_alpha": alpha, "classification_counts": dict(sorted(counts.items())),
        "classification_rule": {
            "REPLICATED": "h2 pass, direction concordant, p below locked Bonferroni alpha",
            "DIRECTIONALLY_CONCORDANT": "h2 pass and direction concordant but p not below locked alpha",
            "UNDERPOWERED": "replication h2 gate failed",
            "FAILED_REPLICATION": "h2 pass but effect direction discordant",
            "NO_INDEPENDENT_DATASET": "completed pre-result search found no qualifying independent source",
        },
        "output": str(args.out), "output_sha256": sha256(args.out),
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    print(f"REPLICATION_RESULTS_OK pairs={len(output)} classes={dict(sorted(counts.items()))}")


if __name__ == "__main__":
    main()
