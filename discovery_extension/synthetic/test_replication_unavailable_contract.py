#!/usr/bin/env python3
"""Exercise locked replication plus an evidence-backed unavailable outcome."""

from __future__ import annotations

import csv
import hashlib
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    panel = read_tsv(ROOT / "discovery_extension/config/candidate_traits.tsv")
    sleep = next(row for row in read_tsv(ROOT / "config/analysis_panel.tsv") if row["domain"] == "sleep")
    with tempfile.TemporaryDirectory(prefix="synthetic_replication_contract_") as temporary:
        work = Path(temporary)
        priority = work / "priority.tsv"
        priority_rows = []
        for index, trait in enumerate(panel[:2]):
            priority_rows.append({
                "pair_id": f"{sleep['trait_id']}__{trait['extension_trait_id']}",
                "sleep_trait": sleep["trait_id"], "extension_trait_id": trait["extension_trait_id"],
                "phenotype_name": trait["phenotype_name"], "rg": "0.2", "se": "0.03",
                "extension_fdr": "0.01", "pair_novelty_class": "APPARENTLY_NOVEL",
                "pair_novelty_strength": "MODERATE", "priority_tier": "A",
                "replication_class": "NOT_YET_ATTEMPTED",
            })
        write_tsv(priority, list(priority_rows[0]), priority_rows)
        queue, candidate_lock = work / "queue.tsv", work / "candidate.lock.json"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/21_prepare_replication_queue.py"),
            "--priority", str(priority), "--panel", str(ROOT / "discovery_extension/config/candidate_traits.tsv"),
            "--core", str(ROOT / "config/analysis_panel.tsv"), "--out", str(queue),
            "--provenance-out", str(work / "queue.json"), "--candidate-lock-out", str(candidate_lock),
        ], check=True)
        rows = read_tsv(queue)
        source = work / "independent.tsv.gz"
        source.write_bytes(b"SYNTHETIC INDEPENDENT GWAS - NOT REAL RESULTS\n")
        source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
        rows[0].update({
            "replication_source_id": "synthetic_independent", "replication_study_accession": "SYNTHETIC001",
            "replication_publication": "synthetic test only", "replication_PMID": "00000000",
            "replication_DOI": "10.0000/synthetic", "replication_source_url": "https://example.org/synthetic.gz",
            "replication_checksum": f"sha256:{source_sha}", "replication_local_path": str(source),
            "replication_phenotype_definition": "synthetic exact phenotype", "phenotype_match_status": "EXACT",
            "ancestry": "EUR", "build": "GRCh37", "sample_size": "50000", "cases": "NA", "controls": "NA",
            "discovery_cohort_relation": "NON_UKB", "participant_overlap_status": "NON_OVERLAPPING_CONFIRMED",
            "participant_overlap_evidence": "synthetic disjoint cohorts", "source_identity_status": "VERIFIED",
            "schema_status": "VERIFIED", "effect_allele_status": "UNAMBIGUOUS",
            "full_resolution_availability": "YES", "source_curation_status": "COMPLETE_BEFORE_RESULTS",
            "replication_search_databases": "GWAS Catalog;PubMed", "replication_search_queries": "synthetic query one",
            "replication_search_date": "2099-01-01", "replication_search_evidence": "synthetic source record",
            "unavailable_reason": "NA",
        })
        rows[1].update({
            "replication_search_databases": "GWAS Catalog;PubMed;consortium repositories",
            "replication_search_queries": "synthetic exact phenotype || synthetic synonym GWAS",
            "replication_search_date": "2099-01-01",
            "replication_search_evidence": "synthetic exhaustive negative search record",
            "unavailable_reason": "No qualifying non-overlapping EUR full-summary-statistics source found in synthetic test",
            "source_curation_status": "NO_INDEPENDENT_DATASET_COMPLETE_BEFORE_RESULTS",
        })
        write_tsv(queue, list(rows[0]), rows)
        manifest, manifest_lock = work / "manifest.tsv", work / "manifest.lock.json"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/22_lock_replication_manifest.py"),
            "--queue", str(queue), "--candidate-lock", str(candidate_lock),
            "--out", str(manifest), "--lock", str(manifest_lock),
        ], check=True)
        h2 = work / "h2.tsv"
        write_tsv(h2, ["replication_source_id", "h2", "h2_se", "h2_z", "LDSC_intercept", "primary_status"], [{
            "replication_source_id": "synthetic_independent", "h2": "0.2", "h2_se": "0.02",
            "h2_z": "10", "LDSC_intercept": "1.0", "primary_status": "PASS",
        }])
        rg = work / "rg.tsv"
        write_tsv(rg, [
            "pair_id", "sleep_trait", "replication_source_id", "rg", "se", "z", "p",
            "cross_trait_LDSC_intercept", "cross_trait_LDSC_intercept_se",
            "snp_overlap_valid_alleles", "ancestry", "analysis_status",
        ], [{
            "pair_id": rows[0]["pair_id"], "sleep_trait": sleep["trait_id"],
            "replication_source_id": "synthetic_independent", "rg": "0.18", "se": "0.04",
            "z": "4.5", "p": "0.01", "cross_trait_LDSC_intercept": "0.0",
            "cross_trait_LDSC_intercept_se": "0.005", "snp_overlap_valid_alleles": "900000",
            "ancestry": "EUR", "analysis_status": "REPLICATION_RG_COMPLETE",
        }])
        results = work / "results.tsv"
        subprocess.run([
            "python3", str(ROOT / "discovery_extension/scripts/23_collate_replication.py"),
            "--manifest", str(manifest), "--lock", str(manifest_lock), "--h2", str(h2), "--rg", str(rg),
            "--out", str(results), "--provenance-out", str(work / "results.json"),
        ], check=True)
        result_rows = read_tsv(results)
        if [row["replication_class"] for row in result_rows] != ["REPLICATED", "NO_INDEPENDENT_DATASET"]:
            raise SystemExit(f"ERROR: synthetic replication classifications are incorrect: {result_rows}")
        if result_rows[1]["analysis_status"] != "NO_INDEPENDENT_DATASET_PRE_RESULT_CLASSIFICATION":
            raise SystemExit("ERROR: unavailable replication outcome lost its pre-result status")
        if float(result_rows[0]["heterogeneity_p"]) <= 0 or float(result_rows[0]["heterogeneity_p"]) > 1:
            raise SystemExit("ERROR: replication heterogeneity was not reported")
    print("REPLICATION_UNAVAILABLE_SYNTHETIC_OK pairs=2 replicated=1 unavailable=1 family_locked=true isolated=true")


if __name__ == "__main__":
    main()
