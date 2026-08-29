#!/usr/bin/env python3
"""Apply result-free source decisions to the frozen 217-pair replication family."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from importlib.util import module_from_spec, spec_from_file_location


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


def source_fields() -> list[str]:
    path = Path(__file__).with_name("21_prepare_replication_queue.py")
    spec = spec_from_file_location("prepare_replication_queue", path)
    if spec is None or spec.loader is None:
        raise SystemExit("ERROR: cannot load replication queue schema")
    module = module_from_spec(spec); spec.loader.exec_module(module)
    return list(module.SOURCE_FIELDS)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", type=Path, default=Path("discovery_extension/results/replication/replication_source_queue.tsv"))
    parser.add_argument("--candidate-lock", type=Path, default=Path("discovery_extension/config/replication_candidate_family.lock.json"))
    parser.add_argument("--sources", type=Path, default=Path("discovery_extension/config/replication_finngen_r13_sources.tsv"))
    parser.add_argument("--mvp-sources", type=Path, default=Path("discovery_extension/config/replication_mvp_sources.tsv"))
    parser.add_argument("--search", type=Path, default=Path("discovery_extension/results/replication/replication_unique_phenotype_search.tsv"))
    parser.add_argument("--out", type=Path, default=Path("discovery_extension/results/replication/replication_source_queue.tsv"))
    parser.add_argument("--provenance-out", type=Path, default=Path("discovery_extension/provenance/replication_sources/curated_queue.json"))
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()

    queue_fields, queue = read_tsv(args.queue)
    _, source_rows = read_tsv(args.sources)
    _, mvp_source_rows = read_tsv(args.mvp_sources)
    _, search_rows = read_tsv(args.search)
    lock = json.loads(args.candidate_lock.read_text(encoding="utf-8"))
    if [row["pair_id"] for row in queue] != lock.get("pair_ids_in_locked_order") or len(queue) != 217:
        raise SystemExit("ERROR: queue membership/order differs from frozen candidate family")
    if lock.get("results_accessed_before_lock") is not False:
        raise SystemExit("ERROR: candidate family was not frozen before results")
    source_by_trait = {row["extension_trait_id"]: row for row in source_rows}
    mvp_by_trait = {row["extension_trait_id"]: row for row in mvp_source_rows}
    if set(source_by_trait) & set(mvp_by_trait):
        raise SystemExit("ERROR: a trait has more than one selected independent source")
    search_by_trait = {row["extension_trait_id"]: row for row in search_rows}
    if len(search_by_trait) != 44 or any(row["search_status"] != "COMPLETE_BEFORE_RESULTS" for row in search_rows):
        raise SystemExit("ERROR: unique-phenotype source search is incomplete")

    if args.validate_only:
        expected_fields = source_fields()
        testable_rows = [row for row in queue if row.get("source_curation_status") == "COMPLETE_BEFORE_RESULTS"]
        unavailable_rows = [row for row in queue if row.get("source_curation_status") == "NO_INDEPENDENT_DATASET_COMPLETE_BEFORE_RESULTS"]
        if queue_fields != expected_fields or len(testable_rows) != 58 or len(unavailable_rows) != 159:
            raise SystemExit("ERROR: curated replication queue schema or frozen counts differ")
        if any(row.get("results_accessed_before_lock") != "NO" for row in queue):
            raise SystemExit("ERROR: replication result-access timing flag differs")
        if any(not row.get("unavailable_reason") or row["unavailable_reason"].startswith("PENDING") for row in unavailable_rows):
            raise SystemExit("ERROR: an unavailable pair lacks its pre-result reason")
        if len({row["replication_source_id"] for row in testable_rows}) != 13:
            raise SystemExit("ERROR: expected 13 independently selected replication sources")
        print(f"REPLICATION_QUEUE_VALID pairs={len(queue)} testable={len(testable_rows)} unavailable={len(unavailable_rows)} sources=13 sha256={sha256(args.queue)}")
        return

    fields = source_fields()
    output: list[dict[str, str]] = []
    testable = overlap_rejected = unavailable = 0
    for original in queue:
        row = {field: original.get(field, "PENDING") for field in fields}
        search = search_by_trait[row["extension_trait_id"]]
        row.update({
            "results_accessed_before_lock": "NO",
            "replication_search_databases": "GWAS Catalog REST API; FinnGen R13 phenotype API and public summary-statistics object store",
            "replication_search_queries": search["search_terms"],
            "replication_search_date": search["search_date"],
            "replication_search_evidence": f"{args.search}#{row['extension_trait_id']}",
        })
        source = source_by_trait.get(row["extension_trait_id"])
        mvp_source = mvp_by_trait.get(row["extension_trait_id"])
        if (source is not None and row["sleep_trait"] != "sleep_apnea") or mvp_source is not None:
            if mvp_source is not None:
                source = mvp_source
                code = source["study_accession"]
                source_id = f"gwas_catalog_{code}"
                url = source["source_url"]
                study_accession = code
                publication = "Diversity and scale: Genetic architecture of 2068 traits in the VA Million Veteran Program"
                pmid, doi = "39024449", "10.1126/science.adj1182"
                generation = "GWAS_CATALOG_ACCESSION_AND_MD5"
                storage_mode = "CHECKSUM_PINNED_REMOTE_STREAMING"
                full_resolution = "YES_CHECKSUM_PINNED_REMOTE"
                cohort_evidence = "Discovery external phenotype is Pan-UK Biobank; MVP is an independent U.S. Department of Veterans Affairs cohort."
            else:
                assert source is not None
                code = source["finngen_phenocode"]
                source_id = f"finngen_r13_{code}"
                url = f"https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_{code}.gz?generation={source['source_generation']}"
                study_accession = f"FINNGEN_R13_{code}"
                publication = "FinnGen flagship study and Release 13 public summary statistics"
                pmid, doi = "36653562", "10.1038/s41586-022-05473-8"
                generation = source["source_generation"]
                storage_mode = "VERSIONED_REMOTE_STREAMING"
                full_resolution = "YES_VERSIONED_REMOTE"
                cohort_evidence = "Discovery external phenotype is Pan-UK Biobank; FinnGen contains Finnish biobank participants and no UK Biobank participants."
            row.update({
                "replication_source_id": source_id,
                "replication_study_accession": study_accession,
                "replication_publication": publication,
                "replication_PMID": pmid, "replication_DOI": doi,
                "replication_source_url": url, "replication_checksum": f"md5:{source['md5_hex']}",
                "replication_source_generation": generation,
                "replication_etag": source["etag"],
                "replication_content_length_bytes": source["content_length_bytes"],
                "replication_storage_mode": storage_mode,
                "replication_local_path": "NOT_RETAINED_VERSIONED_REMOTE_STREAMING",
                "replication_receipt_path": f"discovery_extension/provenance/replication_streaming_receipts/{source_id}.json",
                "replication_munged_path": f"discovery_extension/data/replication/munged/{source_id}.sumstats.gz",
                "replication_phenotype_definition": source["replication_phenotype_definition"],
                "phenotype_match_status": source["phenotype_match_status"], "ancestry": "EUR", "build": "GRCh38",
                "sample_size": source["sample_size"], "cases": source["cases"], "controls": source["controls"],
                "discovery_cohort_relation": "NON_UKB", "participant_overlap_status": "NON_OVERLAPPING_CONFIRMED",
                "participant_overlap_evidence": cohort_evidence,
                "source_identity_status": "VERIFIED", "schema_status": "VERIFIED",
                "effect_allele_status": "UNAMBIGUOUS", "full_resolution_availability": full_resolution,
                "unavailable_reason": "", "source_curation_status": "COMPLETE_BEFORE_RESULTS",
            })
            testable += 1
        else:
            if source is not None:
                reason = "FinnGen R13 external endpoint rejected because the discovery sleep-apnea GWAS is FinnGen R9; participant overlap is non-negligible and violates the independent-replication contract."
                overlap_rejected += 1
            else:
                reason = "Completed GWAS Catalog exact-trait/synonym and FinnGen R13 searches found no non-UKB EUR source with a comparable definition, unambiguous effects, and downloadable full-resolution summary statistics. Controlled-access or unindexed sources may still exist."
                unavailable += 1
            for field in fields:
                if field.startswith("replication_") and field not in {
                    "replication_search_databases", "replication_search_queries", "replication_search_date",
                    "replication_search_evidence", "replication_family_size", "replication_alpha",
                }:
                    row[field] = "NOT_AVAILABLE"
            row.update({
                "unavailable_reason": reason,
                "source_curation_status": "NO_INDEPENDENT_DATASET_COMPLETE_BEFORE_RESULTS",
            })
        output.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(output)
    provenance = {
        "schema_version": "1.0.0", "completed_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": "before_replication_result_access", "replication_results_accessed": False,
        "candidate_family_pair_count": len(output), "testable_pair_count": testable,
        "participant_overlap_rejected_pair_count": overlap_rejected,
        "no_eligible_source_pair_count": unavailable, "bonferroni_alpha": lock["bonferroni_alpha"],
        "candidate_lock_sha256": sha256(args.candidate_lock), "source_config_sha256": sha256(args.sources),
        "mvp_source_config_sha256": sha256(args.mvp_sources),
        "search_evidence_sha256": sha256(args.search), "curated_queue_sha256": sha256(args.out),
        "replacement_after_results_allowed": False,
    }
    args.provenance_out.parent.mkdir(parents=True, exist_ok=True)
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"REPLICATION_QUEUE_CURATED pairs={len(output)} testable={testable} overlap_rejected={overlap_rejected} unavailable_no_source={unavailable} sha256={sha256(args.out)}")


if __name__ == "__main__":
    main()
