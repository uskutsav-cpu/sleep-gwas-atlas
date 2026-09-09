#!/usr/bin/env python3
"""Freeze the extension panel and its pre-analysis statistical contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


HIGH_NOVELTY_DOMAINS = {
    "renal_urinary",
    "pulmonary",
    "gastrointestinal",
    "sensory",
    "musculoskeletal_pain",
    "allergy_infectious",
    "medication_use",
    "oral_dental",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidate-pool",
        type=Path,
        default=Path("discovery_extension/results/candidate_pool.tsv"),
    )
    parser.add_argument(
        "--selection-seed",
        type=Path,
        default=Path("discovery_extension/config/panel_selection_seed.tsv"),
    )
    parser.add_argument(
        "--core-panel",
        type=Path,
        default=Path("config/analysis_panel.tsv"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("discovery_extension/config/candidate_traits.tsv"),
    )
    parser.add_argument(
        "--lock",
        type=Path,
        default=Path("discovery_extension/config/extension_panel.lock.json"),
    )
    args = parser.parse_args()

    pool = read_tsv(args.candidate_pool)
    seed = read_tsv(args.selection_seed)
    core = read_tsv(args.core_panel)
    pool_by_id = {row["extension_trait_id"]: row for row in pool}
    if len(pool_by_id) != len(pool):
        raise SystemExit("ERROR: duplicate ID in candidate pool")
    if not 75 <= len(seed) <= 125:
        raise SystemExit(f"ERROR: locked extension panel must contain 75-125 traits; found {len(seed)}")
    if len({row["extension_trait_id"] for row in seed}) != len(seed):
        raise SystemExit("ERROR: duplicate ID in selection seed")
    missing = [row["extension_trait_id"] for row in seed if row["extension_trait_id"] not in pool_by_id]
    if missing:
        raise SystemExit(f"ERROR: selected IDs absent from candidate pool: {missing}")

    fields = [
        "extension_trait_id",
        "phenotype_name",
        "phenotype_domain",
        "phenotype_definition",
        "source",
        "study_accession",
        "PMID",
        "DOI",
        "ancestry",
        "sample_size",
        "cases",
        "controls",
        "build",
        "source_url",
        "checksum",
        "source_tabix_url",
        "source_tabix_checksum",
        "binary_or_continuous",
        "available_beta",
        "available_se",
        "available_effect_allele",
        "available_other_allele",
        "available_frequency",
        "available_info",
        "prior_sleep_screen_coverage",
        "prior_direct_sleep_genetics_evidence",
        "novelty_priority",
        "selection_reason",
        "source_h2",
        "source_h2_se",
        "source_h2_z",
        "source_LDSC_intercept",
        "PanUKBB_QC_status",
        "in_max_independent_set",
        "source_filename",
        "source_file_size_bytes",
        "source_tabix_filename",
        "source_tabix_size_bytes",
    ]
    output: list[dict[str, str]] = []
    for selected in seed:
        source = pool_by_id[selected["extension_trait_id"]]
        novelty = source["novelty_priority"]
        if novelty == "UNDEREXPLORED" and selected["phenotype_domain"] in HIGH_NOVELTY_DOMAINS:
            novelty = "HIGH_NOVELTY_PRIORITY"
        output.append(
            {
                "extension_trait_id": source["extension_trait_id"],
                "phenotype_name": source["phenotype_name"],
                "phenotype_domain": selected["phenotype_domain"],
                "phenotype_definition": source["phenotype_definition"],
                "source": source["source"],
                "study_accession": source["study_accession"],
                "PMID": source["PMID"],
                "DOI": source["DOI"],
                "ancestry": source["ancestry"],
                "sample_size": source["sample_size"],
                "cases": source["cases"],
                "controls": source["controls"],
                "build": source["build"],
                "source_url": source["source_url"],
                "checksum": source["checksum"],
                "source_tabix_url": source["source_tabix_url"],
                "source_tabix_checksum": source["source_tabix_checksum"],
                "binary_or_continuous": source["binary_or_continuous"],
                "available_beta": source["available_beta"],
                "available_se": source["available_se"],
                "available_effect_allele": source["available_effect_allele"],
                "available_other_allele": source["available_other_allele"],
                "available_frequency": source["available_frequency"],
                "available_info": source["available_info"],
                "prior_sleep_screen_coverage": source["prior_sleep_screen_coverage"],
                "prior_direct_sleep_genetics_evidence": source[
                    "prior_direct_sleep_genetics_evidence"
                ],
                "novelty_priority": novelty,
                "selection_reason": (
                    f"{selected['selection_reason_code']}; passed pre-correlation Pan-UKB "
                    "sample, phenotype-QC, h2-Z, and source-intercept screens; selected "
                    "before any extension genetic correlation"
                ),
                "source_h2": source["h2"],
                "source_h2_se": source["h2_se"],
                "source_h2_z": source["h2_z"],
                "source_LDSC_intercept": source["LDSC_intercept"],
                "PanUKBB_QC_status": source["PanUKBB_QC_status"],
                "in_max_independent_set": source["in_max_independent_set"],
                "source_filename": source["filename"],
                "source_file_size_bytes": source["file_size_bytes"],
                "source_tabix_filename": source["filename_tabix"],
                "source_tabix_size_bytes": source["source_tabix_size_bytes"],
            }
        )
    write_tsv(args.out, fields, output)

    sleep_ids = [row["trait_id"] for row in core if row["domain"] == "sleep"]
    if len(sleep_ids) != 12:
        raise SystemExit(f"ERROR: core panel must expose exactly 12 sleep traits; found {len(sleep_ids)}")
    lock = {
        "schema_version": "1.0.0",
        "panel_name": "novelty-enriched-phenome-discovery-extension-v1",
        "locked_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "selection_timing": "before_any_extension_genetic_correlation",
        "no_result_driven_trait_selection": True,
        "extension_trait_count": len(output),
        "core_sleep_trait_count": len(sleep_ids),
        "planned_raw_rg_test_count": len(output) * len(sleep_ids),
        "extension_trait_ids_in_locked_order": [row["extension_trait_id"] for row in output],
        "core_sleep_trait_ids_in_locked_order": sleep_ids,
        "domain_counts": dict(sorted(Counter(row["phenotype_domain"] for row in output).items())),
        "novelty_priority_counts": dict(
            sorted(Counter(row["novelty_priority"] for row in output).items())
        ),
        "primary_h2_gate": {"h2_z_minimum": 4.0, "LDSC_intercept_maximum": 1.2},
        "primary_rg_gate": {
            "scope": "12 locked core sleep traits crossed with extension traits passing rerun h2 QC",
            "multiple_testing": "Benjamini-Hochberg FDR within extension family only",
            "fdr_threshold": 0.05,
            "effect_size_caution_threshold_abs_rg": 0.15,
            "qc_failed_traits": "excluded from primary family and retained in sensitivity outputs",
        },
        "prioritization_tiers": {
            "A": "extension FDR<0.05, abs(rg)>=0.15, both h2 pass, no major QC flag, pair-level novelty Strong or Moderate",
            "B": "replicates or is strongly reinforced by local/fine-mapping/colocalization evidence",
            "C": "nominal/suggestive, QC-sensitive, weakly replicated, or literature-ambiguous",
        },
        "immutability": {
            "no_trait_replacement_after_lock": True,
            "failed_h2_traits_remain_documented": True,
            "core_panel_unchanged": True,
            "extension_fdr_never_merged_with_core_396_test_family": True,
        },
        "artifact_hashes": {
            str(args.core_panel): sha256(args.core_panel),
            str(args.candidate_pool): sha256(args.candidate_pool),
            str(args.selection_seed): sha256(args.selection_seed),
            str(args.out): sha256(args.out),
        },
    }
    args.lock.parent.mkdir(parents=True, exist_ok=True)
    args.lock.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n")
    print(
        f"EXTENSION_PANEL_LOCKED traits={len(output)} planned_rg={len(output) * len(sleep_ids)} "
        f"candidate_traits_sha256={sha256(args.out)}"
    )


if __name__ == "__main__":
    main()
