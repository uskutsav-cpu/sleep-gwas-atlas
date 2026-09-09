#!/usr/bin/env python3
"""Audit independent replication sources and publish related-phenotype LDSC checks."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


EXPECTED = [
    "longsleep__parental_lifespan",
    "shortsleep__parental_lifespan",
    "sleep_apnea__healthspan",
    "sleep_apnea__parental_lifespan",
    "sleep_efficiency__frailty",
    "snoring__healthspan",
]
AS_OF_DATE = "2026-08-28"


SOURCE_TEMPLATES = {
    "parental_lifespan": [
        {
            "source_id": "timmers_2019_parental_lifespan", "external_phenotype": "Parental survival/lifespan",
            "phenotype_relation": "EXACT", "ancestry": "EUR", "sample_size": "1012240 parents", "pmid": "30642433",
            "doi": "10.7554/eLife.39856", "source_url": "https://www.ebi.ac.uk/gwas/studies/GCST009890",
            "summary_statistics_status": "PUBLIC_FULL_SUMMARY_STATISTICS", "independence_from_discovery_external": "IDENTICAL_DISCOVERY_DATASET",
            "decision": "REJECT_PRIMARY_REPLICATION", "reason": "This is the discovery external GWAS.",
        },
        {
            "source_id": "joshi_2017_lifegen_parental_survival", "external_phenotype": "Parental survival/lifespan",
            "phenotype_relation": "EXACT_OR_NEAR_EXACT", "ancestry": "EUR", "sample_size": "606059 parents", "pmid": "29030599",
            "doi": "10.1038/s41467-017-00934-5", "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5715013/",
            "summary_statistics_status": "PUBLIC", "independence_from_discovery_external": "ANCESTOR_SUBSET_OF_TIMMERS_META_ANALYSIS",
            "decision": "REJECT_PRIMARY_REPLICATION", "reason": "The LifeGen/UKB data contribute to the later Timmers discovery meta-analysis.",
        },
        {
            "source_id": "shadyab_2017_hrs_parental_lifespan", "external_phenotype": "At least one long-lived parent",
            "phenotype_relation": "RELATED_DICHOTOMOUS_PARENTAL_LONGEVITY", "ancestry": "EUR_AFR", "sample_size": "HRS; small discovery",
            "pmid": "27816938", "doi": "10.1093/gerona/glw206", "source_url": "https://pubmed.ncbi.nlm.nih.gov/27816938/",
            "summary_statistics_status": "NO_PUBLIC_FULL_SUMMARY_STATISTICS_IDENTIFIED", "independence_from_discovery_external": "POTENTIALLY_INDEPENDENT_BUT_HRS_MAY_ENTER_OTHER_LONGEVITY_META_ANALYSES",
            "decision": "REJECT_PRIMARY_REPLICATION", "reason": "No analysis-ready genome-wide summary statistics were located; the reported lead did not replicate.",
        },
        {
            "source_id": "kaplanis_2019_ancestrydna_parental_lifespan", "external_phenotype": "Parental lifespan",
            "phenotype_relation": "EXACT_OR_NEAR_EXACT", "ancestry": "EUR", "sample_size": ">300000 discovery; >650000 meta-analysis",
            "pmid": "31484785", "doi": "10.1534/g3.119.400448", "source_url": "https://pubmed.ncbi.nlm.nih.gov/31484785/",
            "summary_statistics_status": "NO_PUBLIC_FULL_SUMMARY_STATISTICS_IDENTIFIED", "independence_from_discovery_external": "ANCESTRYDNA_DISCOVERY_INDEPENDENT;UKB_META_COMPONENT_OVERLAPS",
            "decision": "REJECT_PRIMARY_REPLICATION", "reason": "Genome-wide summary statistics were not publicly available for a reproducible LDSC test.",
        },
        {
            "source_id": "deelen_2019_longevity_90th", "external_phenotype": "Survival beyond cohort-specific 90th percentile",
            "phenotype_relation": "RELATED_EXTREME_LONGEVITY", "ancestry": "EUR", "sample_size": "11262 cases;25483 controls",
            "pmid": "31413261", "doi": "10.1038/s41467-019-11558-2", "source_url": "https://www.ebi.ac.uk/gwas/studies/GCST008598",
            "summary_statistics_status": "PUBLIC_FULL_SUMMARY_STATISTICS_LOCAL_HARMONIZED", "independence_from_discovery_external": "COHORT_OVERLAP_WITH_TIMMERS_POSSIBLE;NOT_AN_EXACT_TRAIT",
            "decision": "SELECT_RELATED_PHENOTYPE_SENSITIVITY", "reason": "Best available public aging endpoint, but not an exact independent parental-survival replication.",
        },
    ],
    "healthspan": [
        {
            "source_id": "zenin_2019_healthspan", "external_phenotype": "Healthspan termination",
            "phenotype_relation": "EXACT", "ancestry": "EUR_UKB", "sample_size": "300447", "pmid": "30729179",
            "doi": "10.1038/s42003-019-0290-0", "source_url": "https://www.ebi.ac.uk/gwas/studies/GCST007406",
            "summary_statistics_status": "PUBLIC_FULL_SUMMARY_STATISTICS", "independence_from_discovery_external": "IDENTICAL_DISCOVERY_DATASET",
            "decision": "REJECT_PRIMARY_REPLICATION", "reason": "This is the discovery external GWAS.",
        },
        {
            "source_id": "no_second_exact_healthspan_gwas", "external_phenotype": "Healthspan termination",
            "phenotype_relation": "EXACT", "ancestry": "NA", "sample_size": "NA", "pmid": "NA", "doi": "NA",
            "source_url": "https://www.ebi.ac.uk/gwas/", "summary_statistics_status": "NO_SECOND_EXACT_PUBLIC_GWAS_LOCATED",
            "independence_from_discovery_external": "NOT_APPLICABLE", "decision": "NO_VALID_DATASET",
            "reason": "Catalog, publication, and targeted web searches identified only the Zenin UKB genome-wide release for the exact endpoint.",
        },
        {
            "source_id": "timmers_2019_parental_lifespan", "external_phenotype": "Parental survival/lifespan",
            "phenotype_relation": "RELATED_AGING_ENDPOINT", "ancestry": "EUR", "sample_size": "1012240 parents", "pmid": "30642433",
            "doi": "10.7554/eLife.39856", "source_url": "https://www.ebi.ac.uk/gwas/studies/GCST009890",
            "summary_statistics_status": "PUBLIC_FULL_SUMMARY_STATISTICS_LOCAL_HARMONIZED", "independence_from_discovery_external": "RELATED_ENDPOINT_WITH_SUBSTANTIAL_UKB_COMPONENT",
            "decision": "SELECT_RELATED_PHENOTYPE_SENSITIVITY", "reason": "Related lifespan endpoint; not an exact healthspan replication.",
        },
        {
            "source_id": "deelen_2019_longevity_90th", "external_phenotype": "Survival beyond cohort-specific 90th percentile",
            "phenotype_relation": "RELATED_EXTREME_LONGEVITY", "ancestry": "EUR", "sample_size": "11262 cases;25483 controls",
            "pmid": "31413261", "doi": "10.1038/s41467-019-11558-2", "source_url": "https://www.ebi.ac.uk/gwas/studies/GCST008598",
            "summary_statistics_status": "PUBLIC_FULL_SUMMARY_STATISTICS_LOCAL_HARMONIZED", "independence_from_discovery_external": "DIFFERENT_ENDPOINT_AND_COHORT_SET;SOME_META_COHORT_OVERLAP_POSSIBLE",
            "decision": "SELECT_RELATED_PHENOTYPE_SENSITIVITY", "reason": "Related aging endpoint; not an exact healthspan replication.",
        },
    ],
    "frailty": [
        {
            "source_id": "atkins_2021_frailty_index", "external_phenotype": "Continuous frailty index",
            "phenotype_relation": "EXACT", "ancestry": "EUR_UKB_TWINGENE", "sample_size": "175226", "pmid": "34431594",
            "doi": "10.1111/acel.13459", "source_url": "https://www.ebi.ac.uk/gwas/studies/GCST90020053",
            "summary_statistics_status": "PUBLIC_FULL_SUMMARY_STATISTICS", "independence_from_discovery_external": "IDENTICAL_DISCOVERY_DATASET",
            "decision": "REJECT_PRIMARY_REPLICATION", "reason": "This is the discovery external GWAS.",
        },
        {
            "source_id": "kwak_2026_korean_frailty", "external_phenotype": "FRAIL scale, CHS index, and frailty index",
            "phenotype_relation": "EXACT_OR_RELATED_ALTERNATE_ANCESTRY", "ancestry": "KOR_EAS", "sample_size": "14664", "pmid": "41733842",
            "doi": "10.1007/s11357-026-02154-6", "source_url": "https://pubmed.ncbi.nlm.nih.gov/41733842/",
            "summary_statistics_status": "NO_PUBLIC_FULL_SUMMARY_STATISTICS_IDENTIFIED", "independence_from_discovery_external": "INDEPENDENT",
            "decision": "REJECT_PRIMARY_REPLICATION", "reason": "No public analysis-ready summary statistics were located; study is small and reports suggestive rather than genome-wide signals.",
        },
        {
            "source_id": "neale_2018_left_grip_strength", "external_phenotype": "Left-hand grip strength",
            "phenotype_relation": "RELATED_PHYSICAL_FUNCTION", "ancestry": "EUR_UKB", "sample_size": "359704", "pmid": "NA", "doi": "NA",
            "source_url": "http://www.nealelab.is/uk-biobank", "summary_statistics_status": "PUBLIC_FULL_SUMMARY_STATISTICS_LOCAL_HARMONIZED",
            "independence_from_discovery_external": "DIFFERENT_ENDPOINT_BUT_OVERLAPS_UKB_SLEEP_SAMPLE",
            "decision": "SELECT_RELATED_PHENOTYPE_SENSITIVITY", "reason": "Available physical-function proxy, explicitly non-independent and not frailty.",
        },
    ],
}


ALTERNATES = {
    "longsleep__parental_lifespan": ["longevity"],
    "shortsleep__parental_lifespan": ["longevity"],
    "sleep_apnea__healthspan": ["parental_lifespan", "longevity"],
    "sleep_apnea__parental_lifespan": ["longevity"],
    "sleep_efficiency__frailty": ["grip_strength"],
    "snoring__healthspan": ["parental_lifespan", "longevity"],
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def table_text(fields: list[str], rows: list[dict[str, object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "NA") for field in fields})
    return buffer.getvalue()


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def bh_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [1.0] * len(values)
    running = 1.0
    size = len(values)
    for rank_index in range(size - 1, -1, -1):
        index = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[index] * size / rank)
        result[index] = min(1.0, running)
    return result


def conceptual_concordance(discovery_external: str, discovery_rg: float, alternate: str, alternate_rg: float) -> bool:
    # Healthspan is coded toward termination/adverse outcome; lifespan, longevity, and grip are favorable.
    discovery_favorable_rg = -discovery_rg if discovery_external in {"healthspan", "frailty"} else discovery_rg
    alternate_favorable_rg = alternate_rg
    return (discovery_favorable_rg > 0) == (alternate_favorable_rg > 0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest_path = root / "results/validation/candidate_validation_manifest.tsv"
    manifest = read_tsv(manifest_path)
    if [row["candidate_id"] for row in manifest] != EXPECTED:
        fail("frozen candidate family/order drifted")
    manifest_by_id = {row["candidate_id"]: row for row in manifest}

    source_rows: list[dict[str, object]] = []
    for candidate_id in EXPECTED:
        external = manifest_by_id[candidate_id]["external_trait"]
        for source in SOURCE_TEMPLATES[external]:
            source_rows.append({"candidate_id": candidate_id, "discovery_external_trait": external, **source, "checked_date": AS_OF_DATE})

    rg_path = root / "results/tables/rg_matrix.tsv"
    rg_rows = read_tsv(rg_path)
    rg_index = {(row["sleep_trait"], row["disease_trait"]): row for row in rg_rows}
    replication: list[dict[str, object]] = []
    for candidate_id in EXPECTED:
        candidate = manifest_by_id[candidate_id]
        replication.append({
            "candidate_id": candidate_id, "analysis_tier": "PRIMARY_EXACT_INDEPENDENT", "replication_status": "NO_VALID_DATASET",
            "sleep_trait": candidate["sleep_trait"], "external_trait_tested": candidate["external_trait"],
            "external_source_id": "NA", "phenotype_relation": "EXACT", "ancestry": candidate["ancestry"],
            "sample_size": "NA", "independence_assessment": "NO_SECOND_EXACT_NONOVERLAPPING_PUBLIC_GWAS",
            "harmonization": "NOT_RUN", "method": "LDSC_NOT_RUN", "rg": "NA", "se": "NA", "p": "NA", "fdr_within_validation_alternates": "NA",
            "discovery_rg": candidate["discovery_rg"], "conceptual_direction_concordant": "NA",
            "calculation_origin": "NA", "calculation_input_sha256": "NA",
            "limitation": "Exact replication dataset unavailable; candidate retained and not replaced.",
        })
        for alternate in ALTERNATES[candidate_id]:
            key = (candidate["sleep_trait"], alternate)
            if key not in rg_index:
                fail(f"missing real LDSC alternate result: {key}")
            row = rg_index[key]
            source_id = {
                "longevity": "deelen_2019_longevity_90th",
                "parental_lifespan": "timmers_2019_parental_lifespan",
                "grip_strength": "neale_2018_left_grip_strength",
            }[alternate]
            source = next(item for item in SOURCE_TEMPLATES[candidate["external_trait"]] if item["source_id"] == source_id)
            relation = source["phenotype_relation"]
            replication.append({
                "candidate_id": candidate_id, "analysis_tier": "RELATED_PHENOTYPE_SENSITIVITY",
                "replication_status": "REAL_LDSC_RESULT_AVAILABLE_NONEXACT",
                "sleep_trait": candidate["sleep_trait"], "external_trait_tested": alternate,
                "external_source_id": source_id, "phenotype_relation": relation, "ancestry": source["ancestry"],
                "sample_size": source["sample_size"], "independence_assessment": source["independence_from_discovery_external"],
                "harmonization": "GRCh37;HapMap3;allele-checked;EUR 1000G LD;canonical atlas munging",
                "method": "LDSC_v1.0.1_CANONICAL_PHASE1_REAL_CALCULATION_REUSED", "rg": row["rg"], "se": row["se"], "p": row["p"],
                "fdr_within_validation_alternates": "PENDING", "discovery_rg": candidate["discovery_rg"],
                "conceptual_direction_concordant": str(conceptual_concordance(candidate["external_trait"], float(candidate["discovery_rg"]), alternate, float(row["rg"]))).upper(),
                "calculation_origin": "results/tables/rg_matrix.tsv", "calculation_input_sha256": sha256(rg_path),
                "limitation": source["reason"],
            })
    actual = [row for row in replication if row["analysis_tier"] == "RELATED_PHENOTYPE_SENSITIVITY"]
    adjusted = bh_adjust([float(row["p"]) for row in actual])
    for row, value in zip(actual, adjusted):
        row["fdr_within_validation_alternates"] = format(value, ".12g")

    source_fields = [
        "candidate_id", "discovery_external_trait", "source_id", "external_phenotype", "phenotype_relation", "ancestry",
        "sample_size", "pmid", "doi", "source_url", "summary_statistics_status", "independence_from_discovery_external",
        "decision", "reason", "checked_date",
    ]
    replication_fields = [
        "candidate_id", "analysis_tier", "replication_status", "sleep_trait", "external_trait_tested", "external_source_id",
        "phenotype_relation", "ancestry", "sample_size", "independence_assessment", "harmonization", "method", "rg", "se", "p",
        "fdr_within_validation_alternates", "discovery_rg", "conceptual_direction_concordant", "calculation_origin",
        "calculation_input_sha256", "limitation",
    ]
    outputs = {
        root / "results/validation/replication_source_search.tsv": table_text(source_fields, source_rows),
        root / "results/validation/replication_ldsc.tsv": table_text(replication_fields, replication),
    }
    if args.validate_only:
        for path, payload in outputs.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                fail(f"output absent or drifted: {path.relative_to(root)}")
    else:
        for path, payload in outputs.items():
            atomic_text(path, payload)
    print(f"REPLICATION_SEARCH_OK candidates={len(EXPECTED)} source_rows={len(source_rows)} exact_no_valid={len(EXPECTED)} related_ldsc={len(actual)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
