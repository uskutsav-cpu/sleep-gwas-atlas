#!/usr/bin/env python3
"""Freeze independent-replication sources before accessing replication results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path


AS_OF_DATE = "2026-09-01"
PAIR_MANIFEST = Path("results/track_b/pair_manifest.tsv")
LITERATURE_LOCK = Path("results/track_b/02_literature_audit.lock.json")
OUT = Path("results/track_b/replication_source_audit.tsv")
REPORT = Path("results/track_b/REPLICATION_SOURCE_AUDIT.md")
LOCK = Path("results/track_b/replication_source.lock.json")

FIELDS = [
    "pair_id", "sleep_trait", "discovery_external_trait", "source_id",
    "external_phenotype", "phenotype_relation", "ancestry", "build",
    "sample_size", "cases", "controls", "publication", "PMID", "DOI",
    "source_url", "content_length_bytes", "source_generation", "etag",
    "md5_hex", "summary_statistics_status", "discovery_external_relation",
    "sleep_dataset_relation", "cohort_overlap_assessment", "independence_status",
    "power_status", "selection_decision", "reason", "checked_date",
    "results_accessed_before_selection",
]


def source(pair_id: str, sleep: str, external: str, source_id: str,
           phenotype: str, relation: str, ancestry: str, build: str,
           sample_size: str, cases: str, controls: str, publication: str,
           pmid: str, doi: str, url: str, size: str, generation: str,
           etag: str, md5: str, stats: str, discovery_relation: str,
           sleep_relation: str, overlap: str, independence: str, power: str,
           decision: str, reason: str) -> dict[str, str]:
    return dict(zip(FIELDS, [
        pair_id, sleep, external, source_id, phenotype, relation, ancestry, build,
        sample_size, cases, controls, publication, pmid, doi, url, size,
        generation, etag, md5, stats, discovery_relation, sleep_relation,
        overlap, independence, power, decision, reason, AS_OF_DATE, "FALSE",
    ]))


def rows() -> list[dict[str, str]]:
    return [
        source(
            "A", "snoring", "parental_lifespan", "timmers_2019_parental_lifespan",
            "Parental survival/lifespan", "EXACT", "EUR", "GRCh37", "1012240 parents",
            "NA", "NA", "Timmers et al. 2019 eLife", "30642433", "10.7554/eLife.39856",
            "https://www.ebi.ac.uk/gwas/studies/GCST009890", "NA", "NA", "NA", "NA",
            "PUBLIC_FULL_SUMMARY_STATISTICS", "IDENTICAL_DISCOVERY_DATASET",
            "DISCOVERY_SNORING_GWAS", "SUBSTANTIAL_UKB_OVERLAP", "NOT_INDEPENDENT",
            "HIGH", "REJECT_PRIMARY_REPLICATION", "The external GWAS is the Pair A discovery source.",
        ),
        source(
            "A", "snoring", "parental_lifespan", "joshi_2017_lifegen_parental_survival",
            "Parental survival/lifespan", "EXACT_OR_NEAR_EXACT", "EUR", "GRCh37",
            "606059 parents", "NA", "NA", "Joshi et al. 2017 Nature Communications",
            "29030599", "10.1038/s41467-017-00934-5",
            "https://pmc.ncbi.nlm.nih.gov/articles/PMC5715013/", "NA", "NA", "NA", "NA",
            "PUBLIC_RELEASE_EXISTS", "ANCESTOR_COMPONENT_OF_TIMMERS_META_ANALYSIS",
            "DISCOVERY_SNORING_GWAS", "UKB_AND_LIFESPAN_SAMPLE_OVERLAP",
            "NOT_INDEPENDENT", "MODERATE_TO_HIGH", "REJECT_PRIMARY_REPLICATION",
            "LifeGen/UKB data contribute to the later discovery meta-analysis.",
        ),
        source(
            "A", "snoring", "parental_lifespan", "tanaka_2017_hrs_parental_lifespan",
            "At least one long-lived parent", "RELATED_DICHOTOMOUS_PARENTAL_LONGEVITY",
            "EUR_AND_AFR", "GRCh37_OR_STUDY_SPECIFIC", "SMALL_HRS_DISCOVERY", "NA", "NA",
            "Tanaka et al. 2017 Journals of Gerontology A", "27816938",
            "10.1093/gerona/glw206", "https://pubmed.ncbi.nlm.nih.gov/27816938/",
            "NA", "NA", "NA", "NA", "NO_PUBLIC_DENSE_SIGNED_SUMMARY_STATISTICS_IDENTIFIED",
            "POTENTIALLY_INDEPENDENT_RELATED_PHENOTYPE", "DISCOVERY_SNORING_GWAS",
            "NO_DIRECT_OVERLAP_IDENTIFIED_BUT_NOT_INDIVIDUAL_LEVEL_VERIFIED",
            "POTENTIALLY_INDEPENDENT", "LOW", "REJECT_PRIMARY_REPLICATION",
            "Phenotype differs, power is low, and no analysis-ready public full statistics were found.",
        ),
        source(
            "A", "snoring", "parental_lifespan", "wright_2019_ancestrydna_parental_lifespan",
            "Parental lifespan", "EXACT_OR_NEAR_EXACT", "EUR", "STUDY_SPECIFIC",
            ">300000 AncestryDNA discovery; >650000 meta-analysis", "NA", "NA",
            "Wright et al. 2019 G3", "31484785", "10.1534/g3.119.400448",
            "https://pubmed.ncbi.nlm.nih.gov/31484785/", "NA", "NA", "NA", "NA",
            "NO_PUBLIC_DENSE_SIGNED_ANCESTRYDNA_ONLY_STATISTICS_IDENTIFIED",
            "ANCESTRYDNA_COMPONENT_INDEPENDENT;META_ANALYSIS_INCLUDES_UKB",
            "DISCOVERY_SNORING_GWAS", "ANCESTRYDNA_COMPONENT_NON_UKB;META_ANALYSIS_OVERLAPS_UKB",
            "INDEPENDENT_COMPONENT_NOT_SEPARATELY_AVAILABLE", "HIGH_IF_COMPONENT_AVAILABLE",
            "REJECT_PRIMARY_REPLICATION", "The independent component is not publicly available for reproducible LDSC; the public paper's meta-analysis includes UKB.",
        ),
        source(
            "A", "snoring", "parental_lifespan", "deelen_2019_longevity_90th",
            "Survival beyond cohort-specific 90th percentile", "RELATED_EXTREME_LONGEVITY",
            "EUR", "GRCh37", "36745", "11262", "25483",
            "Deelen et al. 2019 Nature Communications", "31413261",
            "10.1038/s41467-019-11558-2", "https://www.ebi.ac.uk/gwas/studies/GCST008598",
            "NA", "NA", "NA", "NA", "PUBLIC_FULL_SUMMARY_STATISTICS",
            "DIFFERENT_AGING_ENDPOINT_WITH_POSSIBLE_CONTRIBUTING_COHORT_OVERLAP",
            "DISCOVERY_SNORING_GWAS", "POSSIBLE_LIFESPAN_META_COHORT_OVERLAP;NO_EXACT_UKB_IDENTITY",
            "UNCERTAIN_AND_NONEXACT", "LOW_TO_MODERATE", "RELATED_PHENOTYPE_SENSITIVITY_ONLY",
            "Useful only as a nonexact aging sensitivity; it cannot satisfy exact independent replication.",
        ),
        source(
            "A", "snoring", "parental_lifespan", "no_valid_exact_pair_a_replication",
            "Parental lifespan", "EXACT", "EUR", "NA", "NA", "NA", "NA",
            "Targeted source audit", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
            "NO_VALID_PUBLIC_ANALYSIS_READY_SOURCE", "NO_ELIGIBLE_SOURCE",
            "DISCOVERY_SNORING_GWAS", "NOT_APPLICABLE", "NO_VALID_REPLICATION",
            "NOT_APPLICABLE", "NO_VALID_REPLICATION", "No well-powered, exact, independent, public analysis-ready parental-lifespan GWAS was identified.",
        ),
        source(
            "B", "insomnia", "adhd", "demontis_2023_adhd_eur",
            "Clinical ADHD meta-analysis", "EXACT", "EUR", "GRCh37", "225534",
            "38691", "186843", "Demontis et al. 2023 Nature Genetics", "36702997",
            "10.1038/s41588-022-01285-8", "https://pubmed.ncbi.nlm.nih.gov/36702997/",
            "NA", "NA", "NA", "NA", "PUBLIC_FULL_SUMMARY_STATISTICS",
            "IDENTICAL_DISCOVERY_DATASET", "DISCOVERY_UKB_INSOMNIA_GWAS",
            "EXTERNAL_IDENTICAL;SLEEP_IDENTICAL", "NOT_INDEPENDENT", "HIGH",
            "REJECT_PRIMARY_REPLICATION", "The external GWAS is the Pair B discovery source.",
        ),
        source(
            "B", "insomnia", "adhd", "hernandez_2023_abcd_childhood",
            "Dimensional childhood ADHD symptoms", "RELATED_DIMENSIONAL_PEDIATRIC",
            "EUR_SUBSET", "INDIVIDUAL_LEVEL_STUDY", "4728 genetics", "NA", "NA",
            "Hernandez et al. 2023 Biological Psychiatry: Global Open Science", "36712562",
            "10.1016/j.bpsgos.2021.12.011", "https://pubmed.ncbi.nlm.nih.gov/36712562/",
            "NA", "NA", "NA", "NA", "NO_STANDARDIZED_PUBLIC_GWAS_PAIR_FOR_LDSC",
            "INDEPENDENT_BUT_DIFFERENT_AGE_AND_DIMENSIONAL_PHENOTYPE",
            "CHILDHOOD_INSOMNIA_DIMENSION_NOT_DISCOVERY_GWAS", "ABCD_NON_UKB_AND_NON_DISCOVERY_ADHD",
            "INDEPENDENT_BUT_NOT_COMPARABLE", "LOW", "RELATED_PHENOTYPE_CONTEXT_ONLY",
            "Small developmental covariance result cannot serve as the requested standardized adult exact-pair LDSC replication.",
        ),
        source(
            "B", "insomnia", "adhd", "finngen_r13_F5_ADHD",
            "Disturbance of activity and attention; ICD-10 F90.0 / ICD-9 3140",
            "CLOSELY_MATCHED_REGISTER_ADHD", "FINNISH_EUR", "GRCh38", "495052",
            "5559", "489493", "FinnGen Release 13", "36653562",
            "10.1038/s41586-022-05473-8",
            "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_ADHD.gz?generation=1777989549884404",
            "802195177", "1777989549884404", "e7aedc4602071840bbd5abe62d6fcdcc",
            "e7aedc4602071840bbd5abe62d6fcdcc", "PUBLIC_FULL_SIGNED_SUMMARY_STATISTICS_INFO_GT_0.6",
            "DIFFERENT_FINNISH_BIOBANK_COHORT_FROM_IPSYCH_DECODE_PGC_BY_PUBLISHED_COHORT_DESCRIPTIONS",
            "SAME_DISCOVERY_UKB_INSOMNIA_GWAS_REUSED_AS_ALLOWED_EXTERNAL_PHENOTYPE_REPLICATION_DESIGN",
            "NO_DIRECT_PARTICIPANT_OVERLAP_IDENTIFIED;NOT_INDIVIDUAL_LEVEL_VERIFIED",
            "INDEPENDENT_EXTERNAL_COHORT_WITH_DOCUMENTED_FOUNDER_POPULATION_CAVEAT",
            "UNDERPOWERED_RELATIVE_TO_DISCOVERY_BUT_ANALYZABLE", "SELECT_PRIMARY_REPLICATION",
            "Best public exact/comparable independent ADHD source; Finnish LD/ascertainment and low case count require cautious interpretation.",
        ),
    ]


def tsv_text(output_rows: list[dict[str, str]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_rows)
    return buffer.getvalue()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def validate_inputs() -> None:
    if not PAIR_MANIFEST.is_file() or not LITERATURE_LOCK.is_file():
        raise SystemExit("ERROR: pair manifest or literature lock missing")
    with PAIR_MANIFEST.open(newline="", encoding="utf-8") as handle:
        manifest = {r["pair_id"]: (r["sleep_trait"], r["external_trait"]) for r in csv.DictReader(handle, delimiter="\t")}
    if manifest.get("A") != ("snoring", "parental_lifespan") or manifest.get("B") != ("insomnia", "adhd"):
        raise SystemExit("ERROR: frozen discovery pairs drifted")


def report_text() -> str:
    return """# Track B independent-replication source audit

Frozen: 2026-09-01, before reading any new replication result.

## Pair A — snoring / parental lifespan

`NO_VALID_REPLICATION`. The exact public source is the discovery GWAS. LifeGen is an ancestor component of that meta-analysis. HRS is small, differently defined, and lacks an identified dense public release. The independent AncestryDNA component is not separately available as analysis-ready public summary statistics, while its combined meta-analysis includes UK Biobank. Extreme longevity is retained only as a nonexact sensitivity and cannot earn replication credit.

Mechanistic follow-up may proceed as an unreplicated discovery analysis, but Pair A cannot reach Track B Tier 2 until a truly independent, closely matched parental-lifespan dataset is obtained.

## Pair B — insomnia / ADHD

FinnGen R13 endpoint `F5_ADHD` is frozen as the primary external-phenotype replication source: 5,559 cases and 489,493 controls, Finnish-European, GRCh38, 802,195,177 compressed bytes. The version is pinned by GCS generation `1777989549884404` and MD5/ETag `e7aedc4602071840bbd5abe62d6fcdcc`.

This is independent by published cohort descriptions from the iPSYCH/deCODE/PGC discovery meta-analysis, not by individual-level linkage. It is much smaller in case count and uses Finnish register ascertainment, so a null result can be `UNDERPOWERED`; founder-population LD and phenotype differences must be carried into heterogeneity interpretation.

The UKB insomnia GWAS is intentionally held fixed. Track B asks for an independent external phenotype GWAS whenever feasible; replacing both GWAS would change the scientific question and reduce comparability.

## Frozen execution rule

- Use the same standardized LDSC/HapMap3 pipeline and thresholds as the atlas.
- Do not substitute a different ADHD endpoint after seeing the result.
- Record h2/QC before interpreting rg.
- A well-powered, materially incompatible opposite direction triggers the specified no-go gate.
- A low-power null is not promoted to failed replication.
"""


def build() -> tuple[str, str, str]:
    validate_inputs()
    output_rows = rows()
    if len({r["source_id"] for r in output_rows}) != len(output_rows):
        raise SystemExit("ERROR: duplicate source IDs")
    selected = [r for r in output_rows if r["selection_decision"] == "SELECT_PRIMARY_REPLICATION"]
    if [(r["pair_id"], r["source_id"]) for r in selected] != [("B", "finngen_r13_F5_ADHD")]:
        raise SystemExit("ERROR: selected replication source family drifted")
    table = tsv_text(output_rows)
    report = report_text()
    lock = {
        "schema_version": 1,
        "frozen_date": AS_OF_DATE,
        "selection_timing": "BEFORE_REPLICATION_RESULT_ACCESS",
        "replication_results_accessed_before_selection": False,
        "pair_manifest_sha256": hashlib.sha256(PAIR_MANIFEST.read_bytes()).hexdigest(),
        "literature_lock_sha256": hashlib.sha256(LITERATURE_LOCK.read_bytes()).hexdigest(),
        "source_audit_sha256": sha256_bytes(table.encode()),
        "report_sha256": sha256_bytes(report.encode()),
        "pair_a_primary_replication": "NO_VALID_REPLICATION",
        "pair_b_primary_replication_source": "finngen_r13_F5_ADHD",
        "pair_b_source_generation": "1777989549884404",
        "pair_b_source_md5": "e7aedc4602071840bbd5abe62d6fcdcc",
        "pair_b_source_content_length_bytes": 802195177,
        "post_result_source_replacement": "FORBIDDEN_WITHOUT_A_NEW_EXPLICIT_PROTOCOL_VERSION_AND_NONRESULT_BASED_JUSTIFICATION",
    }
    return table, report, json.dumps(lock, indent=2, sort_keys=True) + "\n"


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(value, encoding="utf-8")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    expected = build()
    paths = [OUT, REPORT, LOCK]
    if args.verify:
        for path, value in zip(paths, expected):
            if not path.is_file() or path.read_text(encoding="utf-8") != value:
                raise SystemExit(f"ERROR: replication source freeze missing or drifted: {path}")
        print("verified Track B replication source freeze: Pair A none; Pair B finngen_r13_F5_ADHD")
        return
    for path, value in zip(paths, expected):
        atomic_text(path, value)
    print("wrote Track B replication source freeze: Pair A none; Pair B finngen_r13_F5_ADHD")


if __name__ == "__main__":
    main()
