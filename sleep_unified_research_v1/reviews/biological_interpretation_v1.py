#!/usr/bin/env python3
"""Independent, lightweight interpretation-evidence inventory; no GWAS estimation.

Run from any directory. Writes only biological_interpretation_v1.* reviewer files.
Table consistency checks do not validate biological or clinical interpretations.
"""
import csv
import hashlib
import itertools
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OLD = ROOT / "discovery_extension/sleep_submission_evidence_v1"
INPUTS = [
    ROOT / "sleep_unified_research_v1/FROZEN_NEW_ANALYSIS_PROTOCOL.md",
    HERE / "statistical_validity_v1.md",
    HERE / "numerical_reproduction_v1.md",
    HERE / "novelty_v1_report.md",
    HERE / "provenance_v1_adversarial_review.md",
    HERE / "provenance_v1_mvp_insomnia_priority.md",
    HERE / "provenance_v1_mvp_insomnia_source_contract.json",
    HERE / "provenance_v1_mvp_insomnia_217_admission.tsv",
    HERE / "sleep_phenotyping_v1.md",
    ROOT / "docs/genomic_sem_validation_results.md",
    ROOT / "docs/full_covariance_results.md",
    ROOT / "discovery_extension/final_report.md",
    ROOT / "discovery_extension/analysis/frozen_extension_audit.md",
    ROOT / "discovery_extension/results/top_novel_discoveries.tsv",
    OLD / "08_SLEEP_PHENOTYPE_AND_CLINICAL_INTERPRETATION.md",
    OLD / "09_NEGATIVE_AND_BLOCKED_RESULTS.md",
    OLD / "10_ADVERSARIAL_SCIENTIFIC_REVIEW.md",
    OLD / "07_SENSITIVITY_RESULTS.md",
    OLD / "tables/replicated_23.tsv",
    OLD / "tables/replication_family_217.tsv",
    OLD / "tables/global_1200.tsv",
    OLD / "tables/full_claim_evidence_ledger.tsv",
    OLD / "tables/sleep_trait_metadata.tsv",
    OLD / "tables/sensitivity_results.tsv",
    ROOT / "sleep_unified_research_v1/native/core_pilot_v1/rg_insomnia__bmi.full_precision.json",
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def pearson(a, b):
    aa, bb = sum(a) / len(a), sum(b) / len(b)
    x, y = [z - aa for z in a], [z - bb for z in b]
    return sum(p * q for p, q in zip(x, y)) / math.sqrt(
        sum(p * p for p in x) * sum(q * q for q in y)
    )


def main():
    rep = rows(OLD / "tables/replicated_23.tsv")
    glob = rows(OLD / "tables/global_1200.tsv")
    ledger = rows(OLD / "tables/full_claim_evidence_ledger.tsv")
    tops = rows(ROOT / "discovery_extension/results/top_novel_discoveries.tsv")
    sens = rows(OLD / "tables/sensitivity_results.tsv")
    metadata = {x["trait_id"]: x for x in rows(OLD / "tables/sleep_trait_metadata.tsv")}
    assert len(rep) == 23 and len(glob) == len(ledger) == 1200 and len(tops) == 23
    assert all(x["sleep_GWAS_reused"] == "True" for x in rep)
    assert all(x["claim_type"] == "OBSERVED_GLOBAL_GENETIC_CORRELATION" for x in ledger)
    assert all(x["allowed_interpretation"] ==
               "Global association only; no causal, local, molecular or clinical benefit inference"
               for x in ledger)
    mechanistic = ["local_loci", "pleiotropic_loci", "colocalized_genes", "strongest_mechanistic_evidence"]
    assert all(x[k] == "NA_BLOCKED_UPSTREAM" for x in tops for k in mechanistic)
    sleep_ids = sorted({x["sleep_trait"] for x in glob})
    outcomes = sorted({x["extension_trait_id"] for x in rep})
    vectors = {(x["sleep_trait"], x["extension_trait_id"]): float(x["rg"]) for x in glob}
    independently_computed_profiles = {
        (a, b): pearson([vectors[(s, a)] for s in sleep_ids],
                        [vectors[(s, b)] for s in sleep_ids])
        for a, b in itertools.combinations(outcomes, 2)
    }
    historical_profiles = {
        tuple(x["scope"].split(" versus ")): float(x["value"])
        for x in sens if x["analysis_id"] == "SLEEP_PROFILE_SIMILARITY"
    }
    assert len(historical_profiles) == len(independently_computed_profiles) == 28
    maximum_error = max(abs(v - historical_profiles[k])
                        for k, v in independently_computed_profiles.items())
    assert maximum_error < 1e-12
    stratum = {}
    for measurement in sorted({x["sleep_measurement"] for x in glob}):
        rr = [x for x in glob if x["sleep_measurement"] == measurement]
        stratum[measurement] = {
            "tests": len(rr),
            "BH_positive": sum(x["extension_fdr_pass_0.05"] == "True" for x in rr),
        }
    pilot = json.loads(INPUTS[-1].read_text())
    profile_values = list(independently_computed_profiles.values())
    report = {
        "review_date_client": "2026-10-08",
        "client_timezone": "America/Chicago",
        "inventory_generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Independent biological interpretation inventory of the frozen package; no new GWAS outcome inspection or estimation",
        "status": "GLOBAL_SHARING_ONLY; NO_LOCAL_MOLECULAR_CAUSAL_OR_CLINICAL_ACTIONABILITY_ADMITTED",
        "historical_positive_pairs": len(rep),
        "distinct_sleep_traits": len({x["sleep_trait"] for x in rep}),
        "distinct_external_outcome_labels": len(outcomes),
        "sleep_trait_positive_counts": dict(Counter(x["sleep_trait"] for x in rep)),
        "outcome_positive_counts": dict(Counter(x["external_phenotype_name"] for x in rep)),
        "historical_positive_sleep_GWAS_reused": sum(x["sleep_GWAS_reused"] == "True" for x in rep),
        "historical_two_trait_independent_positive_results": 0,
        "phenotype_match_labels": dict(Counter(x["phenotype_match_status"] for x in rep)),
        "all_positive_metadata_measurement_categories": dict(Counter(metadata[x["sleep_trait"]]["measurement_category"] for x in rep)),
        "discovery_measurement_strata_descriptive_only": stratum,
        "historical_ledger_claim_types": dict(Counter(x["claim_type"] for x in ledger)),
        "historical_ledger_allowed_interpretation": dict(Counter(x["allowed_interpretation"] for x in ledger)),
        "blocked_followup_counts": {k: dict(Counter(x[k] for x in tops)) for k in mechanistic},
        "outcome_profile_similarities": {
            "pairs": len(profile_values), "traits_per_vector": len(sleep_ids),
            "minimum": min(profile_values), "maximum": max(profile_values),
            "maximum_disagreement_with_prior_descriptive_table": maximum_error,
            "interpretation": "Pearson similarity of sleep-rg point-estimate vectors; not outcome genetic correlations, clusters or independent diseases",
        },
        "native_historical_pilot_boundary": {
            "pair": "insomnia__bmi", "rg_ratio": pilot["estimates"][0]["rg_ratio"],
            "purpose": "computational reproduction control; cannot identify BMI mediation or validate clinical insomnia",
        },
        "new_MVP_insomnia_pair_results_inspected": 0,
        "clinical_actionability": "NOT_TESTED; no risk prediction, severity, treatment or clinical decision endpoint in this audit",
        "biological_nulls": "No absence-of-mechanism inference from blocked analyses or threshold-negative rg",
        "primary_sources": [
            {"url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4797329/", "doi": "10.1038/ng.3406", "read": "primary model and scale sections", "supports": "global rg as standardized common-variant covariance; overlap model conditional on LDSC assumptions"},
            {"url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6451011/", "doi": "10.1038/s41467-019-09576-1", "read": "primary methods and limitation sections", "supports": "UKB device sampling, algorithm-derived sleep, time-window and quiet-wake limitations"},
            {"url": "https://www.med.upenn.edu/syspharmatt/assets/user-content/documents/science.adj1182.pdf", "doi": "10.1126/science.adj1182", "read": "primary PDF, page 10 automated-phenotyping limitations", "supports": "MVP diagnosis-code ascertainment, not clinical intervention utility"},
        ],
        "failed_primary_access": [
            {"url": "https://www.nature.com/articles/s41467-019-09576-1", "status": "web tool internal error; full PMC version successfully read"},
            {"url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10327290/", "status": "reCAPTCHA; exact encounter-rule details rely on separately documented provenance review"},
            {"url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12857194/", "status": "reCAPTCHA; authors' primary Science PDF successfully read for general phenotyping limits"},
        ],
        "input_hashes": [{"path": str(p), "bytes": p.stat().st_size, "sha256": digest(p)} for p in INPUTS],
        "reviewer_outputs": [{"path": str(p), "sha256": digest(p)} for p in
                             [Path(__file__), HERE / "biological_interpretation_v1.md", HERE / "biological_interpretation_v1.tsv"]],
        "warning": "Assertions certify internal evidence inventory only; they are not scientific validation or human review.",
    }
    (HERE / "biological_interpretation_v1.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["historical_positive_pairs", "distinct_sleep_traits", "distinct_external_outcome_labels", "blocked_followup_counts", "outcome_profile_similarities", "status"]}, indent=2))


if __name__ == "__main__":
    main()
