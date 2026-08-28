#!/usr/bin/env python3
"""Build the pair-complete Phase-1 literature audit and replication analyses."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path


PRIMARY = "PRIMARY_PHASE1"
DIRECT_CLASSES = {
    "DIRECT_RG_PREVIOUSLY_REPORTED",
    "DIRECT_RG_REPLICATION_DIFFERENT_DATASET",
}
ALLOWED_CLASSES = DIRECT_CLASSES | {
    "RELATED_GENETIC_EVIDENCE_ONLY",
    "MR_ONLY",
    "OBSERVATIONAL_ONLY",
    "NO_DIRECT_RG_FOUND",
    "APPARENTLY_NOVEL",
    "UNCERTAIN",
}
ADVERSE_SLEEP = {
    "insomnia", "shortsleep", "longsleep", "sleepiness", "napping", "snoring", "sleep_apnea",
}
FAVORABLE_EXTERNAL = {"hdl", "parental_lifespan", "longevity", "telomere_length", "grip_strength"}
ADVERSE_EXTERNAL = {
    "mdd", "scz", "bipolar", "adhd", "ibd", "crohn", "ra", "asthma", "bmi",
    "triglycerides", "cad", "stroke", "atrial_fibrillation", "frailty", "breast_cancer",
    "prostate_cancer", "colorectal_cancer", "lung_cancer", "ovarian_cancer", "parkinson",
}


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_text(rows: list[dict[str, object]], fields: list[str]) -> str:
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


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def fmt(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.12g}"


def semicolon(values: list[object] | set[object]) -> str:
    clean = sorted({str(value) for value in values if str(value) not in {"", "NA", "None"}})
    return ";".join(clean) if clean else "NA"


def pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) < 2 or len(x) != len(y):
        return None
    x_mean = statistics.fmean(x)
    y_mean = statistics.fmean(y)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y))
    denominator = math.sqrt(
        sum((a - x_mean) ** 2 for a in x) * sum((b - y_mean) ** 2 for b in y)
    )
    return numerator / denominator if denominator else None


def normalize_source_evidence(row: dict[str, str]) -> dict[str, str]:
    result = dict(row)
    result["evidence_location"] = (
        f"{row.get('source_table', 'NA')} row {row.get('source_row_number', 'NA')}"
    )
    result["published_external_sample_size"] = row.get("published_trait_sample_size", "NA")
    result["evidence_summary"] = "Direct rg row extracted from a reviewed checksum-pinned publication supplement."
    result["contradiction_note"] = "NONE"
    return result


def normalize_additional_evidence(row: dict[str, str]) -> dict[str, str]:
    result = dict(row)
    result["source_table"] = row.get("evidence_location", "NA")
    result["source_row_number"] = "NA"
    result["published_trait_pmid"] = "NA"
    result["published_trait_sample_size"] = row.get("published_external_sample_size", "NA")
    return result


def evidence_rank(row: dict[str, str]) -> tuple[int, int, int, int, str]:
    exact = 0 if row.get("external_phenotype_match") == "EXACT_CONCEPT" else 1
    numeric = 0 if number(row.get("published_rg")) is not None else 1
    level = 0 if any(token in row.get("extraction_status", "") for token in ("TABLE_LEVEL", "FULLTEXT")) else 1
    preprint = 1 if row.get("paper_doi", "").startswith("10.1101/") else 0
    return exact, numeric, level, preprint, row.get("evidence_source_id", "")


def direct_classification(evidence: list[dict[str, str]]) -> str:
    changed_dataset = any(
        "DIFFERENT_OR_EXPANDED" in row.get("sleep_dataset_relation", "")
        for row in evidence
    )
    return (
        "DIRECT_RG_REPLICATION_DIFFERENT_DATASET"
        if changed_dataset else "DIRECT_RG_PREVIOUSLY_REPORTED"
    )


def direct_confidence(evidence: list[dict[str, str]]) -> str:
    if any(row.get("external_phenotype_match") == "EXACT_CONCEPT" for row in evidence):
        return "HIGH"
    return "MODERATE"


def contradiction_text(atlas_rg: float, evidence: list[dict[str, str]]) -> str:
    notes = [
        row.get("contradiction_note", "")
        for row in evidence
        if row.get("contradiction_note", "") not in {"", "NA", "NONE"}
    ]
    for row in evidence:
        published = number(row.get("published_rg"))
        published_p = number(row.get("published_rg_p"))
        if published is None or published == 0 or atlas_rg == 0:
            continue
        if math.copysign(1, published) != math.copysign(1, atlas_rg):
            label = "PUBLISHED_DIRECTION_DISCORDANT" if published_p is not None and published_p < .05 else "PUBLISHED_NULL_OR_DIRECTION_UNSTABLE"
            notes.append(
                f"{label}: {row.get('evidence_source_id')} published_rg={published:.4g} atlas_rg={atlas_rg:.4g}"
            )
    return " | ".join(dict.fromkeys(notes)) if notes else "NONE_IDENTIFIED"


def observational_context(sleep: str, external: str) -> str:
    if sleep in ADVERSE_SLEEP and external in FAVORABLE_EXTERNAL:
        return "NEGATIVE_DIRECTION_BROADLY_CONSISTENT_WITH_ADVERSE_SLEEP_AND_FAVORABLE_HEALTH_TRAIT_LITERATURE"
    if sleep in {"sleep_efficiency", "sleepdur"} and external in ADVERSE_EXTERNAL:
        return "NEGATIVE_DIRECTION_BROADLY_CONSISTENT_WITH_HEALTHIER_SLEEP_AND_LOWER_DISEASE_BURDEN"
    if sleep == "chronotype" and external in ADVERSE_EXTERNAL:
        return "NEGATIVE_MORNINGNESS_DIRECTION_BROADLY_CONSISTENT_BUT_PHENOTYPE_CODING_REQUIRES_CARE"
    if sleep in ADVERSE_SLEEP and external in ADVERSE_EXTERNAL:
        return "NEGATIVE_DIRECTION_OPPOSES_OR_COMPLICATES_THE_USUAL_POSITIVE_COMORBIDITY_NARRATIVE"
    return "MIXED_OR_NO_CLEAR_OBSERVATIONAL_DIRECTION_EXPECTATION"


def replication_rows(
    significant: list[dict[str, str]], evidence: list[dict[str, str]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    master = {(row["sleep_trait"], row["external_trait"]): row for row in significant}
    output: list[dict[str, object]] = []
    for item in evidence:
        published = number(item.get("published_rg"))
        if published is None:
            continue
        key = (item["sleep_trait"], item["external_trait"])
        atlas = master[key]
        atlas_rg = float(atlas["rg"])
        atlas_se = float(atlas["se"])
        published_se = number(item.get("published_rg_se"))
        exact = item.get("external_phenotype_match") == "EXACT_CONCEPT"
        direction = (
            "CONCORDANT" if published == 0 or atlas_rg == 0 or math.copysign(1, published) == math.copysign(1, atlas_rg)
            else "DISCORDANT"
        )
        uncertainty = "NA"
        if published_se is not None and published_se > 0:
            uncertainty = str(
                abs(atlas_rg - published) <= 1.96 * math.sqrt(atlas_se ** 2 + published_se ** 2)
            )
        output.append({
            "sleep_trait": key[0],
            "external_trait": key[1],
            "external_domain": atlas["external_domain"],
            "atlas_rg": atlas_rg,
            "atlas_se": atlas_se,
            "atlas_p": atlas["p"],
            "atlas_fdr": atlas["fdr"],
            "evidence_source_id": item["evidence_source_id"],
            "paper_pmid": item.get("paper_pmid", "NA"),
            "paper_doi": item.get("paper_doi", "NA"),
            "article_title": item.get("article_title", "NA"),
            "evidence_location": item.get("evidence_location", item.get("source_table", "NA")),
            "published_trait_label": item.get("published_trait_label", "NA"),
            "external_phenotype_match": item.get("external_phenotype_match", "NA"),
            "published_rg": published,
            "published_rg_se": item.get("published_rg_se", "NA"),
            "published_rg_p": item.get("published_rg_p", "NA"),
            "atlas_minus_published_rg": atlas_rg - published,
            "absolute_difference": abs(atlas_rg - published),
            "direction_concordance": direction,
            "consistent_within_combined_95pct_uncertainty": uncertainty,
            "sleep_dataset_relation": item.get("sleep_dataset_relation", "NA"),
            "external_dataset_relation": item.get("external_dataset_relation", "NA"),
            "comparison_independence": item.get("comparison_independence", "NA"),
            "comparator_eligible_for_headline_metrics": str(exact),
            "source_url": item.get("source_url", "NA"),
            "source_sha256": item.get("source_sha256", "NA"),
            "contradiction_note": item.get("contradiction_note", "NONE"),
        })
    output.sort(key=lambda row: (
        row["sleep_trait"], row["external_trait"], row["evidence_source_id"], row["published_rg"]
    ))
    eligible = [row for row in output if row["comparator_eligible_for_headline_metrics"] == "True"]
    x = [float(row["published_rg"]) for row in eligible]
    y = [float(row["atlas_rg"]) for row in eligible]
    uncertainty_rows = [
        row for row in eligible
        if row["consistent_within_combined_95pct_uncertainty"] in {"True", "False"}
    ]
    pair_values: defaultdict[tuple[str, str], list[float]] = defaultdict(list)
    for row in eligible:
        pair_values[(str(row["sleep_trait"]), str(row["external_trait"]))].append(float(row["published_rg"]))
    pair_x, pair_y = [], []
    for key, values in sorted(pair_values.items()):
        pair_x.append(statistics.median(values))
        pair_y.append(float(master[key]["rg"]))
    metrics: dict[str, object] = {
        "all_numeric_direct_evidence_rows": len(output),
        "exact_concept_numeric_comparator_rows": len(eligible),
        "exact_concept_unique_pairs": len(pair_values),
        "row_level_pearson_atlas_vs_published": pearson(x, y),
        "row_level_mean_absolute_difference": statistics.fmean(abs(a - b) for a, b in zip(x, y)) if x else None,
        "row_level_direction_concordance": sum(a == b for a, b in zip([v >= 0 for v in x], [v >= 0 for v in y])) / len(x) if x else None,
        "uncertainty_comparable_row_count": len(uncertainty_rows),
        "proportion_consistent_within_combined_95pct_uncertainty": (
            sum(row["consistent_within_combined_95pct_uncertainty"] == "True" for row in uncertainty_rows) / len(uncertainty_rows)
            if uncertainty_rows else None
        ),
        "pair_level_median_published_pearson": pearson(pair_x, pair_y),
        "pair_level_median_published_mean_absolute_difference": (
            statistics.fmean(abs(a - b) for a, b in zip(pair_x, pair_y)) if pair_x else None
        ),
        "pair_level_median_published_direction_concordance": (
            sum(a == b for a, b in zip([v >= 0 for v in pair_x], [v >= 0 for v in pair_y])) / len(pair_x)
            if pair_x else None
        ),
        "interpretation_limit": "Most published comparisons reuse or overlap public GWAS inputs; these are pipeline positive controls, not independent biological replications.",
    }
    return output, metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    analysis = root / "results/analysis"

    master_all = read_tsv(analysis / "phase1_master_analysis.tsv")
    significant = [
        row for row in master_all
        if row["primary_or_sensitivity"] == PRIMARY and row["locked_primary_significant"] == "True"
    ]
    if len(significant) != 153:
        fail(f"expected 153 locked primary discoveries, observed {len(significant)}")
    significant_keys = {(row["sleep_trait"], row["external_trait"]) for row in significant}
    confidence = {
        (row["sleep_trait"], row["external_trait"]): row
        for row in read_tsv(analysis / "result_confidence.tsv")
    }

    source_evidence = [
        normalize_source_evidence(row)
        for row in read_tsv(analysis / "published_rg_source_evidence.tsv")
    ]
    additional_evidence = [
        normalize_additional_evidence(row)
        for row in read_tsv(root / "config/reviewed_additional_rg_evidence.tsv")
    ]
    evidence = source_evidence + additional_evidence
    if any((row["sleep_trait"], row["external_trait"]) not in significant_keys for row in evidence):
        fail("direct evidence contains a pair outside the 153 primary discoveries")
    evidence_by_pair: defaultdict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in evidence:
        evidence_by_pair[(row["sleep_trait"], row["external_trait"])].append(row)

    reviews = read_tsv(root / "config/literature_pair_reviews.tsv")
    review_by_pair = {(row["sleep_trait"], row["external_trait"]): row for row in reviews}
    if len(review_by_pair) != len(reviews):
        fail("manual literature review keys are duplicated")
    if set(evidence_by_pair) & set(review_by_pair):
        fail("manual non-direct reviews overlap direct-evidence pairs")
    if set(evidence_by_pair) | set(review_by_pair) != significant_keys:
        fail("direct evidence plus manual reviews do not cover exactly 153 discoveries")
    if any(row["novelty_classification"] not in ALLOWED_CLASSES - DIRECT_CLASSES for row in reviews):
        fail("manual literature review contains an invalid classification")

    queries = read_tsv(analysis / "literature_search_queries.tsv")
    candidates = read_tsv(analysis / "literature_search_candidates.tsv")
    query_by_pair: defaultdict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    candidate_by_pair: defaultdict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in queries:
        query_by_pair[(row["sleep_trait"], row["external_trait"])].append(row)
    for row in candidates:
        candidate_by_pair[(row["sleep_trait"], row["external_trait"])].append(row)
    if any(len(query_by_pair[key]) != 4 for key in significant_keys):
        fail("each significant pair must have exactly four complete literature queries")
    if any(row["query_status"] != "COMPLETE" for row in queries):
        fail("literature audit cannot proceed with incomplete queries")

    audit: list[dict[str, object]] = []
    for atlas in significant:
        key = (atlas["sleep_trait"], atlas["external_trait"])
        pair_evidence = evidence_by_pair.get(key, [])
        pair_queries = sorted(query_by_pair[key], key=lambda row: row["evidence_lens"])
        pair_candidates = candidate_by_pair.get(key, [])
        screened_pmids = {row["pmid"] for row in pair_candidates if row["pmid"] != "NA"}
        if pair_evidence:
            classification = direct_classification(pair_evidence)
            classification_confidence = direct_confidence(pair_evidence)
            representative = sorted(pair_evidence, key=evidence_rank)[0]
            in_depth_pmids = {row.get("paper_pmid", "NA") for row in pair_evidence}
            screened_pmids |= in_depth_pmids
            direct_search_result = "EXPLICIT_DIRECT_RG_EVIDENCE_REVIEWED"
            observation = observational_context(*key)
            contradiction = contradiction_text(float(atlas["rg"]), pair_evidence)
            replication_feasibility = "POSSIBLE_BUT_EXISTING_GWAS_OVERLAP_MUST_BE_MANAGED"
            resource_context = "DENSE_SUMMARY_STATISTICS_AVAILABLE_FOR_ATLAS_INPUTS"
            rationale = (
                f"{len(pair_evidence)} direct evidence row(s) were reviewed; classification distinguishes a changed or expanded dataset version from the same-source published analysis."
            )
            phenotype_match = semicolon({row.get("external_phenotype_match", "NA") for row in pair_evidence})
            evidence_pmids = semicolon({row.get("paper_pmid", "NA") for row in pair_evidence})
            evidence_dois = semicolon({row.get("paper_doi", "NA") for row in pair_evidence})
            evidence_sources = semicolon({row.get("evidence_source_id", "NA") for row in pair_evidence})
            rg_values = semicolon([
                f"{row.get('evidence_source_id')}={row.get('published_rg')}"
                for row in pair_evidence if number(row.get("published_rg")) is not None
            ])
            sample_values = semicolon([
                f"{row.get('evidence_source_id')}:sleepN={row.get('published_sleep_sample_size', 'NA')},externalN={row.get('published_external_sample_size', row.get('published_trait_sample_size', 'NA'))}"
                for row in pair_evidence
            ])
            dataset_relation = semicolon({row.get("sleep_dataset_relation", "NA") for row in pair_evidence})
            independence = semicolon({row.get("comparison_independence", "NA") for row in pair_evidence})
            review_status = "REVIEWED_DIRECT_EVIDENCE"
        else:
            review = review_by_pair[key]
            classification = review["novelty_classification"]
            classification_confidence = review["classification_confidence"]
            representative = {
                "paper_pmid": review["representative_pmid"],
                "paper_doi": review["representative_doi"],
                "article_title": review["representative_title"],
                "published_rg": "NA", "published_rg_se": "NA", "published_rg_p": "NA",
                "published_external_sample_size": "NA", "published_sleep_sample_size": "NA",
            }
            in_depth_pmids = set(review["papers_inspected_pmids"].split(";"))
            screened_pmids |= in_depth_pmids
            direct_search_result = review["direct_rg_search_result"]
            observation = review["observational_direction_context"]
            contradiction = review["contradictory_literature"]
            replication_feasibility = review["replication_feasibility"]
            resource_context = review["resource_context"]
            rationale = review["review_rationale"]
            phenotype_match = review["phenotype_match"]
            evidence_pmids = review["representative_pmid"]
            evidence_dois = review["representative_doi"]
            evidence_sources = "NA"
            rg_values = "NA"
            sample_values = "NA"
            dataset_relation = "NA"
            independence = "NA"
            review_status = review["review_status"]
        lens_queries = {row["evidence_lens"]: row for row in pair_queries}
        numeric_direct = sum(number(row.get("published_rg")) is not None for row in pair_evidence)
        conf = confidence[key]
        audit.append({
            "pair_id": f"{key[0]}__{key[1]}",
            "sleep_trait": key[0], "sleep_trait_label": atlas["sleep_trait_label"],
            "external_trait": key[1], "external_trait_label": atlas["external_trait_label"],
            "external_domain": atlas["external_domain"],
            "rg": atlas["rg"], "se": atlas["se"], "p": atlas["p"], "fdr": atlas["fdr"],
            "abs_rg": atlas["abs_rg"], "direction": atlas["direction"],
            "result_confidence_class": conf["confidence_class"],
            "measurement_category": atlas["sleep_measurement_category"],
            "novelty_classification": classification,
            "classification_confidence": classification_confidence,
            "direct_rg_found": str(classification in DIRECT_CLASSES),
            "direct_evidence_count": len(pair_evidence),
            "numeric_direct_evidence_count": numeric_direct,
            "representative_pmid": representative.get("paper_pmid", "NA"),
            "representative_doi": representative.get("paper_doi", "NA"),
            "representative_title": representative.get("article_title", "NA"),
            "representative_published_rg": representative.get("published_rg", "NA"),
            "representative_published_rg_se": representative.get("published_rg_se", "NA"),
            "representative_published_rg_p": representative.get("published_rg_p", "NA"),
            "representative_published_external_sample_size": representative.get("published_external_sample_size", representative.get("published_trait_sample_size", "NA")),
            "representative_published_sleep_sample_size": representative.get("published_sleep_sample_size", "NA"),
            "phenotype_definition_match": phenotype_match,
            "dataset_relation": dataset_relation,
            "comparison_independence": independence,
            "evidence_pmids": evidence_pmids,
            "evidence_dois": evidence_dois,
            "evidence_source_ids": evidence_sources,
            "previous_rg_values": rg_values,
            "previous_sample_sizes": sample_values,
            "papers_screened_pmids": semicolon(screened_pmids),
            "papers_inspected_in_depth_pmids": semicolon(in_depth_pmids),
            "direct_rg_query": lens_queries["DIRECT_RG"]["europe_pmc_query"],
            "ldsc_query": lens_queries["LDSC"]["europe_pmc_query"],
            "related_genetic_query": lens_queries["RELATED_GENETIC"]["europe_pmc_query"],
            "mr_observational_query": lens_queries["MR_OBSERVATIONAL"]["europe_pmc_query"],
            "query_hit_counts": semicolon([
                f"{row['evidence_lens']}={row['hit_count']}" for row in pair_queries
            ]),
            "search_as_of_date": pair_queries[0]["as_of_date"],
            "search_status": "COMPLETE_FOUR_LENS_PLUS_SOURCE_TABLE_REVIEW",
            "direct_rg_search_result": direct_search_result,
            "observational_direction_context": observation,
            "contradictory_literature": contradiction,
            "replication_feasibility": replication_feasibility,
            "resource_context": resource_context,
            "review_rationale": rationale,
            "review_status": review_status,
            "claim_limit": "Apparently novel means apparently unreported after the recorded search; it is not a first-ever claim.",
        })
    audit.sort(key=lambda row: (row["sleep_trait"], row["external_trait"]))
    if len(audit) != 153 or set(row["novelty_classification"] for row in audit) - ALLOWED_CLASSES:
        fail("literature audit invariant failed")

    replication, replication_metrics = replication_rows(significant, evidence)
    negative = read_tsv(analysis / "negative_rg.tsv")
    audit_by_pair = {(row["sleep_trait"], row["external_trait"]): row for row in audit}
    for row in negative:
        key = (row["sleep_trait"], row["external_trait"])
        reviewed = audit_by_pair[key]
        pair_evidence = evidence_by_pair.get(key, [])
        prior_negative = any(
            number(item.get("published_rg")) is not None
            and float(item["published_rg"]) < 0
            and (number(item.get("published_rg_p")) is None or float(item["published_rg_p"]) < .05)
            for item in pair_evidence
        )
        row["prior_observational_expectation"] = reviewed["observational_direction_context"]
        row["direct_negative_rg_previously_reported"] = (
            "YES_EXPLICIT_DIRECT_NEGATIVE_RG" if prior_negative
            else "DIRECT_RG_EXISTS_BUT_NEGATIVE_NOT_CONFIRMED" if reviewed["direct_rg_found"] == "True"
            else "NO_EXPLICIT_DIRECT_NEGATIVE_RG_FOUND"
        )
        row["literature_novelty_classification"] = reviewed["novelty_classification"]
        row["literature_contradiction"] = reviewed["contradictory_literature"]

    counts = Counter(row["novelty_classification"] for row in audit)
    direct_pairs = sum(counts[name] for name in DIRECT_CLASSES)
    no_direct_pairs = len(audit) - direct_pairs
    contradiction_pairs = [
        row for row in audit
        if row["contradictory_literature"] not in {"NONE_IDENTIFIED", "NONE", "No contradiction located."}
        and ("opposite" in str(row["contradictory_literature"]).lower() or "DISCORDANT" in str(row["contradictory_literature"]))
    ]
    summary: dict[str, object] = {
        "schema_version": "sleep-atlas-phase1-literature-audit.1",
        "search_as_of_date": audit[0]["search_as_of_date"],
        "significant_pairs_reviewed": len(audit),
        "four_lens_query_count": len(queries),
        "candidate_records_screened": len(candidates),
        "direct_rg_evidence_pairs": direct_pairs,
        "no_direct_rg_found_pairs": no_direct_pairs,
        "apparently_novel_pairs": counts["APPARENTLY_NOVEL"],
        "classification_counts": dict(sorted(counts.items())),
        "direct_evidence_rows": len(evidence),
        "direct_source_supplement_pairs": len({(row["sleep_trait"], row["external_trait"]) for row in source_evidence}),
        "additional_reviewed_direct_pairs": len({(row["sleep_trait"], row["external_trait"]) for row in additional_evidence}),
        "directionally_contradictory_pair_count": len(contradiction_pairs),
        "directionally_contradictory_pairs": [row["pair_id"] for row in contradiction_pairs],
        "replication_metrics": replication_metrics,
        "classification_rule": {
            "direct": "Explicit numeric or verbal genome-wide rg in a reviewed publication table, full text, or abstract.",
            "apparently_novel": "No direct or pair-specific genetic evidence in all four fixed queries plus targeted/source review, abs_rg>=0.15, FDR<=0.01, and no QC_CAUTION result.",
            "no_direct": "No explicit direct rg located; related genetic, MR, observational, uncertainty, and effect/QC-qualified apparently novel states remain distinct.",
        },
        "claim_limit": "No first-ever claims are made. Search absence is frozen to the recorded date and databases and can be overturned by missed, newly indexed, or future work.",
    }
    expected_apparent = {
        row["pair_id"] for row in audit if row["novelty_classification"] == "APPARENTLY_NOVEL"
    }
    for row in audit:
        if row["novelty_classification"] == "APPARENTLY_NOVEL":
            if float(row["abs_rg"]) < .15 or float(row["fdr"]) > .01 or row["result_confidence_class"] == "QC_CAUTION":
                fail(f"apparently novel rule violation: {row['pair_id']}")
    if len(expected_apparent) != 6:
        fail(f"expected six cautiously apparently-novel pairs, observed {len(expected_apparent)}")

    audit_fields = list(audit[0])
    replication_fields = list(replication[0])
    negative_fields = list(negative[0])
    audit_payload = tsv_text(audit, audit_fields)
    replication_payload = tsv_text(replication, replication_fields)
    negative_payload = tsv_text(negative, negative_fields)
    summary["input_hashes"] = {
        str(path.relative_to(root)): file_sha256(path) for path in [
            analysis / "phase1_master_analysis.tsv",
            analysis / "published_rg_source_evidence.tsv",
            analysis / "literature_search_queries.tsv",
            analysis / "literature_search_candidates.tsv",
            root / "config/reviewed_additional_rg_evidence.tsv",
            root / "config/literature_pair_reviews.tsv",
        ]
    }
    summary["output_content_hashes"] = {
        "literature_novelty_audit.tsv": hashlib.sha256(audit_payload.encode()).hexdigest(),
        "published_rg_replication.tsv": hashlib.sha256(replication_payload.encode()).hexdigest(),
        "negative_rg.tsv": hashlib.sha256(negative_payload.encode()).hexdigest(),
    }
    summary_json = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    summary_md = [
        "# Phase-1 literature novelty audit",
        "",
        f"All **{len(audit)}** locked primary discoveries were reviewed through four fixed Europe PMC query lenses as of **{audit[0]['search_as_of_date']}**, then cross-checked against reviewed sleep-GWAS supplements and selected full-text or main-table sources.",
        "",
        f"- Direct prior genome-wide rg evidence: **{direct_pairs} pairs**",
        f"- No direct prior rg located: **{no_direct_pairs} pairs**",
        f"- Cautiously classified apparently unreported: **{counts['APPARENTLY_NOVEL']} pairs**",
        f"- Search queries / screened candidate records: **{len(queries)} / {len(candidates)}**",
        "",
        "## Classification counts",
        "",
        "| Classification | Pairs |",
        "|---|---:|",
    ]
    summary_md.extend(f"| {name} | {counts[name]} |" for name in sorted(counts))
    summary_md.extend([
        "",
        "## Replication summary",
        "",
        f"Among {replication_metrics['exact_concept_unique_pairs']} pairs with numeric exact-concept comparators, the pair-level correlation between atlas rg and the median published rg was **{fmt(number(replication_metrics['pair_level_median_published_pearson']))}**; mean absolute difference was **{fmt(number(replication_metrics['pair_level_median_published_mean_absolute_difference']))}** and direction concordance was **{fmt(number(replication_metrics['pair_level_median_published_direction_concordance']))}**.",
        "",
        "These comparisons are positive controls, not independent replication: most reuse or overlap public GWAS inputs.",
        "",
        "## Apparently unreported candidates",
        "",
    ])
    summary_md.extend(
        f"- `{row['pair_id']}`: rg={float(row['rg']):+.4f}, FDR={float(row['fdr']):.3g}, {row['result_confidence_class']}"
        for row in sorted(
            (row for row in audit if row["novelty_classification"] == "APPARENTLY_NOVEL"),
            key=lambda item: -float(item["abs_rg"]),
        )
    )
    summary_md.extend([
        "",
        "## Interpretation limits",
        "",
        "Search hits were screened rather than accepted automatically. Co-occurrence is not direct rg evidence. `APPARENTLY_NOVEL` means apparently unreported under this dated, recorded search—not first-ever. MR, observational association, and global genetic correlation are kept distinct, and contradictions are retained in the audit.",
        "",
    ])
    summary_md_payload = "\n".join(summary_md)

    outputs = {
        analysis / "literature_novelty_audit.tsv": audit_payload,
        analysis / "published_rg_replication.tsv": replication_payload,
        analysis / "negative_rg.tsv": negative_payload,
        analysis / "literature_audit_summary.json": summary_json,
        analysis / "literature_audit_summary.md": summary_md_payload,
    }
    if args.validate_only:
        for path, payload in outputs.items():
            if not path.is_file() or path.read_text(encoding="utf-8") != payload:
                fail(f"literature audit output drifted: {path.relative_to(root)}")
        print(
            f"PHASE1_LITERATURE_AUDIT_VALID pairs={len(audit)} direct={direct_pairs} "
            f"no_direct={no_direct_pairs} apparent={counts['APPARENTLY_NOVEL']}"
        )
        return 0
    for path, payload in outputs.items():
        atomic_text(path, payload)
    print(
        f"PHASE1_LITERATURE_AUDIT_BUILT pairs={len(audit)} direct={direct_pairs} "
        f"no_direct={no_direct_pairs} apparent={counts['APPARENTLY_NOVEL']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
