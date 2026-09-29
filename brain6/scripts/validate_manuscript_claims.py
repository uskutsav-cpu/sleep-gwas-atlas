#!/usr/bin/env python3
"""Check key global-map manuscript numbers against frozen result tables."""
from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "brain6"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sci_markup(value: float, decimals: int) -> str:
    mantissa, exponent = f"{value:.{decimals}e}".split("e")
    power = int(exponent)
    sign = "−" if power < 0 else ""
    return f"{mantissa} × 10<sup>{sign}{abs(power)}</sup>"


def count_word(value: int) -> str:
    words = {
        0: "Zero", 1: "One", 2: "Two", 3: "Three", 4: "Four", 5: "Five",
        6: "Six", 7: "Seven", 8: "Eight", 9: "Nine", 10: "Ten",
        11: "Eleven", 12: "Twelve", 13: "Thirteen", 14: "Fourteen",
        15: "Fifteen", 16: "Sixteen", 17: "Seventeen", 18: "Eighteen",
        19: "Nineteen", 20: "Twenty", 21: "Twenty-one", 22: "Twenty-two",
        23: "Twenty-three", 24: "Twenty-four", 25: "Twenty-five",
        26: "Twenty-six", 27: "Twenty-seven", 28: "Twenty-eight",
        29: "Twenty-nine", 30: "Thirty", 31: "Thirty-one", 32: "Thirty-two",
        33: "Thirty-three", 34: "Thirty-four", 35: "Thirty-five",
    }
    try:
        return words[value]
    except KeyError as exc:
        raise ValueError(f"No manuscript word form configured for count {value}") from exc


def validate_citations(manuscript: str, bibliography: str) -> dict[str, int]:
    """Require cited BibTeX keys to resolve uniquely to DOI-bearing entries."""
    cited = set(re.findall(r"@([A-Za-z0-9_:-]+)", manuscript))
    entry_matches = list(re.finditer(r"(?m)^@([A-Za-z]+)\s*\{\s*([^,\s]+)\s*,", bibliography))
    keys = [match.group(2) for match in entry_matches]
    if len(keys) != len(set(keys)):
        raise ValueError("Bibliography contains duplicate citation keys")
    entries = {match.group(2): bibliography[match.start():entry_matches[i + 1].start()
                                               if i + 1 < len(entry_matches) else len(bibliography)]
               for i, match in enumerate(entry_matches)}
    missing = cited - set(entries)
    unused = set(entries) - cited
    if missing:
        raise ValueError(f"Cited keys are missing from references.bib: {sorted(missing)}")
    if unused:
        raise ValueError(f"references.bib contains uncited entries: {sorted(unused)}")
    for key, entry in entries.items():
        match = re.search(r"(?mi)^\s*doi\s*=\s*\{([^}]+)\}", entry)
        if match is None or not re.fullmatch(r"10\.\d{4,9}/\S+", match.group(1).strip()):
            raise ValueError(f"Bibliography entry {key!r} is missing a valid DOI")
    return {"cited_references": len(cited), "doi_bearing_entries": len(entries)}


def _need(text: str, expected: str, label: str, *, occurrences: int = 1) -> None:
    found = text.count(expected)
    if found < occurrences:
        raise ValueError(f"Manuscript claim {label!r} is stale or missing: expected {expected!r} at least {occurrences} time(s), found {found}")


def validate_claims(
    manuscript: str,
    global_rows: list[dict[str, str]],
    similarities: list[dict[str, str]],
    profile_leave_one_out: list[dict[str, str]],
    replication: list[dict[str, str]],
    published_context: list[dict[str, str]],
) -> dict[str, int]:
    if len(global_rows) != 72 or len({(r["sleep_trait"], r["brain_disorder"]) for r in global_rows}) != 72:
        raise ValueError("Frozen global table is not the expected 72 unique pairs")
    significant = [r for r in global_rows if r["significance_under_original_396_family"].lower() == "true"]
    by_disorder = {d: sum(r["brain_disorder"] == d for r in significant)
                   for d in ("adhd", "mdd", "scz", "bipolar", "parkinson", "alz")}

    _need(manuscript, f"The map contained {len(global_rows)} unique pairs. {count_word(len(significant))} were significant under the original atlas's 396-pair FDR procedure", "global map totals")
    _need(
        manuscript,
        "The corrected canonical LAVA v3 family completed its full receipt audit but failed the frozen untested-cell limit, so it is not promoted",
        "canonical LAVA v3 QC disposition",
    )
    _need(
        manuscript,
        "This global-map manuscript does not report local-rg, PLACO, or other downstream results as final analyses.",
        "downstream-analysis evidence boundary",
    )
    disorder_phrase = (
        f"{by_disorder['adhd']} for attention-deficit/hyperactivity disorder (ADHD), "
        f"{by_disorder['mdd']} for major depressive disorder (MDD), "
        f"{by_disorder['scz']} for schizophrenia, {by_disorder['bipolar']} for bipolar disorder, "
        f"{by_disorder['parkinson']} for Parkinson's disease, and {by_disorder['alz']} for Alzheimer's disease."
    )
    _need(manuscript, disorder_phrase, "per-disorder inherited significance counts")

    strongest = max(global_rows, key=lambda row: abs(float(row["rg"])))
    disorder_labels = {
        "adhd": "ADHD", "mdd": "MDD", "scz": "schizophrenia",
        "bipolar": "bipolar disorder", "parkinson": "Parkinson's disease", "alz": "Alzheimer's disease",
    }
    pair = f"{strongest['sleep_trait']}–{disorder_labels[strongest['brain_disorder']]}"
    estimate = (
        f"*r*<sub>g</sub> = {float(strongest['rg']):.4f}, SE = {float(strongest['se']):.4f}, "
        f"*P* = {sci_markup(float(strongest['p']), 3)}, original *q* = {sci_markup(float(strongest['original_BH_FDR_q']), 3)}"
    )
    _need(manuscript, f"largest absolute point estimate in the 72-pair subset was {pair} ({estimate})", "largest absolute global estimate")

    similarity = {(r["brain_disorder"], c): float(r[c])
                  for r in similarities for c in r if c != "brain_disorder"}
    adhd_mdd = similarity[("adhd", "mdd")]
    scz_bipolar = similarity[("scz", "bipolar")]
    profile_phrase = (
        f"ADHD and MDD had a descriptive Pearson correlation of {adhd_mdd:.3f}; "
        f"schizophrenia and bipolar disorder had a correlation of {scz_bipolar:.3f}."
    )
    _need(manuscript, profile_phrase, "abstract disorder-profile similarities")
    results_profile_phrase = (
        f"highest descriptive profile similarities were ADHD–MDD (*r* = {adhd_mdd:.3f}) "
        f"and schizophrenia–bipolar disorder (*r* = {scz_bipolar:.3f})."
    )
    _need(manuscript, results_profile_phrase, "results disorder-profile similarities")

    leave_one_out_summaries: dict[str, tuple[float, float, int]] = {}
    for claim_id in ("profile_adhd_mdd", "profile_scz_bipolar"):
        rows = [row for row in profile_leave_one_out if row["claim_id"] == claim_id]
        if len(rows) != 12 or len({row["omitted_sleep_trait"] for row in rows}) != 12:
            raise ValueError(f"Expected 12 distinct leave-one-sleep-trait-out estimates for {claim_id}")
        values = [float(row["leave_one_out_r"]) for row in rows]
        minimum, maximum = min(values), max(values)
        if any(float(row["leave_one_out_min_r"]) != float(f"{minimum:.12g}") or
               float(row["leave_one_out_max_r"]) != float(f"{maximum:.12g}") for row in rows):
            raise ValueError(f"Leave-one-out summary differs from individual rows for {claim_id}")
        leave_one_out_summaries[claim_id] = (minimum, maximum, sum(value > 0 for value in values))
    adhd_mdd_loo = leave_one_out_summaries["profile_adhd_mdd"]
    scz_bipolar_loo = leave_one_out_summaries["profile_scz_bipolar"]
    loo_phrase = (
        f"Leave-one-sleep-trait-out profile correlations ranged from {adhd_mdd_loo[0]:.3f} to "
        f"{adhd_mdd_loo[1]:.3f} for ADHD–MDD and {scz_bipolar_loo[0]:.3f} to "
        f"{scz_bipolar_loo[1]:.3f} for schizophrenia–bipolar disorder; all 12 omission estimates "
        "remained positive (Supplementary Table S19). This descriptive within-atlas sensitivity "
        "does not provide independent validation or sampling uncertainty."
    )
    _need(manuscript, loo_phrase, "leave-one-sleep-trait-out profile sensitivity")

    adhd_rep = [r for r in replication if r["sleep_trait"] == "insomnia" and r["brain_disorder"] == "adhd"
                and r["status"] == "RESULT_ARCHIVED"]
    if len(adhd_rep) != 1:
        raise ValueError(f"Expected one archived insomnia–ADHD replication row, found {len(adhd_rep)}")
    rep = adhd_rep[0]
    rep_values = (
        f"*r*<sub>g</sub> = {float(rep['replication_rg']):.4f}; SE = {float(rep['replication_SE']):.4f}; "
        f"*P* = {sci_markup(float(rep['replication_P']), 3)}"
    )
    _need(manuscript, rep_values, "abstract insomnia–ADHD replication estimate")
    rep_values_results = rep_values.replace("; SE =", ", SE =").replace("; *P* =", ", *P* =")
    _need(manuscript, rep_values_results, "results insomnia–ADHD replication estimate")

    ctx = [r for r in published_context if r["record_id"] == "austin_zimmerman_2023_longsleep_scz"]
    if len(ctx) != 1:
        raise ValueError(f"Expected one long-sleep–SCZ context row, found {len(ctx)}")
    context = ctx[0]
    context_values = (
        f"*r*<sub>g</sub> = {float(context['published_rg']):.2f}, SE = {float(context['published_se']):.2f}, "
        f"*P* = {sci_markup(float(context['published_p']), 2)}"
    )
    _need(manuscript, context_values, "published long-sleep–SCZ context")

    return {"global_pairs": len(global_rows), "original_family_significant": len(significant),
            "profile_leave_one_out_pairs": len(leave_one_out_summaries),
            "replication_rows": len(adhd_rep), "published_context_rows": len(ctx)}


def validate_local_feasibility_claims(
    manuscript: str,
    canonical_trait_rows: list[dict[str, str]],
    comparison_rows: list[dict[str, str]],
    feasibility: dict,
) -> dict[str, int]:
    """Bind manuscript eligibility claims to receipt-derived power summaries."""
    canonical_cells = sum(int(row["planned_loci"]) for row in canonical_trait_rows)
    canonical_tested = sum(int(row["tested"]) for row in canonical_trait_rows)
    canonical_not_run = sum(int(row["NOT_RUN"]) for row in canonical_trait_rows)
    low_h2 = sum(int(row["low_local_h2_NOT_RUN"]) for row in canonical_trait_rows)
    shared_ref = sum(int(row["shared_reference_minK_NOT_RUN"]) for row in canonical_trait_rows)
    other_min_k = sum(int(row["other_minK_NOT_RUN"]) for row in canonical_trait_rows)
    limit = int(feasibility["family"]["maximum_not_run"])
    if canonical_cells != int(feasibility["family"]["cells"]):
        raise ValueError("Canonical trait audit and family-feasibility denominator disagree")
    if canonical_not_run != int(feasibility["family"]["canonical_not_run"]):
        raise ValueError("Canonical trait audit and family-feasibility NOT_RUN total disagree")
    if canonical_tested + canonical_not_run != canonical_cells:
        raise ValueError("Canonical TESTED and NOT_RUN statuses do not account for every cell")
    if low_h2 + shared_ref + other_min_k != canonical_not_run:
        raise ValueError("Canonical NOT_RUN causes do not account for every NOT_RUN cell")
    _need(
        manuscript,
        f"{canonical_tested:,} of {canonical_cells:,} cells as `TESTED` and "
        f"{canonical_not_run:,} as `NOT_RUN` ({canonical_not_run / canonical_cells:.2%}), "
        f"exceeding the frozen 5% ceiling of {limit:,}",
        "canonical local-analysis family status counts",
    )
    _need(
        manuscript,
        f"{low_h2:,} low-local-h² failures, {shared_ref:,} insufficient-shared-reference-variant "
        f"failures, and {other_min_k:,} other minimum-*K* failures",
        "canonical local-analysis failure causes",
    )

    by_trait = {row["trait_id"]: row for row in canonical_trait_rows}
    screen = {row["metric"]: row for row in comparison_rows}
    longsleep_tested = int(by_trait["longsleep"]["tested"])
    longsleep_not_run = int(by_trait["longsleep"]["NOT_RUN"])
    longsleep_low_h2 = int(by_trait["longsleep"]["low_local_h2_NOT_RUN"])
    candidate_tested = int(float(screen["processable_tested_loci"]["continuous_duration_candidate"]))
    candidate_low_h2 = int(float(screen["low_local_h2_not_run"]["continuous_duration_candidate"]))
    tested_gain = int(float(screen["processable_tested_loci"]["difference"]))
    processability_gain = float(screen["processable_percent"]["difference"])
    low_h2_reduction = float(screen["low_local_h2_relative_reduction_percent"]["continuous_duration_candidate"])
    _need(
        manuscript,
        f"{candidate_tested:,} of 2,495 loci ({candidate_tested / 2495:.2%}), compared with "
        f"{longsleep_tested:,} ({longsleep_tested / 2495:.2%}) for the canonical long-sleep tail: "
        f"an increase of {tested_gain:,} loci or {processability_gain:.2f} percentage points",
        "continuous-duration trait-level processability comparison",
    )
    _need(
        manuscript,
        f"Low-local-h² `NOT_RUN` cells decreased from {longsleep_low_h2:,} to "
        f"{candidate_low_h2:,} ({low_h2_reduction:.2f}%)",
        "continuous-duration low-local-h2 comparison",
    )
    observed_candidate_not_run = int(feasibility["family"]["candidate_family_not_run"])
    both_sleep_lower_bound = int(feasibility["best_case_lower_bounds"]["both_sleep_traits_zero_not_run"]["retained_five_disorder_traits_not_run"])
    _need(
        manuscript,
        f"The substituted two-sleep/five-disorder family still had {observed_candidate_not_run:,} `NOT_RUN` cells, "
        f"above the same {limit:,}-cell ceiling",
        "substituted family gate",
    )
    _need(
        manuscript,
        f"the five disorder inputs alone leave a lower bound of {both_sleep_lower_bound:,} `NOT_RUN` cells",
        "best-case family feasibility lower bound",
    )
    return {"canonical_cells": canonical_cells, "canonical_tested": canonical_tested,
            "canonical_not_run": canonical_not_run, "candidate_family_not_run": observed_candidate_not_run}


def validate_pgc_mdd2025_sensitivity_claims(manuscript: str, rows: list[dict[str, str]]) -> dict[str, str]:
    if len(rows) != 1:
        raise ValueError(f"Expected one PGC MDD2025 sensitivity record, found {len(rows)}")
    row = rows[0]
    expected = {
        "rg": "0.4771", "se": "0.0228", "p": "4.6535e-97",
        "cross_trait_intercept": "0.0087", "cross_trait_intercept_se": "0.0069",
        "hm3_valid_allele_snps": "1061081",
        "status": "SENSITIVITY_ONLY_NOT_INDEPENDENT_REPLICATION",
    }
    if any(row.get(key) != value for key, value in expected.items()):
        raise ValueError("PGC MDD2025 sensitivity row has changed or lost its non-independent label")
    required_claims = (
        "*r*<sub>g</sub> = 0.4771 (SE = 0.0228, *P* = 4.65 × 10<sup>−97</sup>",
        "cross-trait intercept was 0.0087 (SE = 0.0069)",
        "nine PGC cohort entries",
        "not an independent replication",
        "Supplementary Table S34",
    )
    missing = [claim for claim in required_claims if claim not in manuscript]
    if missing:
        raise ValueError(f"Manuscript omits PGC MDD2025 sensitivity details: {missing}")
    return {"status": row["status"], "rg": row["rg"], "se": row["se"], "p": row["p"]}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_grover_sleep_pd_context(manuscript: str, supplement: str) -> dict[str, str]:
    table_path = OUT / "results/novelty/grover_2022_sleep_pd_ldsc_context.tsv"
    provenance_path = table_path.with_suffix(".provenance.json")
    builder_path = OUT / "scripts/build_grover_2022_sleep_pd_ldsc_context.py"
    source_path = OUT / "manifests/gwas_master.tsv"
    rows = read_tsv(table_path)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    source = provenance.get("source", {})
    supplementary = source.get("linked_supplementary_file", {})
    if (provenance.get("status") != "PASS_PRIOR_NEGATIVE_LDSC_CONTEXT" or
            provenance.get("builder_sha256") != sha256_file(builder_path) or
            provenance.get("input_hashes", {}).get(str(source_path.relative_to(OUT.parent))) != sha256_file(source_path) or
            provenance.get("output", {}).get("sha256") != sha256_file(table_path) or
            provenance.get("source_file_hashes", {}).get(
                "brain6/qc/replication_sources/grover_2022_article_fullTextXML.xml") !=
                sha256_file(OUT / "qc/replication_sources/grover_2022_article_fullTextXML.xml") or
            len(rows) != 1):
        raise ValueError("Grover long-sleep–Parkinson context failed its source/code/output hash checks")
    row = rows[0]
    if (row.get("sleep_trait") != "longsleep" or row.get("brain_disorder") != "parkinson" or
            row.get("published_rg") != "NA" or row.get("numeric_pair_estimate_reviewed") != "False" or
            row.get("independent_replication") != "False" or
            source.get("numeric_rg_se_p_in_article_text") is not False or
            source.get("linked_supplementary_workbook_reviewed") is not False or
            supplementary.get("filename") != "Table_1.xlsx" or
            supplementary.get("size_bytes") != 722516 or
            supplementary.get("binary_retrieved") is not False or
            supplementary.get("numeric_pair_estimate_reviewed") is not False):
        raise ValueError("Grover long-sleep–Parkinson context overstates an unreported pair estimate")
    if ("[@grover2022sleeppd]" not in manuscript or "Table S35" not in manuscript or
            "Table S35" not in supplement or "no numeric estimate or exact zero is inferred" not in supplement or
            "`MR_ONLY`" not in manuscript):
        raise ValueError("Manuscript omits the qualified long-sleep–Parkinson prior-context limits")
    return {"status": provenance["status"], "pair": "longsleep-parkinson", "numeric_rg_available": "False"}


def validate_sleepchart_2026_longsleep_context(manuscript: str, supplement: str) -> dict[str, object]:
    table_path = OUT / "results/novelty/sleepchart_2026_long_sleep_rg_context.tsv"
    provenance_path = table_path.with_suffix(".provenance.json")
    builder_path = OUT / "scripts/build_sleepchart_2026_long_sleep_rg_context.py"
    rows = read_tsv(table_path)
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    expected = {"adhd": "0.28", "mdd": "0.29", "scz": "0.28", "bipolar": "0.21"}
    if (provenance.get("status") != "PASS_PRIOR_CONTEXT_ROUNDED_ESTIMATES_NOT_INDEPENDENT_REPLICATION" or
            provenance.get("builder_sha256") != sha256_file(builder_path) or
            provenance.get("output", {}).get("sha256") != sha256_file(table_path) or
            len(rows) != 4 or
            {row.get("brain_disorder") for row in rows} != set(expected)):
        raise ValueError("SleepChart 2026 long-sleep context failed source/code/output checks")
    for relative, expected_hash in provenance.get("input_hashes", {}).items():
        source_path = OUT.parent / relative
        if not source_path.is_file() or sha256_file(source_path) != expected_hash:
            raise ValueError(f"SleepChart input changed: {relative}")
    for row in rows:
        if (row.get("sleepchart_rg_rounded") != expected[row["brain_disorder"]] or
                row.get("direction_concordant") != "True" or
                row.get("independent_replication") != "False" or
                row.get("power_replacement_eligible") != "False" or
                row.get("sleepchart_sleep_n_cases") != "25049" or
                row.get("sleepchart_sleep_n_controls") != "300420"):
            raise ValueError("SleepChart estimate, source sample, direction, or interpretation changed")
    if (provenance.get("scope", {}).get("novelty_or_fdr_changed") is not False or
            provenance.get("source", {}).get("sleep_gwas", {}).get("source") != "UK Biobank field 1160" or
            provenance.get("source", {}).get("sleep_gwas", {}).get("summary_statistics_downloaded") is not False or
            "[@multiconsortium2026sleepchart]" not in manuscript or "Table S40" not in manuscript or
            "not independent replication" not in manuscript or
            "Table S40" not in supplement or "does not improve expected power" not in supplement):
        raise ValueError("SleepChart context is missing or overstates independence/power relevance")
    return {"status": provenance["status"], "rows": len(rows), "directions_concordant": 4,
            "independent_replication": False, "power_replacement_eligible": False}


def main() -> None:
    manuscript = (OUT / "paper/global_map_manuscript.md").read_text(encoding="utf-8")
    citations = validate_citations(manuscript, (OUT / "paper/references.bib").read_text(encoding="utf-8"))
    result = validate_claims(
        manuscript,
        read_tsv(OUT / "results/global/brain6_72_locked.tsv"),
        read_tsv(OUT / "results/global/disorder_profile_similarity.tsv"),
        read_tsv(OUT / "results/global/disorder_profile_leave_one_sleep_trait_out.tsv"),
        read_tsv(OUT / "results/replication/replication_master.tsv"),
        read_tsv(OUT / "results/supplement/table_S18_published_rg_context.tsv"),
    )
    result["pgc_mdd2025_sensitivity"] = validate_pgc_mdd2025_sensitivity_claims(
        manuscript, read_tsv(OUT / "results/replication/pgc_mdd2025_insomnia_rg_sensitivity_v1.tsv"))
    result["grover_2022_sleep_pd_context"] = validate_grover_sleep_pd_context(
        manuscript, (OUT / "paper/supplement.md").read_text(encoding="utf-8"))
    result["sleepchart_2026_longsleep_context"] = validate_sleepchart_2026_longsleep_context(
        manuscript, (OUT / "paper/supplement.md").read_text(encoding="utf-8"))
    feasibility = validate_local_feasibility_claims(
        manuscript,
        read_tsv(OUT / "results/power_optimized_sensitivity_v1/canonical_family_trait_power_audit_v1.tsv"),
        read_tsv(OUT / "results/power_optimized_sensitivity_v1/completed_trait_only_screen_comparison_v1.tsv"),
        json.loads((OUT / "results/power_optimized_sensitivity_v1/sensitivity_family_feasibility_v1.json").read_text(encoding="utf-8")),
    )
    result["local_feasibility"] = feasibility
    print(json.dumps({"validation": "PASS_MANUSCRIPT_CLAIMS", **result, **citations}, indent=2))


if __name__ == "__main__":
    main()
