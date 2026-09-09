#!/usr/bin/env python3
"""Write the evidence-linked Phase-1 deep scientific analysis report."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        fail(f"missing non-empty input: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_text(path: Path, payload: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(path)


def fdr_text(value: str) -> str:
    return f"{float(value):.3g}"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", default=".")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    analysis = root / "results/analysis"

    master = read_tsv(analysis / "phase1_master_analysis.tsv")
    primary = [row for row in master if row["primary_or_sensitivity"] == "PRIMARY_PHASE1"]
    sensitivity = [row for row in master if row["primary_or_sensitivity"] == "QC_FAILED_SENSITIVITY"]
    discoveries = [row for row in primary if row["locked_primary_significant"] == "True"]
    audit = read_tsv(analysis / "literature_novelty_audit.tsv")
    audit_by_pair = {(row["sleep_trait"], row["external_trait"]): row for row in audit}
    hubs = read_tsv(analysis / "sleep_trait_hubs.tsv")
    external_hubs = read_tsv(analysis / "external_trait_hubs.tsv")
    domains = read_tsv(analysis / "domain_architecture.tsv")
    similarity = read_tsv(analysis / "sleep_profile_similarity.tsv")
    measurement = read_tsv(analysis / "measurement_mode_comparison.tsv")
    loo = read_tsv(analysis / "leave_one_out_sensitivity.tsv")
    testing = read_tsv(analysis / "multiple_testing_sensitivity.tsv")
    confidence = read_tsv(analysis / "result_confidence.tsv")
    priorities = read_tsv(analysis / "novel_connection_priorities.tsv")
    local = read_tsv(analysis / "local_followup_priorities.tsv")
    replication = json.loads((analysis / "literature_audit_summary.json").read_text(encoding="utf-8"))["replication_metrics"]
    verification = json.loads((analysis / "input_verification.json").read_text(encoding="utf-8"))

    if (len(master), len(primary), len(sensitivity), len(discoveries), len(audit)) != (396, 372, 24, 153, 153):
        fail("report input invariant failed")

    positive = sum(float(row["rg"]) > 0 for row in discoveries)
    negative = sum(float(row["rg"]) < 0 for row in discoveries)
    strongest_positive = max(discoveries, key=lambda row: float(row["rg"]))
    strongest_negative = min(discoveries, key=lambda row: float(row["rg"]))
    apparent = [row for row in audit if row["novelty_classification"] == "APPARENTLY_NOVEL"]
    strongest_apparent = max(apparent, key=lambda row: float(row["abs_rg"]))
    direct_count = sum(row["direct_rg_found"] == "True" for row in audit)
    biggest_degree = max(int(row["fdr_significant_connections"]) for row in hubs)
    biggest_degree_hubs = [row["sleep_trait"] for row in hubs if int(row["fdr_significant_connections"]) == biggest_degree]
    weighted_hub = max(hubs, key=lambda row: float(row["sum_abs_rg_significant"]))
    specialized = min(hubs, key=lambda row: (int(row["fdr_significant_connections"]), int(row["connected_domain_count"])))
    external_hub = max(external_hubs, key=lambda row: int(row["significant_sleep_connections"]))
    enriched = max(domains, key=lambda row: float(row["significant_proportion"]))
    max_count_domain = max(domains, key=lambda row: int(row["fdr_significant_pairs"]))
    closest = max(similarity, key=lambda row: float(row["pearson_rg_profile"]))
    farthest = min(similarity, key=lambda row: float(row["pearson_rg_profile"]))
    confidence_by_key = {(row["sleep_trait"], row["external_trait"]): row for row in confidence}
    confidence_counts = Counter(
        confidence_by_key[(row["sleep_trait"], row["external_trait"])]["confidence_class"]
        for row in discoveries
    )
    primary_testing = [row for row in testing if row["primary_or_sensitivity"] == "PRIMARY_PHASE1"]
    primary_only_count = sum(row["primary_372_fdr_significant_sensitivity"] == "True" for row in primary_testing)
    bonferroni_count = sum(row["bonferroni_396_significant"] == "True" for row in primary_testing)
    fdr01_count = sum(row["locked_fdr_le_0_01"] == "True" for row in primary_testing)
    fdr001_count = sum(row["locked_fdr_le_0_001"] == "True" for row in primary_testing)
    no_insomnia = next(row for row in loo if row["omission_type"] == "SLEEP_TRAIT" and row["omitted_unit"] == "insomnia")
    largest_subjective_gaps = sorted(
        (row for row in measurement if row["subjective_materially_larger"] == "True"),
        key=lambda row: float(row["objective_minus_subjective_mean_abs_rg"]),
    )

    top_underexplored = sorted(apparent, key=lambda row: -float(row["abs_rg"]))
    additional_underexplored = sorted(
        (
            row for row in audit
            if row["novelty_classification"] not in {
                "APPARENTLY_NOVEL", "DIRECT_RG_PREVIOUSLY_REPORTED",
                "DIRECT_RG_REPLICATION_DIFFERENT_DATASET",
            }
            and row["result_confidence_class"] != "QC_CAUTION"
        ),
        key=lambda row: -float(row["abs_rg"]),
    )
    top_ten = (top_underexplored + additional_underexplored)[:10]

    top_ten_lines = [
        f"{index}. `{row['pair_id']}` — rg={float(row['rg']):+.4f}, FDR={fdr_text(row['fdr'])}, {row['novelty_classification']}."
        for index, row in enumerate(top_ten, 1)
    ]
    apparent_lines = [
        f"- `{row['pair_id']}`: rg={float(row['rg']):+.4f}, FDR={fdr_text(row['fdr'])}, {row['result_confidence_class']}; {row['replication_feasibility']}."
        for row in sorted(apparent, key=lambda item: -float(item["abs_rg"]))
    ]
    negative_lines = [
        f"- `{row['sleep_trait']}__{row['external_trait']}`: rg={float(row['rg']):+.4f}, FDR={fdr_text(row['fdr'])}, {audit_by_pair[(row['sleep_trait'], row['external_trait'])]['novelty_classification']}."
        for row in sorted((row for row in discoveries if float(row["rg"]) < 0), key=lambda item: float(item["rg"]))[:10]
    ]
    tier_lines = [
        f"- Rank {row['priority_rank']}: `{row['pair_id']}` ({row['priority_tier']}), rg={float(row['rg']):+.4f}, FDR={fdr_text(row['fdr'])}."
        for row in priorities[:15]
    ]
    local_category_counts = {
        letter: sum(letter in row["priority_categories"].split(";") for row in local)
        for letter in "ABCDEF"
    }
    subjective_gap_text = ", ".join(
        f"{row['external_trait']} ({float(row['objective_minus_subjective_mean_abs_rg']):.3f})"
        for row in largest_subjective_gaps[:5]
    )

    lines = [
        "# Phase-1 deep analysis of the sleep–disease genetic-correlation atlas",
        "",
        "## Executive scientific summary",
        "",
        f"- **Primary discoveries:** {len(discoveries)} of {len(primary)} primary pairs under the immutable 396-family BH-FDR definition ({positive} positive, {negative} negative).",
        f"- **Biggest sleep hub:** `{weighted_hub['sleep_trait']}` by total significant |rg| ({float(weighted_hub['sum_abs_rg_significant']):.4f}); `{', '.join(biggest_degree_hubs)}` tie by degree at {biggest_degree} connections.",
        f"- **Most specialized sleep phenotype:** `{specialized['sleep_trait']}`, with {specialized['fdr_significant_connections']} significant connection in {specialized['connected_domain_count']} domain.",
        f"- **Largest external hub:** `{external_hub['external_trait']}`, connected to {external_hub['significant_sleep_connections']} sleep traits.",
        f"- **Most interconnected domain by proportion:** `{enriched['external_domain']}` ({enriched['fdr_significant_pairs']}/{enriched['tested_pairs']}, {100*float(enriched['significant_proportion']):.1f}%); permutation P={float(enriched['permutation_enrichment_p']):.4f}, BH q={float(enriched['permutation_enrichment_bh_q']):.4f}. `{max_count_domain['external_domain']}` has the largest raw count ({max_count_domain['fdr_significant_pairs']}).",
        f"- **Strongest positive:** `{strongest_positive['sleep_trait']}__{strongest_positive['external_trait']}`, rg={float(strongest_positive['rg']):+.4f}, FDR={fdr_text(strongest_positive['fdr'])}.",
        f"- **Strongest negative:** `{strongest_negative['sleep_trait']}__{strongest_negative['external_trait']}`, rg={float(strongest_negative['rg']):+.4f}, FDR={fdr_text(strongest_negative['fdr'])}.",
        f"- **Strongest apparently unreported pair:** `{strongest_apparent['pair_id']}`, rg={float(strongest_apparent['rg']):+.4f}, FDR={fdr_text(strongest_apparent['fdr'])}.",
        f"- **Prior direct rg evidence:** {direct_count}/153; **no direct prior rg located:** {len(audit)-direct_count}/153; **conservative apparently-unreported replication candidates:** {len(apparent)}.",
        "",
        "### Top 10 connections to investigate next",
        "",
        *top_ten_lines,
        "",
        "Selection rule: all conservative `APPARENTLY_NOVEL` pairs first by descending |rg|, followed by clean no-direct-rg pairs by descending |rg|. This is a transparent lexicographic rule, not a hidden score.",
        "",
        "## 1. Scope, immutable inputs, and result reconstruction",
        "",
        "**OBSERVED RESULT.** The locked matrix contains 396 unique sleep–external pairs: 372 primary rows and 24 explicitly segregated T2D/melanoma sensitivity rows. All 153 headline discoveries reproduce directly from the locked 396-family FDR column. Recomputing BH only across the 372 primary rows gives 155, which is reported solely as sensitivity and never substituted for the locked definition. The requested legacy `results/tables/phase1_completion.json` is absent; the canonical completion-equivalent provenance is `results/atlas/core.provenance.json`. A stale legacy Phase-1 Markdown summary was not used.",
        "",
        f"**OBSERVED RESULT.** Input verification status is `{verification['verification_status']}`. The master table retains every requested quantitative field, source/QC label, ancestry, sample size, overlap, h² statistic, intercept, and immutable primary/sensitivity status.",
        "",
        "**INTERPRETATION.** The atlas can be analyzed defensibly without changing the underlying Phase-1 results, provided the 153/155 metadata discrepancy and absent legacy completion filename remain explicit.",
        "",
        "## 2. Dominant architecture",
        "",
        f"**OBSERVED RESULT.** Adverse or dysregulated self-reported sleep traits form the dominant positive cardiometabolic–psychiatric–frailty axis. The closest whole-profile pair is `{closest['sleep_trait_1']}__{closest['sleep_trait_2']}` (Pearson r={float(closest['pearson_rg_profile']):.3f}); the most opposed profiles are `{farthest['sleep_trait_1']}__{farthest['sleep_trait_2']}` (r={float(farthest['pearson_rg_profile']):.3f}). PC1 explains 76.6% of descriptive profile variance and PC2 8.5%.",
        "",
        "**INTERPRETATION.** Insomnia, short sleep, long sleep, sleepiness, napping, snoring, and sleep apnea largely share a broad adverse-health rg profile. Continuous sleep duration and objective efficiency often point in the opposite direction, so duration extremes and sleep quality should not be collapsed into one construct.",
        "",
        "**HYPOTHESIS.** Genomic SEM should test a broad adverse-sleep factor plus measurement-specific residual factors. This descriptive PCA/clustering result is not itself a latent biological model.",
        "",
        "## 3. Sleep and external hubs",
        "",
        f"**OBSERVED RESULT.** Insomnia, long sleep, and napping each connect to 21 external traits. Insomnia is the largest weighted hub (sum significant |rg|={float(weighted_hub['sum_abs_rg_significant']):.4f}); sleep apnea is most intense on average among its significant edges. Sleep timing is most specialized, with only the rheumatoid-arthritis connection passing locked FDR. Major depression is the external hub with 10 sleep connections; ADHD, BMI, parental lifespan, and frailty each have nine.",
        "",
        f"**OBSERVED RESULT.** Removing insomnia leaves {no_insomnia['significant_connections']} significant edges, with `{no_insomnia['top_remaining_sleep_hub']}` as the leading degree hub(s) and `{no_insomnia['top_remaining_domain']}` as the leading domain.",
        "",
        "**INTERPRETATION.** Insomnia strengthens the story but does not create it. The architecture remains distributed across long sleep, napping, short sleep, apnea, metabolic traits, psychiatric traits, and aging endpoints.",
        "",
        "## 4. Domain architecture and cross-domain structure",
        "",
        f"**OBSERVED RESULT.** Psychiatric traits have the highest connected proportion ({enriched['fdr_significant_pairs']}/{enriched['tested_pairs']}); aging has the largest count ({max_count_domain['fdr_significant_pairs']}). The psychiatric permutation enrichment is nominal (P={float(enriched['permutation_enrichment_p']):.4f}) but narrowly misses domain-level BH significance (q={float(enriched['permutation_enrichment_bh_q']):.4f}). Cancer and neuro show nominal depletion signals, also not robust after domain-level correction.",
        "",
        "**OBSERVED RESULT.** Data-driven external-trait clustering selects two broad cross-domain profile clusters (mean silhouette 0.585), rather than recovering the seven labels cleanly. Network communities are therefore cross-domain by construction of the observed rg profiles.",
        "",
        "**INTERPRETATION.** Sleep-related genetic sharing cuts across clinical taxonomies. Psychiatric, metabolic, cardiovascular, immune, cancer, and aging labels remain useful summaries, but the rg profile does not respect them as isolated modules.",
        "",
        "**HYPOTHESIS.** PLACO and multivariate fine-mapping should test shared loci at cross-domain bridge edges, while Genomic SEM should test whether a broad adverse-health factor explains their covariance.",
        "",
        "## 5. Objective versus subjective sleep",
        "",
        f"**OBSERVED RESULT.** Subjective mean |rg| is materially larger than objective mean |rg| for nine external traits. The five largest objective-minus-subjective mean-|rg| differences are {subjective_gap_text}. Direction discordance occurs for multiple psychiatric, metabolic, immune, and aging endpoints. Objective sleep efficiency often reverses the adverse self-report pattern; actigraphy sleep duration has only four significant links and sleep timing only one.",
        "",
        "**INTERPRETATION.** Self-report, disease-like apnea, actigraphy duration, efficiency, and timing are not interchangeable exposures. Differences may reflect biology, measurement error, sample size, ascertainment, or all four.",
        "",
        "**HYPOTHESIS.** Measurement-stratified Genomic SEM and local rg should test shared versus method-specific components; the current group means are descriptive and not independence-based tests.",
        "",
        "## 6. Strong negative relationships",
        "",
        *negative_lines,
        "",
        "**INTERPRETATION.** The aging pattern is coherent: adverse sleep tends to correlate negatively with lifespan, telomere, HDL, and grip-strength phenotypes. Continuous duration and objective efficiency can instead correlate negatively with disease/frailty. Immune and cancer negatives are more heterogeneous and sensitive to phenotype/source definitions.",
        "",
        "**HYPOTHESIS.** Signed local analyses should distinguish same-locus antagonistic sharing from mixtures of positive and negative local effects. A negative global rg alone does not establish antagonistic pleiotropy.",
        "",
        "## 7. Published-rg replication",
        "",
        f"**OBSERVED RESULT.** Explicit prior direct-rg evidence was found for {direct_count} pairs. Among {replication['exact_concept_unique_pairs']} pairs with numeric exact-concept comparators, atlas rg versus median published rg has Pearson r={replication['pair_level_median_published_pearson']:.3f}, mean absolute difference={replication['pair_level_median_published_mean_absolute_difference']:.3f}, and direction concordance={100*replication['pair_level_median_published_direction_concordance']:.1f}%. Comparable-SE results are mutually consistent within combined 95% uncertainty in {100*replication['proportion_consistent_within_combined_95pct_uncertainty']:.1f}% of rows.",
        "",
        "**OBSERVED RESULT.** Important retained conflicts include sleep apnea–atrial fibrillation (atlas +0.1304 versus published -0.3056) and sleep apnea–triglycerides (atlas +0.2144 versus preprint -0.0118, nonsignificant). The latter source is a 2025 preprint and is labeled accordingly.",
        "",
        "**INTERPRETATION.** The high agreement supports pipeline validity, not independent biological replication: most sources reuse or overlap public GWAS inputs. Dataset relations and phenotype matches are printed row by row.",
        "",
        "## 8. Literature novelty audit",
        "",
        f"**OBSERVED RESULT.** Four fixed search lenses were run for every discovery (612 queries; 1,195 candidate records), then supplemented with checksum-pinned published tables and targeted full-text review. Classification counts are: {', '.join(f'{key}={value}' for key, value in sorted(Counter(row['novelty_classification'] for row in audit).items()))}.",
        "",
        *apparent_lines,
        "",
        "**INTERPRETATION.** These six are “apparently unreported” as of 2026-08-28, not first-ever claims. The rule requires no located direct or pair-specific genetic evidence, |rg|≥0.15, FDR≤0.01, and no QC_CAUTION result.",
        "",
        "**HYPOTHESIS.** Their main scientific value is as falsifiable replication targets, especially aging endpoints whose independent exact-match GWAS options are currently limited.",
        "",
        "## 9. Priority candidates",
        "",
        *tier_lines,
        "",
        "The full 25-row table preserves five Tier A, ten Tier B, five Tier C, and five Tier D controls. Tier A permits an alternative-dataset replication design but explicitly records when an exact independent parental-lifespan or healthspan GWAS is limited; this caveat should travel with every claim.",
        "",
        "## 10. Local and multivariate follow-up",
        "",
        f"**OBSERVED RESULT.** The local-analysis set contains {len(local)} unique pairs: category A apparently novel={local_category_counts['A']}, B known controls={local_category_counts['B']}, C negatives={local_category_counts['C']}, D measurement discordance={local_category_counts['D']}, E globally null cancellation candidates={local_category_counts['E']}, and F bridge edges={local_category_counts['F']}. Categories can overlap.",
        "",
        "**HYPOTHESIS.**",
        "",
        "1. **LAVA / HDL-L:** test signed local covariance for all Tier A pairs, unusual negatives, and the 16 globally weak/null cancellation candidates.",
        "2. **MiXeR:** test polygenic overlap despite near-zero global rg, particularly where related sleep phenotypes show strong effects.",
        "3. **PLACO:** test pleiotropic loci at broad cross-domain bridge edges and known controls before applying it to novel candidates.",
        "4. **Colocalization / fine-mapping:** restrict mechanistic claims to loci with compatible causal configurations and molecular-QTL support.",
        "5. **Genomic SEM:** compare an adverse-sleep factor, a duration-extremes factor, and objective-measure residuals; confirm that factor interpretation is stable to trait omission.",
        "",
        "No local cancellation, antagonistic pleiotropy, mediation, or causality is claimed in this report.",
        "",
        "## 11. Findings not to emphasize",
        "",
        f"**OBSERVED RESULT.** Significant-result QC classes are {', '.join(f'{key}={value}' for key, value in sorted(confidence_counts.items()))}. T2D and melanoma comprise all {len(sensitivity)} sensitivity rows and are excluded from every headline count, literature classification, priority table, and figure-12 status comparison.",
        "",
        "- Do not headline T2D or melanoma; they failed upstream QC gates and remain sensitivity only.",
        "- Do not treat the FinnGen CAD proxy as identical to narrowly defined CAD, FinnGen asthma/MS substitutions as the preferred discovery endpoints, or the NFE telomere substitution as the full mixed-ancestry release.",
        "- Do not promote sleep-efficiency–ovarian cancer as a clean discovery: its external h² Z is low and the pair is QC_CAUTION.",
        "- Do not call OSA–atrial fibrillation or OSA–triglycerides replicated; their reviewed directions conflict or are unstable.",
        "- Do not interpret network betweenness, clustering, or PCA axes as causal biology.",
        "- Do not infer causal direction from genetic correlation, MR co-occurrence, or observational consistency.",
        "",
        "## 12. Multiple-testing and robustness sensitivity",
        "",
        f"**OBSERVED RESULT.** The immutable primary result remains 153. Sensitivity counts among the 372 primary rows are: primary-only BH={primary_only_count}, Bonferroni across 396={bonferroni_count}, locked FDR≤0.01={fdr01_count}, and locked FDR≤0.001={fdr001_count}. The primary-only BH family adds two calls but does not replace the locked family.",
        "",
        f"**OBSERVED RESULT.** Omitting insomnia removes its 21 discoveries but leaves {no_insomnia['significant_connections']} and preserves aging as the leading domain. Omitting long sleep or napping similarly leaves 132. Domain omission changes the leading label when that domain is removed, but a multi-domain architecture remains.",
        "",
        "**INTERPRETATION.** The scientific story is not a single-trait artifact, although insomnia, long sleep, and napping materially determine edge count. Source-substitution summaries are descriptive because ancestry, endpoint breadth, and ascertainment differ.",
        "",
        "## 13. Figures and output map",
        "",
        "All figures are generated by `scripts/103_plot_phase1_atlas.py`; exact bytes and SHA-256 hashes are recorded in the figure manifest.",
        "",
        "1. `results/figures/phase1_deep_analysis/phase1_01_full_rg_heatmap.png`",
        "2. `results/figures/phase1_deep_analysis/phase1_02_primary_fdr_heatmap.png`",
        "3. `results/figures/phase1_deep_analysis/phase1_03_sleep_trait_hubs.png`",
        "4. `results/figures/phase1_deep_analysis/phase1_04_external_trait_hubs.png`",
        "5. `results/figures/phase1_deep_analysis/phase1_05_domain_connectivity_matrix.png`",
        "6. `results/figures/phase1_deep_analysis/phase1_06_sleep_profile_similarity.png`",
        "7. `results/figures/phase1_deep_analysis/phase1_07_sleep_profile_dendrogram.png`",
        "8. `results/figures/phase1_deep_analysis/phase1_08_sleep_profile_pca.png`",
        "9. `results/figures/phase1_deep_analysis/phase1_09_bipartite_network.png`",
        "10. `results/figures/phase1_deep_analysis/phase1_10_negative_rg_network.png`",
        "11. `results/figures/phase1_deep_analysis/phase1_11_objective_subjective_comparison.png`",
        "12. `results/figures/phase1_deep_analysis/phase1_12_known_vs_novel_effect_size.png`",
        "",
        "Core tables: `phase1_master_analysis.tsv`, `literature_novelty_audit.tsv`, `published_rg_replication.tsv`, `novel_connection_priorities.tsv`, `local_followup_priorities.tsv`, and `top_findings.tsv`. The final table contains five separate 25-row views—effect magnitude, statistical evidence, novelty, confidence, and mechanistic follow-up—with its lexicographic ranking rule printed on every row.",
        "",
        "## 14. Literature sources and claim limits",
        "",
        "Representative primary sources include the published insomnia–frailty analysis ([PMCID PMC10903740](https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/)), chronotype–breast-cancer analysis ([PMID 37380357](https://pubmed.ncbi.nlm.nih.gov/37380357/)), insomnia–telomere analysis ([PMID 36182724](https://pubmed.ncbi.nlm.nih.gov/36182724/)), OSA cardiometabolic study ([PMID 33243845](https://pubmed.ncbi.nlm.nih.gov/33243845/)), OSA–BMI analysis ([PMCID PMC11166156](https://pmc.ncbi.nlm.nih.gov/articles/PMC11166156/)), and the OSA multi-phenotype preprint ([PMCID PMC12642705](https://pmc.ncbi.nlm.nih.gov/articles/PMC12642705/)). The pair-complete evidence, exact phenotype/dataset relation, source URL, PMID/DOI, numeric published rg, and contradictions are in the audit and replication tables.",
        "",
        "Search absence can be overturned by missed indexing, inaccessible supplements, newly published work, or a different phenotype definition. “Apparently unreported” is therefore dated and conditional. Genetic correlation describes genome-wide sharing, not causation or a single mechanism. Every mechanistic proposal above is a hypothesis for downstream analysis.",
        "",
    ]
    payload = "\n".join(lines)
    output = analysis / "PHASE1_DEEP_ANALYSIS.md"
    if args.validate_only:
        if not output.is_file() or output.read_text(encoding="utf-8") != payload:
            fail("Phase-1 deep analysis report drifted")
        print(f"PHASE1_REPORT_VALID discoveries={len(discoveries)} direct={direct_count} apparent={len(apparent)}")
        return 0
    atomic_text(output, payload)
    print(f"PHASE1_REPORT_BUILT discoveries={len(discoveries)} direct={direct_count} apparent={len(apparent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
