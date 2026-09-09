# Phase-1 deep analysis of the sleep–disease genetic-correlation atlas

## Executive scientific summary

- **Primary discoveries:** 153 of 372 primary pairs under the immutable 396-family BH-FDR definition (98 positive, 55 negative).
- **Biggest sleep hub:** `insomnia` by total significant |rg| (4.2123); `insomnia, longsleep, napping` tie by degree at 21 connections.
- **Most specialized sleep phenotype:** `sleep_timing`, with 1 significant connection in 1 domain.
- **Largest external hub:** `mdd`, connected to 10 sleep traits.
- **Most interconnected domain by proportion:** `psychiatric` (34/48, 70.8%); permutation P=0.0080, BH q=0.0559. `aging` has the largest raw count (39).
- **Strongest positive:** `insomnia__frailty`, rg=+0.6405, FDR=1.71e-160.
- **Strongest negative:** `shortsleep__parental_lifespan`, rg=-0.2976, FDR=2.02e-16.
- **Strongest apparently unreported pair:** `sleep_apnea__healthspan`, rg=+0.4005, FDR=3.08e-20.
- **Prior direct rg evidence:** 97/153; **no direct prior rg located:** 56/153; **conservative apparently-unreported replication candidates:** 6.

### Top 10 connections to investigate next

1. `sleep_apnea__healthspan` — rg=+0.4005, FDR=3.08e-20, APPARENTLY_NOVEL.
2. `shortsleep__parental_lifespan` — rg=-0.2976, FDR=2.02e-16, APPARENTLY_NOVEL.
3. `snoring__healthspan` — rg=+0.2927, FDR=1.9e-14, APPARENTLY_NOVEL.
4. `sleep_apnea__parental_lifespan` — rg=-0.2617, FDR=5.58e-12, APPARENTLY_NOVEL.
5. `longsleep__parental_lifespan` — rg=-0.2395, FDR=1.33e-08, APPARENTLY_NOVEL.
6. `sleep_efficiency__frailty` — rg=-0.1520, FDR=0.00012, APPARENTLY_NOVEL.
7. `sleep_apnea__frailty` — rg=+0.4687, FDR=2.32e-46, MR_ONLY.
8. `shortsleep__frailty` — rg=+0.4543, FDR=1.57e-60, MR_ONLY.
9. `shortsleep__healthspan` — rg=+0.3554, FDR=9.26e-16, MR_ONLY.
10. `insomnia__healthspan` — rg=+0.3432, FDR=3.7e-15, MR_ONLY.

Selection rule: all conservative `APPARENTLY_NOVEL` pairs first by descending |rg|, followed by clean no-direct-rg pairs by descending |rg|. This is a transparent lexicographic rule, not a hidden score.

## 1. Scope, immutable inputs, and result reconstruction

**OBSERVED RESULT.** The locked matrix contains 396 unique sleep–external pairs: 372 primary rows and 24 explicitly segregated T2D/melanoma sensitivity rows. All 153 headline discoveries reproduce directly from the locked 396-family FDR column. Recomputing BH only across the 372 primary rows gives 155, which is reported solely as sensitivity and never substituted for the locked definition. The requested legacy `results/tables/phase1_completion.json` is absent; the canonical completion-equivalent provenance is `results/atlas/core.provenance.json`. A stale legacy Phase-1 Markdown summary was not used.

**OBSERVED RESULT.** Input verification status is `PASS_WITH_DOCUMENTED_METADATA_DISCREPANCIES`. The master table retains every requested quantitative field, source/QC label, ancestry, sample size, overlap, h² statistic, intercept, and immutable primary/sensitivity status.

**INTERPRETATION.** The atlas can be analyzed defensibly without changing the underlying Phase-1 results, provided the 153/155 metadata discrepancy and absent legacy completion filename remain explicit.

## 2. Dominant architecture

**OBSERVED RESULT.** Adverse or dysregulated self-reported sleep traits form the dominant positive cardiometabolic–psychiatric–frailty axis. The closest whole-profile pair is `insomnia__shortsleep` (Pearson r=0.953); the most opposed profiles are `sleepdur__shortsleep` (r=-0.822). PC1 explains 76.6% of descriptive profile variance and PC2 8.5%.

**INTERPRETATION.** Insomnia, short sleep, long sleep, sleepiness, napping, snoring, and sleep apnea largely share a broad adverse-health rg profile. Continuous sleep duration and objective efficiency often point in the opposite direction, so duration extremes and sleep quality should not be collapsed into one construct.

**HYPOTHESIS.** Genomic SEM should test a broad adverse-sleep factor plus measurement-specific residual factors. This descriptive PCA/clustering result is not itself a latent biological model.

## 3. Sleep and external hubs

**OBSERVED RESULT.** Insomnia, long sleep, and napping each connect to 21 external traits. Insomnia is the largest weighted hub (sum significant |rg|=4.2123); sleep apnea is most intense on average among its significant edges. Sleep timing is most specialized, with only the rheumatoid-arthritis connection passing locked FDR. Major depression is the external hub with 10 sleep connections; ADHD, BMI, parental lifespan, and frailty each have nine.

**OBSERVED RESULT.** Removing insomnia leaves 132 significant edges, with `longsleep;napping` as the leading degree hub(s) and `aging` as the leading domain.

**INTERPRETATION.** Insomnia strengthens the story but does not create it. The architecture remains distributed across long sleep, napping, short sleep, apnea, metabolic traits, psychiatric traits, and aging endpoints.

## 4. Domain architecture and cross-domain structure

**OBSERVED RESULT.** Psychiatric traits have the highest connected proportion (34/48); aging has the largest count (39). The psychiatric permutation enrichment is nominal (P=0.0080) but narrowly misses domain-level BH significance (q=0.0559). Cancer and neuro show nominal depletion signals, also not robust after domain-level correction.

**OBSERVED RESULT.** Data-driven external-trait clustering selects two broad cross-domain profile clusters (mean silhouette 0.585), rather than recovering the seven labels cleanly. Network communities are therefore cross-domain by construction of the observed rg profiles.

**INTERPRETATION.** Sleep-related genetic sharing cuts across clinical taxonomies. Psychiatric, metabolic, cardiovascular, immune, cancer, and aging labels remain useful summaries, but the rg profile does not respect them as isolated modules.

**HYPOTHESIS.** PLACO and multivariate fine-mapping should test shared loci at cross-domain bridge edges, while Genomic SEM should test whether a broad adverse-health factor explains their covariance.

## 5. Objective versus subjective sleep

**OBSERVED RESULT.** Subjective mean |rg| is materially larger than objective mean |rg| for nine external traits. The five largest objective-minus-subjective mean-|rg| differences are frailty (-0.252), healthspan (-0.185), adhd (-0.171), mdd (-0.150), parental_lifespan (-0.137). Direction discordance occurs for multiple psychiatric, metabolic, immune, and aging endpoints. Objective sleep efficiency often reverses the adverse self-report pattern; actigraphy sleep duration has only four significant links and sleep timing only one.

**INTERPRETATION.** Self-report, disease-like apnea, actigraphy duration, efficiency, and timing are not interchangeable exposures. Differences may reflect biology, measurement error, sample size, ascertainment, or all four.

**HYPOTHESIS.** Measurement-stratified Genomic SEM and local rg should test shared versus method-specific components; the current group means are descriptive and not independence-based tests.

## 6. Strong negative relationships

- `shortsleep__parental_lifespan`: rg=-0.2976, FDR=2.02e-16, APPARENTLY_NOVEL.
- `sleep_apnea__hdl`: rg=-0.2749, FDR=1.19e-26, DIRECT_RG_REPLICATION_DIFFERENT_DATASET.
- `insomnia__parental_lifespan`: rg=-0.2740, FDR=7.98e-13, DIRECT_RG_REPLICATION_DIFFERENT_DATASET.
- `sleep_apnea__parental_lifespan`: rg=-0.2617, FDR=5.58e-12, APPARENTLY_NOVEL.
- `longsleep__parental_lifespan`: rg=-0.2395, FDR=1.33e-08, APPARENTLY_NOVEL.
- `sleep_efficiency__ovarian_cancer`: rg=-0.2289, FDR=0.00402, NO_DIRECT_RG_FOUND.
- `sleepdur__frailty`: rg=-0.2220, FDR=7.71e-15, MR_ONLY.
- `snoring__hdl`: rg=-0.1893, FDR=5.03e-22, DIRECT_RG_PREVIOUSLY_REPORTED.
- `napping__hdl`: rg=-0.1892, FDR=3.83e-26, DIRECT_RG_PREVIOUSLY_REPORTED.
- `napping__parental_lifespan`: rg=-0.1793, FDR=4.01e-10, DIRECT_RG_PREVIOUSLY_REPORTED.

**INTERPRETATION.** The aging pattern is coherent: adverse sleep tends to correlate negatively with lifespan, telomere, HDL, and grip-strength phenotypes. Continuous duration and objective efficiency can instead correlate negatively with disease/frailty. Immune and cancer negatives are more heterogeneous and sensitive to phenotype/source definitions.

**HYPOTHESIS.** Signed local analyses should distinguish same-locus antagonistic sharing from mixtures of positive and negative local effects. A negative global rg alone does not establish antagonistic pleiotropy.

## 7. Published-rg replication

**OBSERVED RESULT.** Explicit prior direct-rg evidence was found for 97 pairs. Among 86 pairs with numeric exact-concept comparators, atlas rg versus median published rg has Pearson r=0.822, mean absolute difference=0.070, and direction concordance=88.4%. Comparable-SE results are mutually consistent within combined 95% uncertainty in 81.6% of rows.

**OBSERVED RESULT.** Important retained conflicts include sleep apnea–atrial fibrillation (atlas +0.1304 versus published -0.3056) and sleep apnea–triglycerides (atlas +0.2144 versus preprint -0.0118, nonsignificant). The latter source is a 2025 preprint and is labeled accordingly.

**INTERPRETATION.** The high agreement supports pipeline validity, not independent biological replication: most sources reuse or overlap public GWAS inputs. Dataset relations and phenotype matches are printed row by row.

## 8. Literature novelty audit

**OBSERVED RESULT.** Four fixed search lenses were run for every discovery (612 queries; 1,195 candidate records), then supplemented with checksum-pinned published tables and targeted full-text review. Classification counts are: APPARENTLY_NOVEL=6, DIRECT_RG_PREVIOUSLY_REPORTED=77, DIRECT_RG_REPLICATION_DIFFERENT_DATASET=20, MR_ONLY=30, NO_DIRECT_RG_FOUND=7, OBSERVATIONAL_ONLY=7, RELATED_GENETIC_EVIDENCE_ONLY=5, UNCERTAIN=1.

- `sleep_apnea__healthspan`: rg=+0.4005, FDR=3.08e-20, MODERATE_CONFIDENCE; LIMITED_BY_HEALTHSPAN_PHENOTYPE_AND_SOURCE_OVERLAP.
- `shortsleep__parental_lifespan`: rg=-0.2976, FDR=2.02e-16, HIGH_CONFIDENCE; LIMITED_EXACT_PARENTAL_LIFESPAN_ALTERNATIVES.
- `snoring__healthspan`: rg=+0.2927, FDR=1.9e-14, HIGH_CONFIDENCE; LIMITED_BY_HEALTHSPAN_PHENOTYPE_AND_UKB_OVERLAP.
- `sleep_apnea__parental_lifespan`: rg=-0.2617, FDR=5.58e-12, MODERATE_CONFIDENCE; LIMITED_EXACT_PARENTAL_LIFESPAN_ALTERNATIVES.
- `longsleep__parental_lifespan`: rg=-0.2395, FDR=1.33e-08, HIGH_CONFIDENCE; LIMITED_EXACT_PARENTAL_LIFESPAN_ALTERNATIVES.
- `sleep_efficiency__frailty`: rg=-0.1520, FDR=0.00012, HIGH_CONFIDENCE; LIMITED_BY_OBJECTIVE_SLEEP_SAMPLE.

**INTERPRETATION.** These six are “apparently unreported” as of 2026-08-28, not first-ever claims. The rule requires no located direct or pair-specific genetic evidence, |rg|≥0.15, FDR≤0.01, and no QC_CAUTION result.

**HYPOTHESIS.** Their main scientific value is as falsifiable replication targets, especially aging endpoints whose independent exact-match GWAS options are currently limited.

## 9. Priority candidates

- Rank 1: `sleep_apnea__healthspan` (TIER_A), rg=+0.4005, FDR=3.08e-20.
- Rank 2: `shortsleep__parental_lifespan` (TIER_A), rg=-0.2976, FDR=2.02e-16.
- Rank 3: `snoring__healthspan` (TIER_A), rg=+0.2927, FDR=1.9e-14.
- Rank 4: `sleep_apnea__parental_lifespan` (TIER_A), rg=-0.2617, FDR=5.58e-12.
- Rank 5: `longsleep__parental_lifespan` (TIER_A), rg=-0.2395, FDR=1.33e-08.
- Rank 6: `sleep_apnea__frailty` (TIER_B), rg=+0.4687, FDR=2.32e-46.
- Rank 7: `shortsleep__frailty` (TIER_B), rg=+0.4543, FDR=1.57e-60.
- Rank 8: `shortsleep__healthspan` (TIER_B), rg=+0.3554, FDR=9.26e-16.
- Rank 9: `insomnia__healthspan` (TIER_B), rg=+0.3432, FDR=3.7e-15.
- Rank 10: `longsleep__frailty` (TIER_B), rg=+0.2925, FDR=6.38e-14.
- Rank 11: `sleepiness__frailty` (TIER_B), rg=+0.2779, FDR=1.32e-23.
- Rank 12: `longsleep__healthspan` (TIER_B), rg=+0.2677, FDR=4.12e-10.
- Rank 13: `snoring__frailty` (TIER_B), rg=+0.2516, FDR=1.72e-21.
- Rank 14: `napping__frailty` (TIER_B), rg=+0.2506, FDR=1.37e-26.
- Rank 15: `sleepdur__frailty` (TIER_B), rg=-0.2220, FDR=7.71e-15.

The full 25-row table preserves five Tier A, ten Tier B, five Tier C, and five Tier D controls. Tier A permits an alternative-dataset replication design but explicitly records when an exact independent parental-lifespan or healthspan GWAS is limited; this caveat should travel with every claim.

## 10. Local and multivariate follow-up

**OBSERVED RESULT.** The local-analysis set contains 45 unique pairs: category A apparently novel=6, B known controls=6, C negatives=12, D measurement discordance=12, E globally null cancellation candidates=16, and F bridge edges=12. Categories can overlap.

**HYPOTHESIS.**

1. **LAVA / HDL-L:** test signed local covariance for all Tier A pairs, unusual negatives, and the 16 globally weak/null cancellation candidates.
2. **MiXeR:** test polygenic overlap despite near-zero global rg, particularly where related sleep phenotypes show strong effects.
3. **PLACO:** test pleiotropic loci at broad cross-domain bridge edges and known controls before applying it to novel candidates.
4. **Colocalization / fine-mapping:** restrict mechanistic claims to loci with compatible causal configurations and molecular-QTL support.
5. **Genomic SEM:** compare an adverse-sleep factor, a duration-extremes factor, and objective-measure residuals; confirm that factor interpretation is stable to trait omission.

No local cancellation, antagonistic pleiotropy, mediation, or causality is claimed in this report.

## 11. Findings not to emphasize

**OBSERVED RESULT.** Significant-result QC classes are HIGH_CONFIDENCE=73, MODERATE_CONFIDENCE=79, QC_CAUTION=1. T2D and melanoma comprise all 24 sensitivity rows and are excluded from every headline count, literature classification, priority table, and figure-12 status comparison.

- Do not headline T2D or melanoma; they failed upstream QC gates and remain sensitivity only.
- Do not treat the FinnGen CAD proxy as identical to narrowly defined CAD, FinnGen asthma/MS substitutions as the preferred discovery endpoints, or the NFE telomere substitution as the full mixed-ancestry release.
- Do not promote sleep-efficiency–ovarian cancer as a clean discovery: its external h² Z is low and the pair is QC_CAUTION.
- Do not call OSA–atrial fibrillation or OSA–triglycerides replicated; their reviewed directions conflict or are unstable.
- Do not interpret network betweenness, clustering, or PCA axes as causal biology.
- Do not infer causal direction from genetic correlation, MR co-occurrence, or observational consistency.

## 12. Multiple-testing and robustness sensitivity

**OBSERVED RESULT.** The immutable primary result remains 153. Sensitivity counts among the 372 primary rows are: primary-only BH=155, Bonferroni across 396=102, locked FDR≤0.01=129, and locked FDR≤0.001=103. The primary-only BH family adds two calls but does not replace the locked family.

**OBSERVED RESULT.** Omitting insomnia removes its 21 discoveries but leaves 132 and preserves aging as the leading domain. Omitting long sleep or napping similarly leaves 132. Domain omission changes the leading label when that domain is removed, but a multi-domain architecture remains.

**INTERPRETATION.** The scientific story is not a single-trait artifact, although insomnia, long sleep, and napping materially determine edge count. Source-substitution summaries are descriptive because ancestry, endpoint breadth, and ascertainment differ.

## 13. Figures and output map

All figures are generated by `scripts/103_plot_phase1_atlas.py`; exact bytes and SHA-256 hashes are recorded in the figure manifest.

1. `results/figures/phase1_deep_analysis/phase1_01_full_rg_heatmap.png`
2. `results/figures/phase1_deep_analysis/phase1_02_primary_fdr_heatmap.png`
3. `results/figures/phase1_deep_analysis/phase1_03_sleep_trait_hubs.png`
4. `results/figures/phase1_deep_analysis/phase1_04_external_trait_hubs.png`
5. `results/figures/phase1_deep_analysis/phase1_05_domain_connectivity_matrix.png`
6. `results/figures/phase1_deep_analysis/phase1_06_sleep_profile_similarity.png`
7. `results/figures/phase1_deep_analysis/phase1_07_sleep_profile_dendrogram.png`
8. `results/figures/phase1_deep_analysis/phase1_08_sleep_profile_pca.png`
9. `results/figures/phase1_deep_analysis/phase1_09_bipartite_network.png`
10. `results/figures/phase1_deep_analysis/phase1_10_negative_rg_network.png`
11. `results/figures/phase1_deep_analysis/phase1_11_objective_subjective_comparison.png`
12. `results/figures/phase1_deep_analysis/phase1_12_known_vs_novel_effect_size.png`

Core tables: `phase1_master_analysis.tsv`, `literature_novelty_audit.tsv`, `published_rg_replication.tsv`, `novel_connection_priorities.tsv`, `local_followup_priorities.tsv`, and `top_findings.tsv`. The final table contains five separate 25-row views—effect magnitude, statistical evidence, novelty, confidence, and mechanistic follow-up—with its lexicographic ranking rule printed on every row.

## 14. Literature sources and claim limits

Representative primary sources include the published insomnia–frailty analysis ([PMCID PMC10903740](https://pmc.ncbi.nlm.nih.gov/articles/PMC10903740/)), chronotype–breast-cancer analysis ([PMID 37380357](https://pubmed.ncbi.nlm.nih.gov/37380357/)), insomnia–telomere analysis ([PMID 36182724](https://pubmed.ncbi.nlm.nih.gov/36182724/)), OSA cardiometabolic study ([PMID 33243845](https://pubmed.ncbi.nlm.nih.gov/33243845/)), OSA–BMI analysis ([PMCID PMC11166156](https://pmc.ncbi.nlm.nih.gov/articles/PMC11166156/)), and the OSA multi-phenotype preprint ([PMCID PMC12642705](https://pmc.ncbi.nlm.nih.gov/articles/PMC12642705/)). The pair-complete evidence, exact phenotype/dataset relation, source URL, PMID/DOI, numeric published rg, and contradictions are in the audit and replication tables.

Search absence can be overturned by missed indexing, inaccessible supplements, newly published work, or a different phenotype definition. “Apparently unreported” is therefore dated and conditional. Genetic correlation describes genome-wide sharing, not causation or a single mechanism. Every mechanistic proposal above is a hypothesis for downstream analysis.
