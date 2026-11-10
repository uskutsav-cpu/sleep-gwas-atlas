# Conditional score-weight sensitivity: public All of Us fatigue associations

29 October 2026. Using the [public EUR-group malaise/fatigue candidate associations](stage6_aou_eur_fatigue_candidate_association.csv) and existing European reference LD, we applied the project's previously documented [weighted marginal-Z approximation](stage6_aou_eur_weighted_z_sensitivity_review.md) to all 15 and 22 candidate variants. The outcome `SS_807` is a broad EHR-coded malaise/fatigue phenotype; it is not FI, FFS or a joint sleep–frailty endpoint.

| Research score construction | Variants | Conditional weighted Z | Nominal two-sided p |
|---|---:|---:|---:|
| 15-marker FI weights | 15 | 3.217 | 0.00129 |
| 15-marker inferred FFS weights | 15 | 3.182 | 0.00146 |
| 15-marker insomnia weights | 15 | 1.963 | 0.0496 |
| 22-marker equal-weight count | 22 | 4.659 | 0.00000317 |
| 22-marker FI weights | 22 | 4.713 | 0.00000244 |
| 22-marker inferred FFS weights | 22 | 3.412 | 0.000645 |
| 22-marker insomnia weights | 22 | 4.573 | 0.00000481 |

These are **not regressions of scores on outcomes in people**, and the coefficients are not odds ratios, prediction accuracy, or evidence that one score is better. The scores, trait and analysis were selected after extensive project exploration; the listed p-values are unadjusted for the many tested phenotypes, score variants and methods. They should not be treated as prespecified replication. Fatigue can directly accompany insomnia and appears in physical-frailty constructs, so a signal here does not establish independent shared insomnia–frailty mechanisms. The earlier [variant-level max-p screen](stage6_aou_insomnia_fatigue_eur_conjunction_review.md) found no candidate passing BH adjustment across 33 candidates. It asks a different question from this weighted summary; these results do not contradict it.

A [two-locus removal sensitivity](stage6_aou_fatigue_two_locus_removal_review.md) shows that all seven conditional statistics remain positive after removing the insomnia-associated loci present in each score. This is a postselected change in score definition, not evidence of better prediction or a causal mechanism.

A later [AFR-group sensitivity](stage6_aou_afr_fatigue_sensitivity_review.md) uses 661 AFR-labeled reference individuals. Its conditional statistics have different uncertainty and sample sizes; no formal between-group test or individual-score validation was performed.

The approximation is `Z = aᵀz / sqrt(aᵀRa)` with `a = weight × sqrt(2f(1−f))`; `z` is signed using the source beta and p, `f` is the source allele frequency, and `R` is 503-person 1000 Genomes EUR dosage correlation with cross-chromosome correlations set to zero. It assumes the source beta and frequency refer to ALT, a biallelic diploid contrast, Hardy–Weinberg genotype variance, comparable effective sample sizes and representative target LD. The [pipeline/browser audit](stage6_hunt_aou_allele_direction_review.md) supports the ALT coding convention but does not independently authenticate every source record. Several rsIDs have multiallelic reference representations; FFS effect orientation is supported by source cross-checks but remains inferred. The reference and the target may differ in LD, and the source cohort or relatives may overlap with discovery/selection data. The coefficient of variation of `SE × HWE dosage SD` was 0.54% in the 15 and 0.23% in the 22 set; this is only a diagnostic and does not prove calibration.

The [seven results](stage6_aou_eur_fatigue_weighted_z_sensitivity.csv), [per-marker terms](stage6_aou_eur_fatigue_weighted_z_components.csv), and [receipt](stage6_aou_eur_fatigue_weighted_z_sensitivity_receipt.json) retain leave-one-marker-out ranges and provenance. Those ranges are post hoc sensitivity summaries, not confidence intervals. Reproduce using `work/stage4_env/bin/python work/aou_complete/eur_fatigue_weighted_z_sensitivity.py`. An independently genotyped cohort with appropriate sleep and frailty outcomes remains necessary for any individual-score evaluation.
