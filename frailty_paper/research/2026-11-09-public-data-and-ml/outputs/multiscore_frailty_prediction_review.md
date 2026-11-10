# Comparing frailty polygenic scores and multivariate models

The best existing comparison uses the **same observed FI outcome within each cohort**, rather than comparing an FI score's R² in one cohort with an FFS score's AUC in another. In [Foote et al.](https://www.nature.com/articles/s41588-025-02269-0), published Supplementary Table 44 reports the following incremental FI variance explained above covariates; values are rounded as in that table.

| Evaluation cohort | FI-derived score | Fried-phenotype-derived score | Seven-factor multi-PRS model |
|---|---:|---:|---:|
| ELSA | 0.02 | 0.01 | 0.02 |
| PISA | 0.02 | 0.02 | 0.03 |
| LBC1936 Wave 1 | 0.02 | 0.03 | 0.02 |

These are **R² values**, so 0.02 denotes approximately 2% of observed FI variance in the specified evaluation model, not a 2% chance of frailty. The study's main text describes the combined model as strongest in ELSA and PISA and comparable to FI-PRS in LBC1936. The rounded table values alone do not demonstrate a precise ranking: equal displayed values can conceal differences, and without a paired comparison they do not establish statistically significant superiority. A previously documented discrepancy between some PGS Catalog values and Table 44 remains unresolved; this table reproduces the primary paper's supplement, not the catalog. The article's 14 November 2025 change history describes a license change, not an explanation for the discrepancy.

The multi-PRS combines signals from seven latent frailty factors using multiple regression. The paper says the general factor includes insomnia as one of its indicators. These analyses show that information from several genetic dimensions can be associated with a broad FI; they do not demonstrate that adding an insomnia score specifically improves prediction of a sleep-excluded frailty outcome. The excerpted methods do not establish that a single frozen multi-score model was evaluated unchanged in a wholly new independent cohort, so do not equate these R² values with fully held-out deployment performance. PISA and ELSA FI outcomes include sleep-related items; LBC1936 had no explicit sleep keyword among the inspected FI labels, though indirect overlap is possible. The three cohorts and published estimates do not validate this project's 15-marker candidate sums.

The [nine-row source transcription](foote_same_outcome_multiscore_comparison.csv) retains coefficients, SEs, confidence intervals, q-values, source row numbers and the reported R². Both FI and Fried-derived scores have positive associations with observed FI in all three Table 44 cohorts. For a physical frailty outcome, the same paper reports frail versus non-frail ORs per SD of 1.47 (FI score) and 1.46 (Fried score) in LBC1936; these are not a paired test of equality and are not absolute risks. Cohort-specific outcome definitions, sampling and known sample-metadata discrepancies are described in `fi_ffs_pgs_published_validation_audit.md`.

## A separate multi-score research example

A [Flint et al. preprint, version posted July 2025](https://doi.org/10.1101/2024.05.31.24308260) considered 26 **polygenic scores**, not 26 SNPs, to model frailty in ELSA and LBC1936. It describes elastic-net modeling with repeated cross-validation and reports up to 4.7% variance explained compared with a best single score of 2.5%. Its traits included chronic pain, depressive symptoms, educational attainment and others; the preprint reports that the selected scores differed in their associations with frailty. These are descriptive, model-dependent findings from a non-peer-reviewed preprint, not evidence that changing those traits alters frailty or that insomnia-specific information improves the model. The two studies reuse ELSA and LBC1936, so their reports must not be counted as independent cohort replications of one another. It is also necessary to distinguish selecting or refitting predictors within the evaluation cohort from applying an unchanged model trained elsewhere.

## Implication for the insomnia–frailty project

A defensible assessment would use the same independent participants for a covariate-only baseline and frozen FI, FFS, insomnia and sparse-candidate score models. For a predeclared FI outcome and an FFS outcome, assess each score and the added contribution of the insomnia score; evaluate a sleep- and fatigue-excluded FI as a separate sensitivity. Training, threshold/weight selection and tuning must precede a held-out test and account for overlap and relatedness with discovery GWAS. Compare model performance with paired uncertainty, using out-of-sample R² for continuous outcomes and discrimination and calibration for prespecified binary or time-to-event outcomes. Do not add raw coefficients across differently scaled trait-specific scores. None of these outcome-based evaluations was performed on this project's candidate prototype.

## Traceability

`foote_same_outcome_multiscore_comparison.csv` was transcribed programmatically from Supplementary Table 44, rows 11–13, 21–23 and 31–33. The publisher workbook SHA-256 is `d170d2596ed43cb168d4463a722894969198aa6d8e493f73cdf7cd6d5637b5b4`; the resulting CSV SHA-256 is `875dc0a2ac12798615e17500d5e826a072e4e09e5417486f46eb6677b1487c83`. We did not access individual-level cohort data or refit the published models.
