# Illustrative sample-size sensitivity for a genetic-score association

This calculation asks a narrow planning question: if a continuous, approximately normal frailty outcome were regressed on one fixed genetic score and ten other covariates in an unrelated complete sample, how often might a two-sided 5% test detect an incremental association? It does **not** estimate the actual effect of either project score, power to establish clinical utility, calibration, binary-outcome AUC, a joint sleep–frailty outcome, or a paired comparison of scores.

| Analyzed participants | Assumed partial R² 0.1% | Assumed partial R² 0.5% | Assumed partial R² 1% | Partial R² for ~80% power |
|---:|---:|---:|---:|---:|
| 401 | 9.6% | 28.7% | 50.8% | 1.98% |
| 922 | 15.9% | 57.0% | 85.8% | 0.86% |
| 4,006 | 51.6% | 99.4% | >99% | 0.19% |

The sample sizes are **illustrative**, chosen to show how possible cohort size affects precision. The RADC availability report flags 401 MAP participants with a sleep questionnaire and at least one genetic dataset, and 922 with actigraphy and one genetic dataset; these are not verified analytic sample sizes and overlap each other. The NHATS documentation reports 4,006 genotyped participants; it does not establish 4,006 with every desired phenotype, marker and permission. Missing outcomes, genetic QC, relatedness, selection, unequal weights, clustering, ancestry adjustment, testing multiple outcomes, variant uncertainty and lower outcome variance can materially reduce effective sample size or power. NHATS is controlled access; cost and eligibility for any specific project must be confirmed.

Partial R² here is the proportion of remaining outcome variation after covariate adjustment accounted for by the additional score, in this idealized model. The true value is unknown; these scenarios are **not** effect estimates or a forecast that a cohort will validate the score. Nor does 80% power guarantee a precise estimate or useful discrimination. A small cohort with an uncertain result is not evidence that the score has no effect, and a significant association does not establish calibrated or improved individual prediction. Selection of candidate markers and overlap with target cohorts remain separate threats even at large sample sizes.

For reproducibility, `Rscript work/stage6_score_validation_idealized_precision.R` calculates the one-degree-of-freedom noncentral-F test at alpha 0.05 with residual df = n−10−2 and illustrative noncentrality (n−10−1)R²/(1−R²). It finds the scenario R² giving 80% power numerically; [results](stage6_score_validation_idealized_power.csv). This is a simplified planning approximation, not a cohort-specific power analysis. [MAP flags](stage6_rosmap_sleep_frailty_genetics_feasibility.md); [NHATS access and counts](stage6_nhats_genetic_validation_feasibility.md).
