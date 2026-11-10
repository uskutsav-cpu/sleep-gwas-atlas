# Two public resources: evidence for the sleep–frailty project

Checked 29 October 2026. Neither resource provides a validated evaluation of the project's 15- or 22-variant scores.

## Public FinnGen PGS benchmark and SJLIFE cross-context

[Kolosov et al., 2026](https://www.nature.com/articles/s41467-026-74461-7) evaluated 3,168 existing PGS models in FinnGen and made aggregate results available. Its [supplementary workbook](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-026-74461-7/MediaObjects/41467_2026_74461_MOESM3_ESM.xlsx) supplies directly inspectable summary data without use of the browser. Local source: `pgs_browser_2026_supplement.xlsx`, SHA-256 `10041fe7667e1f606d808a0a9d22742ff8fcf91b3faf5a6be1272babb9913dfe`.

In supplementary Data 5, row 405, the same published insomnia score **PGS001932** as in our SJLIFE grip exploration has a reported raw ROC AUC of **0.541550385** for FinnGen's `F5_INSOMNIA` endpoint, with **5,763 cases and 445,483 controls**. This is discrimination for the insomnia diagnosis in that study, not for grip or frailty. The paper describes Data 5 as classification measures at the 90th percentile; use its ROC AUC as a score-level metric and do not treat the paper's accompanying threshold-specific measures as applying to our SJLIFE analysis. The Data 1 'FinnGen samples' column is marked FALSE for this score; its description says the flag refers to Finnish samples in score development. A negative flag does not independently prove full sample and relatedness non-overlap; the paper warns of hidden overlap in some scores. Selection, endpoint ascertainment and differences between self-reported sleep and coded insomnia are important.

This result is distinct from our [SJLIFE exploratory analysis](sjlife_insomnia_pgs_grip_exploration.md), which used PGS001932 with a measured grip outcome in cancer survivors; it was not significant and did not test insomnia diagnosis. The two analyses cannot be pooled or directly compared as estimates of the same phenotype. Neither tested our custom 15- or 22-marker score. Do not describe either result as score validation for our project.

The workbook's SD1 has 13 rows whose reported/mapped traits or endpoint text contain `insomnia`, and none matching `frail` or `frailty` (literal case-insensitive search of the first five columns). This does **not** rule out indirectly relevant endpoints, phenotypes described without those words, or associations within the full PheWAS atlas. The public supplement contains matched target endpoint information, not an exhaustive table proving no insomnia-score/frailty association.

The [interactive PGS Browser](https://pgs.nchigm.org/) offers PGS-centric and endpoint-centric PheWAS searches and downloads, but presented a Terms of Use acceptance dialog; it was not accepted. The website displays CC BY-NC-ND 4.0 and asks users to contact the authors about commercial use or adaptations. If the account holder accepts the terms, search PGS001932 and other insomnia scores against specified frailty-related FinnGen endpoints, record endpoint definitions and unfiltered result counts, account for multiple testing and source overlap, and interpret as existing-score association evidence. The paper states individual-level FinnGen data are restricted; the available interface does not establish that custom score weights can be submitted for evaluation in FinnGen.

## WebGWAS: approximate phenotype GWAS, not held-out validation

[WebGWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC13170032/) is an open web tool described as using precomputed GWAS, phenotype covariances and anonymous phenotype data to estimate GWAS statistics for user-defined phenotypes. Its published implementation is limited to White British UK Biobank data, available input phenotypes and selected variants. Results for nonlinear phenotype definitions are approximations; phenotypic goodness of fit should be checked, and the method does not replace a direct GWAS. The public web fetcher returned HTTP 403, but a subsequent browser check succeeded. An exploratory broad sleep-disorder and fatigue phenotype was run; see [the screening report](webgwas_sleep_fatigue_candidate_screen.md) for its limited approximation fit, partial marker coverage, and unresolved sample-size discrepancy.

A future exploratory use would be to define a sleep-excluded frailty proxy from genuinely available components and contrast its approximate associations with a definition containing sleep, reporting fit and errors. This cannot supply independent confirmation because UK Biobank overlaps core discovery data and was used in candidate selection, and it cannot measure individual score discrimination, calibration or incremental prediction. It is not equivalent to reconstructing an externally validated FI or FFS.

## Actionable distinction

A free, approved cohort with linked genotype and phenotypes is still needed to calculate the fixed project score and estimate out-of-sample performance. The no-fee [CLS genetics route](https://cls-genetics.github.io/docs/intro.html) requires an application and approval; feasibility for the 1958 cohort is in `stage6_ncds_sleep_frailty_genetics_feasibility.md`. Continue to distinguish free controlled access from open, unrestricted participant-level data.

## Additional public PheWAS

A separate public MGI scan of PGS003326 provides insomnia-score associations with coded gait, weakness and fatigue phenotypes, but no measured frailty score or test of our prototypes. The nominal chronic-fatigue result does not survive correction for the full scan. [Full context and extracted rows](stage6_exprsweb_insomnia_frailty_phenowas_review.md).
