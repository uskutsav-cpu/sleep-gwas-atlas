# Public CHARLS psychological-frailty dataset: reproduction and scope

The authors' [Sleep-PFI CHARLS dataset, version 2](https://doi.org/10.17632/76zz3d74r6.2), is marked CC BY 4.0 and includes a processed individual-level R dataset, accompanying analysis script and supplementary document. The [published study](https://pubmed.ncbi.nlm.nih.gov/42873651/) describes baseline self-reported sleep in 2011 and psychological frailty in 2018. The data release contains 9,476 records and 29 columns. It is a subset of CHARLS processed for the published analysis, not an independently collected replication of that analysis. It contains no genotype columns and does not validate either project score. It does not provide the individual 15 psychological-frailty components needed to remove a sleep item or reconstruct the authors' original index; PFI is already binary in the released R object. Its construction and exact temporal details should be verified against the full study before new analysis.

## Reproduction

We checked the author's R script and fitted a base-R unweighted logistic model using its Qnight categories and listed Model 3 covariates. The 6–8 hour category is the reference. For comparability, the primary and sensitivity models use the same 7,037 complete records and all covariates including biomarkers; the precomputed outcome was positive in 5,191. Standard errors and confidence intervals are ordinary model-based Wald estimates; survey design and selection must be considered before population inference.

| Comparison | Adjust for depression score | Odds ratio (95% Wald CI) |
|---|---|---|
| <6 vs 6–8 hours | Yes | 1.275 (1.077–1.508) |
| <6 vs 6–8 hours | No | 1.617 (1.376–1.900) |
| >8 vs 6–8 hours | Yes | 0.901 (0.701–1.159) |
| >8 vs 6–8 hours | No | 0.930 (0.725–1.192) |

The first estimate reproduces the paper's reported short-sleep estimate of 1.27 (1.08–1.51), to rounding. The supplement reports 1.62 (1.38–1.90) for the no-depression-score short-sleep sensitivity; our estimate agrees to rounding. For the long-sleep sensitivity, the supplement reports 0.94 (0.73–1.20), while the matched-sample estimate is 0.93 (0.72–1.19). The difference is explained by the analysis sample: allowing the 24 records with missing CES-D but complete values for the remaining variables back into the model yields n=7,061 and long-sleep OR 0.938 (0.731–1.202), p=0.612; short-sleep OR 1.621 (1.380–1.905). Both match the supplement to rounding. This reproduces the published numbers but does not prove the authors used precisely the same rows. The authors' supplement notes that depressive symptoms are part of psychological frailty and treats removing the depression-score covariate as an overadjustment sensitivity. The two estimates show sensitivity to adjustment, but do not identify a causal effect or prove adjustment is correct or incorrect. The no-score comparison here deliberately keeps the same complete-case sample; it is not a sleep-item-excluded outcome analysis.

The dataset is selected, observational and overlaps other CHARLS analyses; its nighttime duration is not an insomnia diagnosis, and psychological frailty is neither a standard general FI nor Fried physical frailty. It is useful as a reproducible longitudinal observational comparison and as a demonstration of model sensitivity, not independent validation of the prior CHARLS finding or genetic score prediction. No individual rows or identifiers are reproduced in the results.

Reproduction: `work/stage6_charls_pfi_reproduction.R`; aggregate [results](stage6_charls_pfi_reproduction_results.csv); downloaded R data checksum SHA-256 `97afcc8845f543d5fcc0c2d617949099e3da6812009cec1bb2626a501e61689a`, verified against repository metadata. The source materials remain in `work/sleep_pfi_charls/`.

The available-case sensitivity is reproducible with `work/stage6_charls_pfi_complete_case_sensitivity.R`; aggregate output: [available-case estimates](stage6_charls_pfi_available_case_sensitivity.csv). The matched-sample comparison remains useful for isolating the effect of covariate adjustment from the effect of changing the sample.
