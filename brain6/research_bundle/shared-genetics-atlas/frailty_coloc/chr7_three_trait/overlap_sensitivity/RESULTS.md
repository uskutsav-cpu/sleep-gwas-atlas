# Chromosome 7 three-trait correlation stress test

The previously reported FI–insomnia cluster is sensitive to the assumed correlation between the FI and insomnia association statistics. In the selected 1,055-variant window, illustrative positive correlations of 0.1, 0.3 and 0.5 returned an insomnia–ADHD cluster instead. **These correlations are hypothetical, not estimates of the actual overlap**, so the exercise does not identify which cluster is correct. It strengthens the case for treating the original independence-based FI–insomnia result as provisional.

| Assumed effective FI–insomnia correlation | Returned cluster | Model-conditional posterior | Candidate variant |
|---:|---|---:|---|
| −0.3 | FI, insomnia | 0.9858 | rs12671330 |
| −0.1 | FI, insomnia | 0.9805 | rs12671330 |
| −0.001 | FI, insomnia | 0.8849 | rs12671330 |
| 0, independence baseline | FI, insomnia | 0.8822 | rs12671330 |
| +0.001 | FI, insomnia | 0.8795 | rs12671330 |
| +0.1 | insomnia, ADHD | 0.6622 | rs9969232 |
| +0.3 | insomnia, ADHD | 0.6622 | rs9969232 |
| +0.5 | insomnia, ADHD | 0.6622 | rs9969232 |

Posteriors and candidates are package output conditional on a selected region, the chosen model and priors, not calibrated probabilities of causal variants or corrected discoveries. Positive settings yield the same reported insomnia–ADHD result because the returned pair excludes FI and the assumed ADHD correlations are zero; they are not three independent replications. The two near-zero scenarios approach the baseline, checking that the switch to the correlated-trait code path does not by itself create a large jump near zero. The actual switch point was not estimated, and no conclusion should be drawn from the negative values about the real sign.

The package multiplies `trait.cor` by `sample.overlap`; we set the latter to a matrix of ones and specify the working product directly. This does not establish complete sample overlap. It leaves ADHD-related residual correlations at zero and uses the package default identity LD, original single-signal model, fixed prior, and the prior selected variant support. FI includes sleep-related items, and actual cohort overlap, residual covariance, population matching and the FI meta-analysis allele convention need investigation. A more credible analysis requires cohort and phenotype information or an appropriate estimate of the summary-statistic covariance, followed by method validation and sensitivity to LD and multiple signals. The test does not resolve the failed LAVA gate, identify sleep-independent frailty, or establish causation.

The frozen input on the SSD matched the parent SHA-256 `9f9e631c4220437561695c2aef28ff174907b432e6919506e974f33fac19fcf6`. The [protocol](PROTOCOL.md), [run code](run.R), [log](run.log), [results](clusters.tsv), [serialized model output](sensitivity_results.rds), and [verification](verify.R) are provided. Verification reproduces the parent baseline and checks the near-zero comparison; it is a software and continuity check, not statistical calibration. See the [parent pilot](../RESULTS.md) and [HyPrColoc method](https://doi.org/10.1038/s41467-020-20885-8).

A later [genome-wide LDSC follow-up](../../../global_audit/brain6_frailty_covariance/RESULTS.md) estimated a positive original insomnia–FI cross-trait intercept. That intercept is not a validated local residual-correlation parameter and should not be substituted directly into these illustrative models.
