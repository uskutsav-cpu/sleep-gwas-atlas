# FI and FFS: paired genetic correlations across 12 sleep traits

These are post-hoc comparisons of archived LDSC estimates, with FI and FFS matched to each sleep trait on the same SNP support and 200 jackknife blocks. Each row uses its own support and phenotype sources; the rows are not an independent cohort or a direct comparison among sleep traits. The difference is **FFS minus FI**. The confidence intervals are nominal 95% intervals, not multiplicity-adjusted.

| Sleep trait | FI r_g (SE) | FFS r_g (SE) | Difference (95% CI) | Raw p | Padded BH q396 |
|---|---:|---:|---:|---:|---:|
| Insomnia | 0.640 (0.024) | 0.483 (0.026) | -0.156 [-0.213, -0.100] | 6.21e-08 | 2.46e-05 |
| Sleep duration | -0.222 (0.028) | -0.162 (0.026) | +0.060 [+0.009, +0.111] | 0.0204 | 1 |
| Short sleep | 0.451 (0.026) | 0.440 (0.027) | -0.011 [-0.059, +0.037] | 0.646 | 1 |
| Long sleep | 0.290 (0.038) | 0.355 (0.031) | +0.065 [-0.005, +0.135] | 0.0679 | 1 |
| Chronotype | -0.039 (0.025) | -0.046 (0.025) | -0.008 [-0.062, +0.046] | 0.779 | 1 |
| Daytime sleepiness | 0.277 (0.026) | 0.227 (0.028) | -0.049 [-0.102, +0.003] | 0.0648 | 1 |
| Napping | 0.248 (0.022) | 0.233 (0.022) | -0.015 [-0.061, +0.031] | 0.528 | 1 |
| Snoring | 0.252 (0.025) | 0.335 (0.026) | +0.083 [+0.029, +0.137] | 0.00273 | 0.54 |
| Sleep apnea | 0.471 (0.033) | 0.513 (0.028) | +0.042 [-0.023, +0.108] | 0.203 | 1 |
| Sleep efficiency | -0.153 (0.037) | -0.178 (0.039) | -0.026 [-0.102, +0.051] | 0.512 | 1 |
| Accelerometer sleep duration | -0.051 (0.033) | -0.015 (0.033) | +0.037 [-0.031, +0.104] | 0.289 | 1 |
| Sleep timing | 0.069 (0.042) | -0.042 (0.039) | -0.111 [-0.204, -0.018] | 0.0195 | 1 |

Only insomnia crosses both the archived padded-family BH and Bonferroni thresholds of 0.05. Snoring has raw p=0.00273 but padded BH q=0.540, so it is not a corrected difference. Non-significant contrasts do not establish equal genetic correlations. The q-values apply to **differences**, not to tests that either individual correlation is nonzero. The archived correction used the 12 observed p-values plus 384 p=1 placeholders, not 396 observed contrast tests; it is a conservative padded family, not a joint FDR over all exploratory analyses.

FI includes insomnia and tiredness, and FFS includes exhaustion; overlapping UK Biobank participants and different constructs remain relevant. These are genetic correlations, not phenotypic correlations, estimates of causal effects or item-excluded analyses. HFRS is not included in this particular twelve-trait comparison; the [three-frailty insomnia sensitivity](../../global_audit/brain6_frailty_covariance/three_frailty_common_support/RESULTS.md) uses a different shared support.

The [source table](../../evidence/ALL12_FI_FFS_CONTRAST_RECOMPUTED.tsv), [global audit](../../global_audit/GLOBAL_RECHECK.md), [paired comparison context](../FRAILTY_COMPARISON.md) and [script](build.py) preserve the input and arithmetic. Verification recalculates the difference, normal-tail p-values and padded BH values from the stored results; it does not independently rerun LDSC or validate cohort overlap and phenotypes.
