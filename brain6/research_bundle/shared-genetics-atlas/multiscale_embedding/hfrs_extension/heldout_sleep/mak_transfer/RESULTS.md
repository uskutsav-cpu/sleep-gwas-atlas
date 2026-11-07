# Stress test: transfer to another frailty definition

The existing leave-insomnia-row-out ridge model produced UKB-insomnia-based predictions of **0.590 for FI** and **0.557 for FFS** from six sleep–brain genetic correlations. Without refitting, we compared those saved numbers with the later Mak FI and modified physical-frailty (FP) correlations. Using an FFS prediction as a proxy for FP is an explicit construct transfer; these scores are related but not identical.

| Target insomnia source | Mak FI r_g (SE) | Mak FP r_g (SE) | Mean absolute error of the two ridge proxy predictions |
|---|---:|---:|---:|
| UKB symptoms | 0.641 (0.019) | 0.493 (0.026) | 0.057 |
| FinnGen R13 registry | 0.401 (0.028) | 0.237 (0.035) | 0.254 |
| MVP clinical | 0.536 (0.029) | 0.404 (0.034) | 0.103 |

The same two predictions are reused in all three rows: the model has not been given six R13 or MVP sleep–brain correlations. The R13 comparison is therefore a **source-proxy stress**, not an external model fit. Its physical-phenotype target of 0.237 is well below the reused 0.557 prediction. The result suggests that substituting a UKB-insomnia representation across source definitions is unreliable here; it does not identify why.

As descriptive context, the same two-outcome MAEs for the saved three-neighbor baseline are 0.174, 0.085 and 0.077, respectively; the training-mean baseline gives 0.419, 0.171 and 0.322. With only two dependent targets per row, prior model selection, repeated predictions, shared GWAS samples and different support, these figures do not provide a robust model ranking or a significance test. Actual test cohorts and per-variant inputs are not demonstrably independent; even UKB source byte identity is not fully proven. No model predicts individual disease risk, and no AlphaGenome outputs were used to train it.

The full six-row comparisons, exact source targets, baselines and checksums are in `transfer_comparison.tsv`, `source_summary.tsv`, and `receipt.json`. The [protocol](PROTOCOL.md) describes the comparison and the [original held-out evaluation](../RESULTS.md) explains the ridge and its earlier selection. A useful next evaluation would obtain six sleep–brain correlations from a genuinely independent insomnia source, independently ascertained frailty outcomes and matched LD/alleles, freeze training and baselines, and evaluate a sufficiently large set of traits and cohorts.
