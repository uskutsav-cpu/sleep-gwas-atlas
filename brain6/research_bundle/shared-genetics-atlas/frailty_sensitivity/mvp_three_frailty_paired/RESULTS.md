# Clinical-insomnia correlations across three frailty measures: paired comparisons

On 991,707 common markers, the MVP insomnia genetic correlations were 0.5425 with FI, 0.3949 with FFS, and 0.5766 with HFRS. Using the same 200 genetic blocks to estimate uncertainty in their differences gives:

| Contrast | Difference | Paired SE | Approximate 95% interval | Bonferroni p, three contrasts |
| --- | ---: | ---: | ---: | ---: |
| FI − FFS | 0.1476 | 0.0355 | 0.0780 to 0.2172 | 9.75×10⁻⁵ |
| FI − HFRS | −0.0341 | 0.0417 | −0.1159 to 0.0477 | 1.000 |
| FFS − HFRS | −0.1817 | 0.0407 | −0.2614 to −0.1020 | 2.39×10⁻⁵ |

The two contrasts involving FFS pass the stated three-contrast correction; FI versus HFRS does not. These are **post-hoc** normal-approximation tests, and the three-test correction does not cover the broader atlas search. The non-significant FI–HFRS contrast does not establish equivalence. Point differences use the four-decimal values printed by LDSC, so p-values and intervals are approximate. The three original measures differ in ascertainment and cohort composition as well as construct: FI contains insomnia and tiredness, FFS contains exhaustion, and HFRS contains a somnolence-related diagnosis category. Consequently, the differences cannot be attributed to sleep items or physical frailty alone, and they are not evidence of causality or a shared causal allele. Unknown participant overlap, approximate sample-size conventions, ancestry/reference mismatch and elevated MVP and HFRS intercepts remain material.

All three inputs and the six other feature files were independently verified to have identical ordered SNP and allele signatures; all LDSC pairs retained the same 991,707 markers. We reconstructed each correlation's 200 delete-block values from the saved h² and genetic-covariance deletions, matched individual SEs to logged precision, and confirmed identical MVP h² delete vectors. The [protocol](PROTOCOL.md), [analysis code](analyze.py), [contrasts](contrasts.tsv), [input hashes](receipt.json), [source log and delete vectors](../../multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run/jackknife_uncertainty/RESULTS.md), and [input verification](../../multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/common_support_run/input_verification.json) document the calculation. These checks support paired computation, not source validity or cohort independence. No AlphaGenome outputs were used in model training.

The earlier [three-frailty MVP sensitivity](../../replication/mvp_frailty_genomewide/RESULTS.md) used 1,103,239 markers and did not calculate paired contrasts. The present, smaller common support was originally assembled for the ML feature transfer and should not be confused with an independently prespecified frailty-only analysis.
