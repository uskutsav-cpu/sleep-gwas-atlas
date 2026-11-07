# Source-exclusion sensitivity for the bounded ML model

Removing the long-sleep trait from both training and evaluation increased pooled RMSE on the remaining 33 cells from **0.0956 to 0.1157**. Removing the FinnGen sleep-apnea trait instead changed pooled RMSE from **0.1017 to 0.1001** on its own 33 retained cells, with changes in different directions across the three outcomes.

| Excluded trait | Target | Original bounded RMSE on same 11 traits | Refit exclusion RMSE | Refit training-mean RMSE |
|---|---|---:|---:|---:|
| Long sleep | FI | 0.073 | 0.114 | 0.299 |
| Long sleep | FFS | 0.131 | 0.142 | 0.286 |
| Long sleep | HFRS | 0.070 | 0.083 | 0.234 |
| Sleep apnea | FI | 0.083 | 0.080 | 0.282 |
| Sleep apnea | FFS | 0.140 | 0.123 | 0.266 |
| Sleep apnea | HFRS | 0.068 | 0.092 | 0.199 |

These are post-hoc comparisons using the same small, overlapping dataset; the two exclusion scenarios have different test sets and are not a head-to-head comparison. Removing a trait also changes the training sample, some group sizes and nested lambda selection. The changes therefore do **not** identify an error in either source, estimate its causal contribution, establish statistical significance or prove that excluding it improves generalization. Removing apnea also removes it from evaluation. Other source overlaps remain, and the long-sleep released-model and per-variant sample-size semantics remain unresolved. The six existing sleep–brain correlations are still required for every test trait.

The [protocol](PROTOCOL.md), [66 paired predictions](predictions.tsv), [summary](summary.tsv), [code](run.py) and [input hashes](receipt.json) preserve the grouped-holdout comparison. An independent augmented least-squares implementation reproduced all 66 refits and nested selections (maximum difference 2.53×10⁻¹⁵; [verification](verification.json)). This verifies arithmetic, not source validity, independent-cohort performance or clinical utility. No AlphaGenome outputs were used for training. The earlier [unbounded row-holdout exclusions](../../../source_sensitivity/RESULTS.md) used a different split and should not be interpreted as the same experiment.
