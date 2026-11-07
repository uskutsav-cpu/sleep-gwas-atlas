# Source-related sensitivity of held-out-sleep predictions

We repeated the leave-one-sleep-trait-out ridge evaluation after excluding one additional sleep trait from all training and evaluation: long sleep, whose source/model semantics remain unresolved, or FinnGen sleep apnea, whose cohort overlaps FinnGen HFRS. Each scenario tests the other 11 sleep traits and trains on 10 per outer fold. Nested regularization and preprocessing are repeated from scratch. The original model's performance is recalculated on those **same 11 test traits** for a paired descriptive comparison.

| Excluded trait | Outcome | Original-model RMSE on same 11 traits | Exclusion-model RMSE | Exclusion training-mean RMSE |
|---|---|---:|---:|---:|
| Long sleep | FI | 0.071 | 0.084 | 0.292 |
| Long sleep | FFS | 0.129 | 0.133 | 0.273 |
| Long sleep | HFRS | 0.062 | 0.071 | 0.225 |
| Sleep apnea | FI | 0.079 | 0.086 | 0.277 |
| Sleep apnea | FFS | 0.131 | 0.125 | 0.255 |
| Sleep apnea | HFRS | 0.051 | 0.038 | 0.192 |

On this small, reused set the exclusion-model RMSE remains below its mean baseline; changes relative to the original model vary by outcome. Excluding apnea also removes apnea from evaluation, so the lower HFRS error does not show that the model handles apnea better. Neither exclusion eliminates shared participants from other studies, establishes compatible GWAS inputs or resolves long-sleep semantics. Hyperparameters and these sensitivities were chosen after earlier analyses; no statistical significance, independent-cohort generalization, causal mechanism or clinical prediction is established. The held-out sleep trait still supplies six observed brain-correlation features, and no AlphaGenome outputs are used in training.

`predictions.tsv`, `nested_scores.tsv` and `summary.tsv` retain all 66 paired comparisons and 264 tuning scores. An independently expressed augmented least-squares check (`verify.py`) reproduced every nested score and prediction; `verification.json` records hashes. The [protocol](PROTOCOL.md) and [original held-out-sleep analysis](../RESULTS.md) describe source limitations and selection.
