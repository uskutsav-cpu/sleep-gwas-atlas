# Held-out sleep-trait prediction from brain correlations

This exploratory test holds out all three frailty correlations for one sleep trait at a time. A separate ridge regression for FI, FFS and HFRS learns from the other 11 sleep traits, using the held-out trait's **six already-known sleep–brain genetic correlations** as input features. The model's regularization is selected within the 11 training rows by a second leave-one-out loop; feature scaling is refitted within every inner loop. This is a stronger row-holdout check than masking a single matrix cell, but it still requires six measured associations for the test trait.

| Outcome | Training-mean RMSE | Three-neighbor RMSE | Ridge RMSE |
|---|---:|---:|---:|
| FI | 0.279 | 0.191 | 0.076 |
| FFS | 0.266 | 0.192 | 0.132 |
| HFRS | 0.217 | 0.139 | 0.066 |
| All 36 trait–outcome predictions | 0.255 | 0.175 | 0.096 |

The three-neighbor baseline was added after inspecting the initial ridge and mean results. All comparisons remain post hoc and use only **12 sleep traits**; their GWAS sources and participants overlap across features and targets, and source support, precision and covariance vary. They do not demonstrate performance in independent cohorts, generalization to a phenotype lacking brain-correlation inputs, causal mediation, or individual disease prediction. The model has no [−1,1] constraint, though no prediction fell outside that range in these 36 cases. It uses no AlphaGenome output for training.

Errors vary by trait. For example, the snoring–FFS correlation was 0.335 and the held-out prediction was 0.121; sleep-apnea–HFRS was 0.490 with a prediction of 0.334. Sleep apnea and HFRS involve overlapping FinnGen releases. These examples and the pooled RMSE cannot be interpreted as a model's clinical errors or as evidence that the source estimates are causal.

`predictions.tsv` includes all held-out outcomes, selected λ values and nearest neighbors; `nested_scores.tsv` preserves all 144 inner-grid scores. An independent augmented least-squares implementation in `verify.py` reproduced all 36 outer predictions, all inner scores, nearest-neighbor choices and summary RMSEs; `verification.json` records hashes. The [protocol](PROTOCOL.md) specifies design decisions and the source limitations. A next evaluation needs independently defined sleep and frailty cohorts, source-matched LD and covariance-aware modeling, plus a genuinely external prospective test.

A [row-exclusion sensitivity](source_sensitivity/RESULTS.md) separately removes long sleep and sleep apnea and retrains the model. Errors change by outcome but remain below a training-mean baseline on the same 11 evaluated traits; neither exclusion resolves source or sample-overlap issues.

A later [Mak frailty and source-proxy stress test](mak_transfer/RESULTS.md) reuses the saved insomnia predictions without refitting. Error increases substantially when the UKB-insomnia feature representation is used as a proxy for the FinnGen registry source; the small, overlapping comparison is not independent validation.

A [group-holdout stress test](group_holdout/RESULTS.md) holds related duration and accelerometer measures out together and changes the nested selection folds. Ridge RMSE rises relative to the earlier row holdout, although cohorts still overlap and the test is not independent validation.
