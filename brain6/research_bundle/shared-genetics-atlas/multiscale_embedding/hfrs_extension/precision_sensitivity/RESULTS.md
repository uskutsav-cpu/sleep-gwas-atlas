# Precision-weighted embedding sensitivity

Using the same 50 masks and the same five initialization seeds per fold as the unweighted 12×9 experiment, we fitted a rank-two model with weights derived from the reported training-cell standard errors. The weights were clipped to [0.25, 4]; this is an exploratory sensitivity, not a validated uncertainty model. The final same-seed comparison was specified after an initial run with different seeds had been inspected.

| Masked outcome | Unweighted rank-two RMSE | Precision-weighted rank-two RMSE |
|---|---:|---:|
| All nine outcomes | 0.06976 | 0.06776 |
| FI | 0.06906 | 0.06505 |
| FFS | 0.07754 | 0.08143 |
| HFRS | 0.06530 | 0.06576 |

The overall change is small; FFS and HFRS do not improve under this sensitivity. On the same masks the weighted additive baseline has pooled RMSE 0.14648 versus 0.14628 unweighted, and weighted column means have RMSE 0.18747 versus 0.18763 unweighted. These are dependent repeated reconstructions of only 108 source estimates; neither the pooled differences nor the apparent rank-two advantage is an independently validated model comparison. No AlphaGenome outputs were used as model features.

For context, the root mean squared error divided by each held-out estimate's SE is 2.45 for the unweighted rank-two fit and 2.43 for the weighted fit. This is a descriptive scale comparison, not a chi-square statistic or calibrated uncertainty measure. The estimates share cohorts and have unmodeled covariance, use heterogeneous source methods and support, and some displayed SEs are rounded. Diagonal weighting cannot correct these limitations, nor can it validate a causal pathway, individual risk prediction or generalization to unseen traits.

`predictions.tsv`, `summary.tsv`, `training_objectives.tsv` and `fold_weights.tsv` preserve the paired comparison. `verify.py` checks mask identity, best training-objective selection and metric arithmetic; `verification.json` and `receipt.json` provide hashes and settings. The [protocol](PROTOCOL.md) and [unweighted experiment](../RESULTS.md) record the design and open QC issues.
