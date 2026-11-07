# Held-out groups of related sleep measures

When related sleep measures were held out together, the ridge model's pooled reconstruction RMSE rose from **0.0959** in the earlier leave-one-row-out test to **0.1108** on the same 36 trait–frailty outcomes. This is a post-hoc stress test, not a formal estimate of a leakage penalty: both the outer training set and the inner regularization-selection folds changed. Ridge still had lower RMSE than the training-mean and three-nearest-neighbor baselines in this small dataset.

| Frailty target | Training mean RMSE | Three-neighbor RMSE | Group-held-out ridge RMSE | Earlier row-held-out ridge RMSE |
|---|---:|---:|---:|---:|
| FI | 0.2855 | 0.1841 | 0.0908 | 0.0763 |
| FFS | 0.2776 | 0.1819 | 0.1494 | 0.1319 |
| HFRS | 0.2245 | 0.1346 | 0.0792 | 0.0662 |
| All 36 values | 0.2639 | 0.1684 | 0.1108 | 0.0959 |

Each group-held-out model requires the held-out sleep measure's six **already observed** sleep–brain genetic correlations as features. The split holds out the three self-reported duration phenotypes as a group and the three accelerometer-derived phenotypes as another; other measures are singleton groups. The models were retrained and regularization was selected by nested leave-one-training-group-out validation. There are only 12 sleep measures, the source studies and participants overlap across groups, and the grouping was chosen after prior results. This is **not** independent-cohort or external validation, individual risk prediction, proof of sleep–frailty biology, or evidence of benefit from AlphaGenome features. The model does not account for the covariance or uncertainty of the source r_g estimates; long-sleep coding and HFRS ancestry/reference caveats remain.

The prior row-held-out ridge used inner leave-one-row-out selection; the group-held-out model uses inner leave-one-group-out selection. The change in RMSE therefore combines different training and model-selection regimes. No significance test or interpretation of a precisely measured improvement or loss is warranted. Error estimates reuse correlated sleep traits and target datasets, and the tiny sample makes performance unstable. None of the 36 predictions was outside [−1,1], although the model is not constrained to that interval.

The [protocol](PROTOCOL.md), [code](run.py), [36 predictions](predictions.tsv), [nested scores](nested_scores.tsv), [summary](summary.tsv), and [input receipt](receipt.json) preserve the analysis. An [independent augmented least-squares implementation](verify.py) reproduced all outer predictions, inner scores, neighbor selections and summary RMSEs; the [verification receipt](verification.json) records hashes. These are implementation checks, not evidence of independent generalization. The next meaningful validation requires source-audited independent cohorts or populations, adequately sized training data, feature uncertainty and covariance, and prespecified baselines and evaluation.

A later [feature-ablation stress test](feature_ablation/RESULTS.md) shows that pooled error changes when individual brain-correlation inputs are removed. Its small, post-selected evaluation does not assign biological importance to a feature or establish external generalization.
