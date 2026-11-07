# Feature-ablation stress test of the frailty-correlation model

The existing ridge prototype is **sensitive to which sleep–brain correlation features it receives**. Repeating the same held-out-sleep-group evaluation after dropping one feature at a time gives the following pooled RMSE across 36 FI, FFS and HFRS correlation estimates:

| Dropped feature | Pooled RMSE | Change from all six features |
|---|---:|---:|
| None | 0.1108 | — |
| ADHD | 0.1195 | +0.0087 |
| Depression | 0.1003 | −0.0105 |
| Schizophrenia | 0.1347 | +0.0238 |
| Bipolar disorder | 0.1184 | +0.0076 |
| Alzheimer disease | 0.1080 | −0.0028 |
| Parkinson disease | 0.0846 | −0.0263 |

Each model was retrained, with regularization selected inside the same training-group folds; the all-feature predictions reproduced the earlier evaluation. Dropping Parkinson disease reduced this particular pooled error, while dropping schizophrenia increased it. **These are post-hoc observations, not evidence that either disease causes or protects against frailty or that one feature generalizes better.** The inputs are correlated, model selection changes when a feature is removed, and only 12 sleep measures contribute the 36 dependent outcome values. Selecting a model based on these outer results would reuse the evaluation data. Per-target RMSEs and all 252 predictions are provided in the [summary](summary.tsv) and [prediction table](predictions.tsv).

Every prediction still requires the held-out sleep trait's other measured sleep–brain genetic correlations; it does not predict an individual's disease risk or discover a molecular interaction. The datasets and participants can overlap across held-out groups, measurement-error covariance is not modeled, and trait definitions and long-sleep coding remain limitations. No AlphaGenome output was used to train the model. The [protocol](PROTOCOL.md), [code](run.py), [input receipt](receipt.json) and [independent least-squares verification](verification.json) document the computation. A useful next evaluation requires more independently sourced traits or cohorts, prespecified feature choices and baselines, and uncertainty-aware held-out tests.

A [post-hoc error decomposition](error_decomposition/RESULTS.md) shows that the apparent pooled improvement from omitting Parkinson disease is concentrated in the two related-measure groups; this is not independent validation or an argument to select that feature set.
