# Which held-out groups drive the ablation result?

The previously observed lower pooled error after omitting the Parkinson disease feature is **concentrated in two related-measure groups**. The self-reported duration group and accelerometer group account for about **91% of the net reduction in summed squared error** across the 36 held-out FI, FFS and HFRS values. This percentage is a descriptive ratio of signed error differences; positive and negative contributions can cancel, so it is not a variance-explained measure.

| Evaluated cells | All six features RMSE | Parkinson feature omitted RMSE |
|---|---:|---:|
| All 36 | 0.1108 | 0.0846 |
| The 18 cells outside the duration and accelerometer groups | 0.0860 | 0.0806 |

The second row reuses the same already-fitted outer predictions; it is **not** a new model fit or independent validation. The six features and groups are correlated, the trait datasets may share participants, and the duration group includes long sleep, whose source-model semantics remain unresolved. The feature and groups were inspected before this decomposition, and the small number of related traits precludes a reliable importance ranking or statistical test. The result does not imply Parkinson biology is unhelpful, establish a sleep–frailty mechanism, or justify choosing a model on these evaluation results.

The [protocol](PROTOCOL.md), [all cell contributions](cell_deltas.tsv), [group totals](group_deltas.tsv), [analysis](analyze.py), [receipt](receipt.json) and [verification](verification.json) record the post-hoc arithmetic. The prior [feature-ablation report](../RESULTS.md) describes the nested training procedure and broader limitations. No AlphaGenome outputs were used in model training.
