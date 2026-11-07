# Equal-group evaluation of the sleep–frailty prototype

The prior held-out evaluation pooled 12 sleep measures: two related groups contain three measures each and the other six groups contain one each. Giving each of the **eight groups equal weight** changes the pooled error values but preserves the descriptive ranking on the same saved predictions.

| Model | Original pooled RMSE | Equal-group RMSE |
|---|---:|---:|
| Training mean | 0.264 | 0.253 |
| Three nearest neighbors | 0.168 | 0.159 |
| Unbounded ridge | 0.111 | 0.099 |
| Bounded ridge | 0.103 | 0.094 |

For each target, equal-group RMSE is the square root of the mean of eight group-level mean squared errors; the overall value averages the 24 group–target mean squared errors. This is a **post-hoc reweighting of the same 36 predictions**, not a new train/test split or independent evaluation. It does not remove cohort overlap, make the groups independent, quantify uncertainty, or establish generalization. Model and grouping choices were already informed by these data, and each model requires six measured sleep–brain correlations to produce a prediction. The independent [training-group deletion](../training_group_sensitivity/RESULTS.md) and [input-uncertainty stress tests](../input_uncertainty/RESULTS.md) show important instability; equal weighting does not resolve it.

The [per-group losses](group_loss.tsv), [summary](summary.tsv), [protocol](PROTOCOL.md), [script](run.py) and [receipt](receipt.json) preserve the computation. An independent inverse-group-size weighted calculation reproduced all 16 summaries (maximum difference 2.78×10⁻¹⁷; [verification](verification.json)). This is an arithmetic check, not validation of the source estimates. No AlphaGenome predictions were used to train a model.
