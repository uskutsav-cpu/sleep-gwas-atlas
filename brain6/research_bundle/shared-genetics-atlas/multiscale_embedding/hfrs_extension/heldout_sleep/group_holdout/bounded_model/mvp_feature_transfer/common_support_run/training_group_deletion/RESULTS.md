# Clinical-insomnia ML: training-group deletion sensitivity

The fixed clinical-insomnia feature vector was evaluated after removing one of seven related sleep-measure groups from training at a time, always excluding the original questionnaire-insomnia group and retuning the ridge penalty within the remaining groups. The previously observed MVP target estimates were not used for tuning.

| Target | Full-training prediction | Range across seven deletion refits | Largest absolute shift | Observed MVP r_g |
| --- | ---: | ---: | ---: | ---: |
| FI | 0.654 | 0.640–0.673 | 0.019 | 0.543 |
| FFS | 0.608 | 0.498–0.625 | 0.111 | 0.395 |
| HFRS | 0.535 | 0.419–0.599 | 0.117 | 0.577 |

For FFS, all seven deletion predictions remained above the observed point estimate, although some were closer. For HFRS, the predictions spanned the observed point estimate. Removing the duration group changed FFS and HFRS chosen penalties from 1 to 10; removing the accelerometer group changed the FFS penalty to 10. These changes are a property of this small model and selected data, not proof that any group is incorrect or causally responsible. Group removal changes sample size, normalization, tuning and training support at the same time.

The seven refits reuse most training rows and are not independent observations. Their range is not a confidence or prediction interval, and it ignores uncertainty and covariance in the inputs, targets and source models. The original outcomes were already inspected; selecting a deletion based on its error would be data-driven. FI and FFS contain sleep/tiredness or exhaustion, HFRS contains a somnolence-related code, and source overlap and elevated LDSC intercepts remain unresolved. These results do not demonstrate independent generalization, a causal relationship or individual clinical prediction.

The [protocol](PROTOCOL.md), [all 21 predictions](predictions.tsv), [code](analyze.py), [input hashes and summary](receipt.json), and [independent augmented least-squares verification](verification.json) preserve the calculation. The separate implementation reproduced all 21 predictions and selected penalties, maximum absolute discrepancy 7.78×10⁻¹⁶. This checks arithmetic, not scientific validity. No AlphaGenome output was used for training.
