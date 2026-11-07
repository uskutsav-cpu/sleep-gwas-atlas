# Joint embedding stability and error audit

The initial rank-two model's pooled RMSE of 0.1884 hid both relatively small typical errors and severe failures. Median absolute error was 0.0412, but **one of 50 splits** (repeat 9, fold 0) contributed 83.7% of its total squared error. All four predictions outside [−1, 1] arose in that split; their worst predicted value was −2.56 for a source estimate of +0.290. Repeated splits reuse source data, so these are diagnostics, not independent tests.

A post-hoc multi-start experiment reused the same masks and factorization settings, fit five deterministic initializations per split, and selected the lowest **training** ridge objective without using held-out values for selection. It reduced pooled masked-cell RMSE to 0.0661 and produced no out-of-range predictions in these particular splits:

| Method | Pooled RMSE | Predictions outside [−1,1] |
|---|---:|---:|
| Outcome-column mean | 0.1883 | 0 |
| Original single-start rank two | 0.1884 | 4 |
| Additive sleep-plus-outcome | 0.1523 | 0 |
| Five-start rank two, training-objective selection | 0.0661 | 0 |

For FI and FFS separately, the multi-start RMSEs were 0.0757 and 0.0850, compared with 0.1943 and 0.1738 for the additive baseline. These 120 predictions per outcome reuse just 12 correlations each. The result suggests that optimization instability affected the original factorization; it does **not** demonstrate independent-cohort performance, a validated biology embedding or generalization to new sleep traits. The improvement and this remedy were examined after viewing the original results, and no new held-out source or prospective evaluation was used. Even when no predictions exceeded the range in this particular run, the algorithm has no mathematical range constraint and does not guarantee a valid positive-semidefinite genetic correlation matrix. Missing harmonization, uncertain long-sleep modeling, sample overlap and ignored covariance remain.

As a separate post-hoc descriptive check, clipping the original rank-two predictions to [−1,1] reduces pooled RMSE to 0.1301. Clipping changes the predictor after its errors were seen; it should not be presented as a preregistered improvement. For full details on ordinary and SE-scaled errors, including six brain outcomes, FI, FFS, and excluding long sleep from evaluation, see `error_audit.tsv`. Dividing residuals by reported source SE does not include prediction uncertainty or cross-trait covariance and is not a calibrated statistical test. The exclusion analysis removes long-sleep cells **from evaluation only**; the original training folds still contain other long-sleep cells, so it is not a refit without long sleep.

`MULTISTART_PROTOCOL.md` documents the pre-run selection rule; `multistart_training_objectives.tsv` contains all 250 training fits and choices, `multistart_predictions.tsv` contains all 960 held-out predictions, and `multistart.py` reproduces them. `ERROR_AUDIT_PROTOCOL.md`, `audit_errors.py`, and `error_audit_receipt.json` bind the separate error analysis to its inputs. No AlphaGenome outputs were used for model training. For a stronger next step, freeze source identity and variant support, use suitable uncertainty and covariance, constrain the model appropriately, and test on independent study estimates or entire held-out traits with features that can support prediction for unseen traits.
