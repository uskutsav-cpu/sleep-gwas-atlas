# Exploratory embedding with FI, FFS and HFRS

This extension adds 12 sleep–HFRS correlations to the previous 12×8 sleep–brain–frailty matrix, yielding 108 observed sleep–outcome entries. For all 12, the HFRS rows' recorded sleep dataset version and sample size match the existing source crosswalk. The full original input bytes, per-pair SNP support, underlying covariances and cohort overlap remain incompletely established. HFRS values come from separate source-specific LDSC runs, not the FI/FFS matched-support calculation. The FinnGen apnea and FinnGen HFRS samples overlap; the FI, FFS and several sleep sources involve UK Biobank. HFRS univariate LDSC intercept and the long-sleep model semantics require further review. These entries should be interpreted as an exploratory heterogeneous matrix.

A two-dimensional SVD accounts for 97.04% of the matrix's **in-sample squared numerical magnitude**, with in-sample reconstruction RMSE 0.0357. This is not a percentage of disease risk or biological variance explained. `coordinates.tsv` contains in-sample sleep and outcome coordinates; axis sign and rotation are arbitrary, and coordinate proximity does not identify a causal pathway. `reconstruction.tsv` retains all 108 in-sample residuals.

In a post-hoc ten-repeat, five-fold cell-masking diagnostic, the five-initialization rank-two model had lower pooled RMSE than either simple baseline:

| Model | All 1,080 masked-cell predictions | FI | FFS | HFRS |
|---|---:|---:|---:|---:|
| Outcome-column training mean | 0.1876 | 0.2848 | 0.2608 | 0.2159 |
| Additive sleep-plus-outcome model | 0.1463 | 0.1713 | 0.1566 | 0.1103 |
| Rank-two, best of five training initializations | 0.0698 | 0.0691 | 0.0775 | 0.0653 |

Each outcome column contributes 120 predictions, but the 10 masks reuse each of its 12 input estimates. Sleep traits and outcomes appear in the training folds, and model hyperparameters were informed by earlier experiments. Therefore this is **cell reconstruction**, not 1,080 independent tests, validation on an unseen trait or cohort, clinical prediction, or evidence of a mechanistic interaction. It uses no AlphaGenome outputs for training. None of the predictions exceeded [−1,1] in these masks, but the model has no such constraint.

Performance varies across cells. One held-out long-sleep–schizophrenia estimate of +0.2889 was predicted as −0.2483 in one fold; the long-sleep source has an unresolved model/metadata caveat. For insomnia–HFRS, ten dependent predictions averaged 0.441 against the stored estimate 0.439, with a 0.369–0.491 range; agreement on this one cell does not establish out-of-sample generalization. The original 12×8 model and this 12×9 version use different cells and masks and should not be compared as a controlled incremental evaluation of HFRS.

`source_matrix.tsv`, `masked_predictions.tsv`, `model_summary.tsv`, `training_objectives.tsv`, `coordinates.tsv`, `reconstruction.tsv` and `receipt.json` preserve data and results. `verify.py` checks mask coverage, training-objective selection and RMSE arithmetic; the [protocol](PROTOCOL.md) describes the choices. Before predictive use, establish source identity and ancestry-matched support, model standard errors and covariance, freeze decisions, benchmark against simple models, and test new traits or genuinely independent cohorts. The embedding does not repair either project's failed local-correlation QC.

A [precision-weighting sensitivity](precision_sensitivity/RESULTS.md) modestly changes pooled RMSE and worsens FFS reconstruction. It uses diagonal SE-based weights and does not model covariance or provide independent validation.

A subsequent [held-out-sleep test](heldout_sleep/RESULTS.md) predicts all three frailty correlations from six known brain correlations for each excluded sleep trait. Ridge outperformed training-mean and three-neighbor baselines on the 12 traits, but cohorts overlap and the additional baseline was selected post hoc; this is not independent-cohort validation.

A [related-measure group holdout](heldout_sleep/group_holdout/RESULTS.md) provides a stricter post-hoc split for the ridge model; error rises relative to the previous row holdout, and cohort overlap remains.
