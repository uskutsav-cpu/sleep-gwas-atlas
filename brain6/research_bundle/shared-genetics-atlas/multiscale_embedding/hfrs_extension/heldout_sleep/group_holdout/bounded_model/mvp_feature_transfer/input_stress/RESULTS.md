# Clinical-insomnia feature uncertainty: fixed-model stress test

Holding the trained model and its selected hyperparameters fixed, we varied all six MVP–disease genetic-correlation inputs within their reported estimate ± one LDSC standard error and calculated the extrema of the bounded prediction:

| Outcome | Point prediction | Hypothetical one-SE input-box range | Observed target point estimate |
| --- | ---: | ---: | ---: |
| Frailty index | 0.659 | 0.573–0.731 | 0.553 |
| Fried frailty score | 0.607 | 0.545–0.662 | 0.400 |
| Hospital frailty risk score | 0.522 | 0.420–0.610 | 0.570 |

Under these deliberately simultaneous bounds, the FFS prediction remains above its observed point estimate. Mathematically, reaching that fixed target with the same fitted model would require a common box of approximately **3.02 times the reported input SEs**; FI would require 1.21 and HFRS 0.52. These quantities do not include the target's own SE. The selected extreme directions for these calculations leave all six input correlations within [−1,1].

These are **stress ranges, not confidence intervals** or probabilities. The six input estimates can be correlated, and all components moving to their individually most favorable directions may be implausible. They ignore uncertainty in the outcome, training data, model coefficients and lambda choice. The results do not establish significance or prove the reason for prediction error. The training matrix is small and overlapping, and other model and source limitations remain in the [transfer analysis](../RESULTS.md).

The [protocol](PROTOCOL.md) was written before calculation; [ranges](ranges.tsv), [code](run.py) and [input hashes](receipt.json) preserve the calculation. All 64 corners per outcome independently reproduce the analytic extrema to within 1.2×10⁻¹⁶. The saved point predictions were independently reproduced to 1e-10. These checks verify arithmetic, not statistical calibration or clinical utility. No AlphaGenome output was used for training.
