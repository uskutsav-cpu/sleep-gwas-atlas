# PLACO+ numerical and null-model stress checks

These are simulated data, not observations or replication. We used the official PLACO+ function and the nuisance parameters estimated from the 1,138,262-marker common support. Each pair/scenario has 1,000 independent simulations, using a two-sided 0.05 threshold. The correlation between simulated Z errors was set to the empirical PLACO+ correlation parameter for that pair. This setup does not reproduce LD, genome-wide genetic architecture, real sample overlap, phenotype definitions, or full parameter-estimation uncertainty. No conclusions about calibration at the genome-wide threshold can be drawn from 1,000 simulations.

| Simulated setting | FI: fraction p < 0.05 | FFS: fraction p < 0.05 |
|---|---:|---:|
| Both traits unit-variance, mean zero | 0.015 | 0.023 |
| Sleep Z has inflated variance matching fitted VarZ; other unit variance | 0.038 | 0.034 |
| Frailty Z has inflated variance matching fitted VarZ; other unit variance | 0.034 | 0.041 |
| Sleep has fixed Z mean +2; frailty has zero mean | 0.177 | 0.173 |
| Frailty has fixed Z mean +2; sleep has zero mean | 0.187 | 0.136 |

The first three rows match elements of the method's zero-mean, inflated-variance approximation; the last two deliberately challenge it with a fixed effect for one trait. They show elevated rejection in that stress case at 0.05, not that PLACO+ is generally invalid or that these empirical rates apply to the real data. The [method paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12547305/) reports calibration under its study designs and models and discusses the approximation for nonzero single-trait effects. It should be consulted before interpreting these simulations. The inflated-variance controls and the fixed-mean stress test together indicate that any scan hits require further assessment of null modeling and marginal associations, rather than unconditional promotion as discoveries.

Each scenario used fresh draws with deterministic seeds in `simulate_null.R` and `simulate_model.R`; all 12,000 official numerical integrations returned finite in-range p-values. Full draws and p-values, Wilson intervals, and counts are in `simulated_null.tsv`, `simulated_model_null.tsv`, and `simulation_summary.tsv`. Independent estimates of performance at ~2.2×10⁻⁸ would require vastly more appropriate simulations or a theoretically justified approach, not extrapolation from this test. Rounded source Z values and the absence of most genome-wide variants further constrain the live scan.

A subsequent [conditional-integration stress test](tail_stress/RESULTS.md) evaluates the actual stringent screening threshold under specified fixed-effect one-trait-only Gaussian models. Some scenarios have per-marker rejection substantially above the nominal threshold; this is a model-specific diagnostic, not an empirical family-wise error rate or a general judgment of PLACO+.
