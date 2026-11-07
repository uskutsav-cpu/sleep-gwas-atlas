# Max-p under an extreme one-trait-only null

We simulated 1,000 new datasets per residual-correlation setting, with 44,000 markers in correlated blocks. Every marker had a strong simulated trait-1 effect (mean Z=12) but no trait-2 effect, so **every conjunction null was true**. This deliberately extreme configuration makes the conservative max-p test approach a single-trait test and exercises the boundary more strongly than the previous sparse-signal simulation.

| Correction | Residual correlation | Runs with at least one rejection | Estimated rate (95% exact interval) |
| --- | ---: | ---: | ---: |
| Bonferroni | 0 | 49/1,000 | 4.9% (3.65%–6.43%) |
| Bonferroni | 0.3 | 44/1,000 | 4.4% (3.21%–5.86%) |
| BY | 0 | 8/1,000 | 0.8% (0.35%–1.57%) |
| BY | 0.3 | 3/1,000 | 0.3% (0.06%–0.87%) |
| BH | 0 | 50/1,000 | 5.0% (3.73%–6.54%) |
| BH | 0.3 | 46/1,000 | 4.6% (3.39%–6.09%) |

All rejections are false in this simulation, so each run's false-discovery proportion is one when any marker is rejected and zero otherwise. The event rate is therefore an estimate of marker-level FDR **for this configuration**. The two correlation settings are paired using the same seed and draws; differences are not an experiment on real sample overlap. Intervals are marginal intervals for each setting and method, not simultaneous intervals or proof of equivalence with 5%. BH is descriptive: this simulation does not establish its validity under arbitrary dependence. BY was more conservative in this test.

The uniformly huge trait-1 effect and simple block correlations are not a realistic sleep–frailty model; there are no shared true effects or power estimates. With valid marginal p-values, the max-p intersection–union test is valid despite dependence between the two traits; Bonferroni and BY make different multiplicity guarantees. Real-data calibration, original GWAS confounding, ancestry and LD mismatch, phenotype construction and selected families remain important. This is not a test of the enrichment-based conjunction-FDR implementation and does not resolve the PLACO or LAVA concerns. The earlier [sparse-signal negative control](../maxp_null_calibration/RESULTS.md) found no events in 500 datasets per setting, illustrating how the null architecture changes observed conservativeness.

The [pre-result protocol](PROTOCOL.md), [code](run.R), [6,000 run summaries](runs.tsv), [aggregate table](summary.tsv), and [log](run.log) preserve the run. Sorted-threshold calculations checked the adjustment results for the first seed, and an independent Python read verified unique family keys, aggregate counts and adjusted-value ranges. These are numerical checks, not a validation of the source GWAS p-values.
