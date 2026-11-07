# HFRS–insomnia covariance-intercept stress test

The original unconstrained same-support analysis estimates cross-trait covariance intercept 0.1790 (SE 0.0072) and r_g=0.6171 (SE 0.0326). To show how a fixed-intercept assumption changes the calculation, we re-ran the same pair while leaving both single-trait intercepts free:

| Cross-trait intercept assumption | Resulting r_g (SE) | Interpretation |
|---|---:|---|
| Estimated from data | 0.6171 (0.0326) | Original exploratory estimate. |
| Fixed at 0.165 | 0.6709 (0.0331) | Artificial fixed-value sensitivity. |
| Fixed at 0.193 | 0.5637 (0.0315) | Artificial fixed-value sensitivity. |
| Fixed at 0 | Out of bounds; LDSC reports `nan` | Deliberate stress case; zero is not justified by the available overlap information. |

The values 0.165 and 0.193 are approximately 0.179 ± 0.014, chosen after seeing the main intercept. They do **not** form a confidence interval for r_g or a validated uncertainty range: the intercept and genetic-covariance slope are jointly estimated and correlated, and fixing the intercept changes the model. The SEs in the fixed runs condition on those values. The zero setting is particularly unrealistic as a default for two FinnGen-based outcomes; the out-of-bounds result demonstrates sensitivity to that assumption, not evidence that the true genetic correlation is above one. The intercept may reflect sample overlap, correlated phenotypic errors and model/reference mismatch; its value alone cannot separate these causes. Per-variant N and true participant overlap remain unknown.

The [FinnGen LDSC workflow](https://github.com/FINNGEN/LDSC/blob/master/README.md) documents prebuilt Finnish and European LD-score options, but its listed default Finnish resource is in a Google Cloud location that returned a VPC service-controls HTTP 403 on 27 October 2026. No protected access was bypassed. We therefore cannot run that reference as an ancestry sensitivity from this environment. The original [LDSC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4495769/) discusses how population LD-score mismatch can bias intercepts; it does not establish that mismatch explains this dataset's high HFRS intercept. The 99-person Finnish subset of 1000 Genomes is not substituted as an adequately powered genome-wide reference.

The [pre-result protocol](PROTOCOL.md), `receipt.json`, run script and complete logs preserve the assumptions and output. The unconstrained main estimate remains the project result; this stress test does not validate PLACO, local correlations, causality or independent replication.
