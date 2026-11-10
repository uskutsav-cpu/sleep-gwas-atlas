# HFRS–insomnia LDSC: source and intercept audit

I inspected the project's saved native LDSC results for three insomnia inputs against the same Finnish hospital frailty risk score (HFRS) GWAS and independently recomputed each SE from its 200 stored delete-block estimates. All three recalculated SEs agree with the stored values to numerical precision. The complete values and source files are in `hfrs_insomnia_ldsc_diagnostics.csv`; this verifies arithmetic for the stored runs, not the underlying model assumptions.

| Insomnia source | HFRS rg (approx. 95% CI) | HFRS single-trait intercept in that paired run | Cross-trait intercept |
|---|---:|---:|---:|
| UK Biobank | 0.439 (0.388–0.490) | 1.190 (SE 0.014) | −0.0004 (SE 0.0060) |
| FinnGen R13 | 0.626 (0.563–0.689) | 1.189 (SE 0.014) | 0.1767 (SE 0.0070) |
| MVP | 0.566 (0.502–0.630) | 1.190 (SE 0.014) | 0.0028 (SE 0.0071) |

Each estimate uses a different insomnia source, definition and SNP support, while reusing the **same HFRS GWAS**. The differences have not been tested with a paired covariance-aware comparison and should not be ranked by their intervals alone. R13 insomnia and HFRS R12 overlap in FinnGen participants; its sizeable cross-trait intercept is consistent with overlap and/or correlated biases, but does not quantify the exact overlap. The other pairs are described as different source cohorts, not verified as having zero shared participants. A cross-trait intercept near zero does not prove absence of all bias.

The standalone HFRS LDSC run reported h²=0.0536 (SE 0.0024), intercept 1.1908 (SE 0.0121), mean chi-square 1.7137 and an LDSC ratio of 0.2673 (SE 0.017). The intercept is well above one and appears similarly in all three paired runs. This deserves statistical review, but it is **not** a direct measurement that exactly 26.7% of the result is confounding. Polygenic architecture, residual structure, LD-reference mismatch and other model violations may contribute; the present results cannot identify their proportions.

The stored pipeline used European-reference LD scores for Finnish HFRS, rather than a Finnish-matched reference. Its HFRS preprocessing assumed a constant N=500,737 because per-variant N was unavailable and did not invent an imputation INFO field. Those are material limitations. The prior project records an unsuccessful anonymous attempt to obtain the intended Finnish-matched LD reference (HTTP 403); no bypass was attempted. A valid follow-up requires an authorized, suitable Finnish LD reference and an explicitly designed comparison, along with review of sample-size assumptions and the reference's applicability to each cross-ancestry pair. Constraining the intercept to one would not by itself resolve this issue.

All three correlations are exploratory genome-wide associations with HFRS, a distinct diagnosis-based phenotype that includes categories encompassing drowsiness or exhaustion. They are not proof of causal mediation, an FI-versus-FFS correlation, an individual prediction, or a local shared variant. The project's failed FI×sleep local-correlation QC gates remain failed.
