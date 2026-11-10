# autoMR frailty example data: score-validation feasibility

Checked 2026-10-28. The freely downloadable [autoMR 1.2.1 CRAN source package](https://cran.r-project.org/src/contrib/autoMR_1.2.1.tar.gz) includes files named `fi_49item`, `fried_frailty`, and `merged_data`. The [package manual](https://cran.r-universe.dev/autoMR/doc/manual.html) gives the first two phenotype-sounding descriptions, but inspection of the actual R objects shows that they are SNP-level Mendelian randomization example data, not individual participant genotypes or frailty measurements.

- `fi_49item`: 10 rows and 13 columns, one row per instrument, with exposure and outcome alleles, beta coefficients, standard errors and frequencies. The exposure label is `Zn` and the outcome label is `fi_49item`.
- `fried_frailty`: 10 rows and 13 analogous columns, outcome label `fried_frailty`.
- `merged_data`: 20 rows and 13 columns. Its row count is consistent with the two examples combined; it is not 20 participants.
- Neither the 15-marker candidate list nor the 22-marker candidate list has an exact rsID match to the 10 instruments in `fi_49item` (comparison with `integrated_15_candidate_review.csv` and `22_candidate_risk_direction_table.csv`). Exact rsID comparison does not assess LD proxies.

These examples cannot calculate or validate our scores: they do not contain participant-level genotypes linked to measured FI or Fried frailty, and they do not cover our candidate variants. The package documentation inspected does not establish the underlying cohort, independence from the score-development data, genome build, or phenotype provenance; the example coefficients should not be treated as independent replication. No inferential comparison of FI and Fried frailty was performed.

Reproducibility: source archive SHA-256 `927b82e9ac5215746e4b8d1e0ff3321d9c360fc7f941a0a0a800b6bd73ec9c03`. Loaded the `.rda` files in base R, inspected dimensions and column names, and intersected `Instrument` against the `variant` column in each existing candidate CSV. The archive is cached at `work/automr/autoMR.tar.gz`.
