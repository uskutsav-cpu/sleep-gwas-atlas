# How missing a marker changes the 22-marker weighted research sums

The current scoring tool requires all 22 variants. To help evaluate a candidate validation cohort's marker coverage, we asked how the full score and a score omitting one marker would correlate in the 503 European-labeled 1000 Genomes reference participants used in the existing variance analysis. This is a **genotype-only technical sensitivity**, not a disease association or validation. Full and reduced scores share 21 markers, so high correlations are expected.

| Weight set | Single omission with lowest correlation | Correlation with full sum | Reduced/full score SD |
|---|---|---:|---:|
| FI | rs12671330 | 0.951 | 0.924 |
| FFS | rs28528997 | 0.860 | 0.837 |
| Insomnia | rs1909650 | 0.961 | 0.938 |

For example, if rs28528997 were unavailable, the reduced FFS-weighted sum would be noticeably different in this reference panel. This does **not** establish that the marker predicts frailty, that a reduced score performs worse in patients, or that any particular target cohort lacks the marker. Correlations and SD ratios depend on allele frequency and covariance, ancestry, and sampling; this pooled panel includes multiple subpopulations and may contain relatives. Other missingness patterns can have different effects. The FFS effect orientation remains inferred and the variants' exact target allele representations need verification.

The calculation reuses the saved [22-marker reference variance decomposition](reference_22_score_variance_review.md), without retrieving individual genotypes. If `a = Cov(wᵢXᵢ, S)/Var(S)` and `d = Var(wᵢXᵢ)/Var(S)`, the reduced/full correlation is `(1−a)/sqrt(1−2a+d)` and the SD ratio is `sqrt(1−2a+d)`. The script checks this identity using synthetic data. [All 66 results](stage6_weighted22_missing_marker_sensitivity.csv), [input/output hashes](stage6_weighted22_missing_marker_sensitivity.json), and [code](../work/weighted22_missing_marker_sensitivity.py) are available. The synthetic check verifies the formula, not the original genotype extraction or source-allele coding.

For a real target cohort, check the build, alleles, call/imputation quality and coverage before scoring, and prespecify any missing-variant policy before examining outcomes. Do not substitute these reference correlations for actual target-cohort validation or silently treat an incomplete sum as the original score.
