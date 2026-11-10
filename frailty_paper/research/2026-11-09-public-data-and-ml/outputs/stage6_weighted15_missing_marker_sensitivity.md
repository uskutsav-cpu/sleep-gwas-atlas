# Missing-marker sensitivity for the 15-marker research scores

We used the existing 1000 Genomes Phase 3 European-labeled reference covariance decomposition to quantify how the complete 15-marker research sums change if one marker is omitted. For each weight set we removed each marker in turn and calculated the correlation between the resulting 14-marker sum and the original sum, and the ratio of their standard deviations. These are **reference-genotype diagnostics**, not validation of score associations or predictive performance.

| Weight set | Omission yielding lowest reduced/full correlation | Correlation | Reduced/full SD |
|---|---|---:|---:|
| FI | rs8067974 | 0.927 | 0.967 |
| FFS | rs323509 | 0.879 | 0.857 |
| Insomnia | rs113851554 | 0.911 | 0.878 |

This means, for example, that the score calculated without rs323509 differs materially from the full FFS-weighted sum in this reference panel. It does **not** establish that the marker is predictive or that a target cohort lacks it. The full and reduced sums share 14 markers, so correlation is expected; different missingness patterns, allele frequencies, LD, ancestry and sampling can change the result. Pooled samples may be related and reflect population structure. FFS effect orientation remains inferred, and some markers have multiallelic representations needing source- and target-specific checks. Neither high correlation nor available imputation establishes valid prediction.

The input contains 503 reference individuals' aggregate dosage-covariance decomposition, not measured patient outcomes. For score S and a marker contribution Xᵢ, define `a = Cov(Xᵢ,S)/Var(S)` and `d = Var(Xᵢ)/Var(S)`. Then the reduced/full correlation is `(1−a)/sqrt(1−2a+d)` and the SD ratio is `sqrt(1−2a+d)`. These identities include the source sample's covariance. The script checks the correlation formula with synthetic data; that tests arithmetic, not the original genotypes, allele coding, biological validity, or sensitivity to multiple missing markers.

[All 45 omission results](stage6_weighted15_missing_marker_sensitivity.csv), the [checksum receipt](stage6_weighted15_missing_marker_sensitivity.json), [source decomposition](reference_15_score_variance_decomposition.md), and [script](../work/weighted15_missing_marker_sensitivity.py) support reproducibility. Run `work/stage4_env/bin/python work/weighted15_missing_marker_sensitivity.py` from the project root. A real validation study should freeze the score and missing-marker rule before outcomes are examined, check target build and allele contrast and genotype quality, and report completeness; an incomplete sum must not be silently labeled the full score.
