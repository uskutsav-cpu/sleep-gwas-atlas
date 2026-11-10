# How closely do direct-array partial sums track the full prototypes?

In the existing 503-person European 1000 Genomes reference, I computed each research score twice: first using all score variants, and then using only the variants with exact rsID-and-position matches in the public Human660W-Quad v1_H manifest. The result is the Pearson correlation between the two genetic sums in those same reference individuals.

| Prototype | Weight set | Direct manifest markers | Correlation with full sum |
|---|---|---:|---:|
| 15-marker | FI | 4/15 | 0.468 |
| 15-marker | FFS, inferred direction | 4/15 | 0.634 |
| 15-marker | Insomnia | 4/15 | 0.488 |
| 22-marker weighted | FI | 5/22 | 0.570 |
| 22-marker weighted | FFS, inferred direction | 5/22 | 0.424 |
| 22-marker weighted | Insomnia | 5/22 | 0.474 |
| 22-marker equal count | Shared-direction count | 5/22 | 0.543 |

These modest reference correlations suggest that a direct-marker-only sum could differ substantially from the full prototype. They are **not** estimates of phenotype variance explained, risk prediction, calibration, or performance in ACT. The ACT release's actual array revision, probe alleles, strand, post-QC availability, and individual genotypes are not established. This calculation assumes all name-and-position matches have the correct counted-allele contrast, an optimistic assumption. It does not include imputed variants or the nearby tags identified in the separate [tag screen](stage6_act_array_reference_tag_review.md); ACT sequencing and imputed data could change marker coverage.

Reference dosages were decoded from the existing MAGMA EUR Phase 3 PLINK panel after exact GRCh37 position and unordered allele-pair matching. All 503 reference individuals had complete genotypes for each scored set. Sums use the original coefficients and counted alleles, including the distinct allele conventions of the count and weighted 22-marker files. The FFS weights' direction is inferred, and the prototypes are postselected and unvalidated. Pooled reference correlations may reflect LD, population structure, relatedness or sampling noise, and should not be extrapolated across populations without testing. No participant outcomes were analyzed. Squared correlation in the accompanying output is the fraction of *variation in the full genetic sum linearly associated with the partial sum in this reference sample*, not variance in insomnia or frailty.

Reproduce with `python3 work/act_array_partial_score_reference.py`. [Full table](stage6_act_direct_partial_reference_scores.csv) and [method receipt](stage6_act_direct_partial_reference_scores.json) provide the unrounded values.
