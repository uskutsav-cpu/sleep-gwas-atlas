# Low grip-strength GWAS: candidate comparison

[Jones et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7844411/) meta-analyzed low grip strength in 256,523 older adults from 22 cohorts. Low grip is one physical component of frailty, not the full Fried frailty score or the deficit-accumulation frailty index. The study includes UK Biobank, which also contributes to the candidate-discovery sources; this is **not an independent replication** of the project scores.

The original [GWAS Catalog file GCST90007526](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90007001-GCST90008000/GCST90007526/) was fully streamed and its MD5 matched the publisher's checksum. All 15 markers in the weighted prototype and all 22 markers in the separate count had one exact GRCh37 coordinate and allele-pair match. Four of 15 and five of 22 had nominal p < 0.05; these overlapping, postselected counts are descriptive, not multiple-testing-adjusted evidence. Matched effects and heterogeneity statistics are in `stage6_grip_gwas_candidate_aligned.csv`.

A postselected, summary-statistic weighted-Z comparison of the 15-marker prototype, using 1000 Genomes EUR503 reference LD, yielded the following approximate results:

| Prototype weight | Approximate Z | Approximate two-sided p | Leave-one-marker-out Z range |
|---|---:|---:|---:|
| FI | 4.04 | 5.4 × 10⁻⁵ | 2.97 to 4.50 |
| FFS (inferred orientation) | 4.54 | 5.7 × 10⁻⁶ | 3.49 to 4.79 |
| Insomnia | 3.23 | 0.00124 | 2.34 to 3.79 |

These are weighted combinations of marginal SNP Z statistics, not estimated score coefficients or predictive performance. They assume the target GWAS Z-score correlation is approximated by reference LD, Hardy–Weinberg genotype scaling from the target allele frequency, and sufficiently comparable contributing samples and errors across variants. The meta-analysis contribution pattern is not identical across variants, and reference-panel LD may not match target LD, so the p-values are approximate and their calibration is unverified. Three weights and leave-one-out analyses are postselected; no prospective multiplicity correction is claimed. The FFS weight orientation has additional unresolved source issues noted in the main synthesis.

These signals are consistent with some relationship between the selected variants and grip weakness under the stated assumptions, but do not establish a shared causal mechanism, individual score calibration or discrimination, or predictive accuracy for FI or FFS. They cannot satisfy independent validation because of study overlap. The separate 22-marker count was later assessed in a [full-marker summary sensitivity](stage6_grip_shared22_weighted_z_review.md), subject to the same overlap limitation.

Reproducibility: `stage6_grip_gwas_extract_receipt.json` records the verified full-file MD5 and row counts; `stage6_grip_gwas_candidate_align_receipt.json` records the matches. `stage6_grip_sparse_score_summary.py`, its CSV and receipt contain the aggregate calculation and leave-one-out analyses. Source allele frequencies and effects refer to the published meta-analysis and require the caveats above.

## Influence and cohort-coverage audit

Of the 15 markers, ten source `Direction` strings have no question marks, four have one and one has two. Under [METAL's documented convention](https://csgwiki.sph.umich.edu/wiki/METAL_Quick_Start), question marks indicate that a contributing study lacks a marker result. The strings do not provide the per-marker sample sizes or the pairwise covariance needed to verify the aggregate Z test. This observed difference in cohort participation reinforces the approximate nature of its p-values.

For the FI and insomnia weights, omitting rs8067974 produces the largest absolute change in the aggregate Z (4.04 to 2.97, and 3.23 to 2.34 respectively). For the inferred FFS weights, omitting rs323509 has the largest change (4.54 to 3.49). These exploratory leave-one-out checks describe sensitivity to individual markers; they do not establish replication or causality. Reproduction: `stage6_grip_influence_audit.py` and `stage6_grip_influence_audit.json`.
