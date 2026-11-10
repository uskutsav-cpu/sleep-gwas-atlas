# Insomnia–frailty correlations in two FI and physical-frailty GWAS pairs

The higher insomnia genetic correlation with FI persists when the FI and physical-frailty GWAS come from the same Mak et al. study. Using the same 1,129,146 aligned SNPs and FinnGen R13 insomnia GWAS, the estimated correlations are:

| Frailty GWAS | Insomnia genetic correlation r_g | SE |
|---|---:|---:|
| Mak FI | 0.401 | 0.028 |
| Mak modified Fried physical phenotype (FP) | 0.237 | 0.035 |
| Atkins FI meta-analysis | 0.429 | 0.039 |
| Ye Fried frailty score (FFS) | 0.244 | 0.036 |

In the Mak pair, **FP minus FI = −0.164** (paired block-jackknife SE 0.026; approximate 95% CI −0.216 to −0.112). In the original pair, **FFS minus FI = −0.185** (SE 0.040; CI −0.264 to −0.107). The Mak-minus-original difference of these two contrasts is +0.021 (SE 0.027; CI −0.031 to +0.074; unadjusted p = 0.426). This does not establish that the contrasts are equal, nor that any difference between them is absent. All comparisons are post hoc in a broader sequence of analyses; nominal p-values are not multiplicity-adjusted.

The result strengthens the observation that the FI–physical-phenotype difference is not confined to the particular Atkins/Ye GWAS pair. It does **not** establish why. [Mak et al.](https://doi.org/10.1111/acel.13868) describe their FI and modified FP summary statistics as each drawn from 403,041 UK Biobank participants, but aggregate sample counts do not prove that their individual sets are identical. The Mak and original pair reuse or may reuse UK Biobank participants and are not independent replications. Mak's FI includes sleep-related content and its physical score includes exhaustion; removing these items was not tested. Constructs, definitions, age, ascertainment, association pipelines and sample overlap may all matter. These genetic correlations are not measured clinical correlations, absolute risks, or causal effects.

## LDSC diagnostics and limits

We aligned the five input files (four frailty GWAS and FinnGen insomnia) on exactly the same SNP IDs and allele pairs; all four native LDSC runs report the same ordered SNP hash and 200 jackknife block boundaries. The insomnia h² point estimate and all delete-block values are identical across the four models. Recalculation reproduced every reported correlation SE; paired contrast SEs use differences of the corresponding delete-block estimates. This accounts for covariance in the block-jackknife calculations under LDSC assumptions, not all sources of systematic bias.

On this support, the FinnGen insomnia observed-scale h² estimate was 0.01996 (SE 0.00155), with univariate intercept 1.077. Mak FI and FP intercepts were 1.064 and 1.042, respectively. Precision is adequate for the project's h²/SE ≥4 screening rule but is not proof of unbiasedness. The LD reference is general European rather than Finnish-specific, FinnGen total N is used in the existing sleep preprocessing instead of verified per-variant N, and the original FI input uses approximated study sample size. Exact cross-source participant overlap is unknown. The data were previously restricted by HFRS and UKB-only FI availability, in addition to the current four frailty files; this is a matched-support sensitivity, not a separate discovery cohort.

[Full numerical receipt and contrasts](mak_r13_matched_contrasts.json) and [input hashes](mak_r13_matched_support.json) accompany the native exact LDSC outputs and logs (`mak_r13_old_FI`, `mak_r13_old_FFS`, `mak_r13_mak_FI`, `mak_r13_mak_FP`). [The Mak FI–FP analysis](mak2023_same_cohort_fi_fp_ldsc_review.md) provides original-source checksums, definitions and preparation. No individual genotype or phenotype data were analyzed, and no genetic risk score was evaluated.
