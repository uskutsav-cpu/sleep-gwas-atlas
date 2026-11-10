# FI, FFS and HFRS: same-SNP genetic correlation matrix

Using the same 1,137,927 aligned SNPs and European LD reference in all three pairwise LDSC runs, the exploratory genetic correlations are:

| | FI | FFS | HFRS |
|---|---:|---:|---:|
| **FI** | 1 | 0.757 | 0.640 |
| **FFS** | 0.757 | 1 | 0.550 |
| **HFRS** | 0.640 | 0.550 | 1 |

The FI–FFS standard error is 0.0212 (approximate 95% CI 0.715–0.799); it closely matches the independently located earlier FI–FFS result on a slightly different support. This is genome-wide genetic overlap between these source-defined traits, not a correlation of observed clinical scores, percentage of shared variants, or individual risk prediction.

I verified all three runs used the same ordered SNP hash and 200 block boundaries. Where the same trait occurs in two models, its heritability and all 200 deletion estimates match exactly. The three point-estimate and 200 delete-block genetic correlation matrices are positive definite (smallest eigenvalue across them 0.224); the joint jackknife sampling covariance is positive definite (smallest eigenvalue 0.000340). These numerical checks support consistency of the calculation but do not establish model validity. The matrix and covariance are in `frailty_three_way_common_support_ldsc.json`, with exact native results in the three companion `.exact.json` files.

The FFS–HFRS minus FI–HFRS paired contrast is −0.090 (SE 0.0265; approximate 95% CI −0.142 to −0.038; unadjusted p=0.00065). It is post hoc and exploratory; details are in `hfrs_fi_ffs_paired_contrast.md`. The original common FI/FFS support was restricted by availability in a FinnGen-insomnia input and was further intersected with HFRS; this choice does not itself estimate a sleep association.

The HFRS LDSC intercept remains elevated (~1.188), and the HFRS GWAS is Finnish while LD reference scores are general European. Per-SNP sample sizes were unavailable and total N was assumed in preprocessing. These concerns remain for the HFRS correlations and their contrast even though numerical checks pass. FI and FFS share UK Biobank participants; their HFRS cohort is different by description but exact intersections are not verified. Definitions, ages and ascertainment differ across the GWAS. The result does not establish why the correlations differ, demonstrate causal variants or replace independent outcome validation.
