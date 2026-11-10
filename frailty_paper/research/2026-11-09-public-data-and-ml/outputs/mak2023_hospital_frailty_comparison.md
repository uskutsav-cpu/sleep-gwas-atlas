# FI and physical frailty versus hospital frailty

In the additional Mak et al. frailty GWAS pair, Finnish hospital frailty risk score (HFRS) has a higher estimated genetic correlation with FI than with the modified Fried physical phenotype (FP). On the same 1,129,146 SNPs, the result also holds for the original Atkins FI and Ye FFS pair:

| Frailty GWAS | HFRS genetic correlation r_g (SE) |
|---|---:|
| Mak FI | 0.644 (0.018) |
| Mak modified Fried physical phenotype | 0.537 (0.022) |
| Atkins FI meta-analysis | 0.642 (0.026) |
| Ye Fried frailty score | 0.547 (0.023) |

For the Mak pair, FP minus FI is **−0.107** (paired SE 0.018; approximate 95% CI **−0.143 to −0.072**). For the original pair, FFS minus FI is **−0.094** (SE 0.027; CI −0.147 to −0.042). The difference between those two contrasts is −0.013 (SE 0.017; CI −0.047 to +0.021, unadjusted p = 0.452); this does not establish equality. All analyses are post hoc and not adjusted for the broader sequence of comparisons.

This pattern provides context for the project's insomnia result: a higher FI correlation also appears for another broad frailty-related phenotype. It does not show that the insomnia difference is nonspecific, caused by FI's items, or mediated by a particular disease. FI and HFRS both capture broad health deficits or clinical diagnoses, whereas the physical score emphasizes functioning. HFRS includes broad diagnostic codes that may include drowsiness or exhaustion, so it is not a verified sleep-free proxy. The sources' UK Biobank samples overlap or may overlap, and this is not an independent replication. Correlations are genome-wide genetic estimates, not individual predictions, correlations of measured scores, or causal effects.

## Major limitations and checks

The HFRS univariate LDSC intercept is **1.187** on this support, and the analysis uses a general European LD reference for a Finnish GWAS. HFRS per-variant sample sizes were unavailable in the existing preprocessing, which used total N = 500,737. These issues may affect all four estimates; an h² precision check and matched numerical outputs do not rule out bias. A prior attempt to obtain an appropriate Finnish LD reference was unsuccessful. See [the HFRS intercept audit](hfrs_insomnia_ldsc_intercept_audit.md). The exact overlap between cohorts, source phenotype coding and selection remain uncertain.

The four fits used the same matched SNP and allele order and the same 200 block boundaries. HFRS h² and its delete-block values were identical across the fits; each correlation SE and the paired contrasts were recomputed from the native deletion estimates. The support was already restricted by prior FinnGen-insomnia and HFRS availability as well as the four frailty files. Same block jackknife covariance helps compare fitted correlations under LDSC assumptions; it cannot identify the cause of the differences or correct systematic bias.

[Input receipt](mak_hfrs_matched_support.json), [exact contrasts and diagnostics](mak_hfrs_matched_contrasts.json) and the four `mak_hfrs_*.exact.json` files and logs preserve the analysis. The [Mak FI–FP review](mak2023_same_cohort_fi_fp_ldsc_review.md) records the source statistics, checksum verification and phenotype definitions. Scripts: `work/mak2023_hfrs_common.py` and `work/mak_hfrs_contrasts.py`. No individual participant data or score prediction was evaluated.
