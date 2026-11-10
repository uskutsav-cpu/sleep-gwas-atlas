# UK Biobank-only FI sensitivity against FinnGen insomnia

On an identical set of 1,134,953 variants, the UK Biobank-only FI and FI meta-analysis have very similar genetic correlations with FinnGen clinical insomnia. Both estimates are higher than the FFS estimate in this exploratory comparison.

| Frailty GWAS | Insomnia genetic correlation (rg) | SE |
|---|---:|---:|
| Atkins UK Biobank-only FI | 0.419 | 0.037 |
| Atkins FI meta-analysis, UK Biobank plus TwinGene | 0.426 | 0.039 |
| Ye Fried Frailty Score (FFS) | 0.244 | 0.037 |

The paired difference **UKB-only FI minus meta-analysis FI** is −0.0074 (SE 0.0072; 95% CI −0.0215 to +0.0066; two-sided p=0.300). The paired difference **FFS minus UKB-only FI** is −0.1744 (SE 0.0385; 95% CI −0.2497 to −0.0990; unadjusted p=5.77×10⁻⁶). Confidence intervals use the normal approximation. Both contrasts are post hoc, and the second p-value is unadjusted for the broader set of analyses.

This sensitivity suggests that including the Swedish TwinGene cohort is not the main reason that the original FI and FFS studies differ in their observed genetic correlations with this specific insomnia GWAS. It does **not** identify the reason for the difference. The UKB FI still contains sleeplessness and tiredness-related items and the FFS includes exhaustion; FI and FFS also differ in construct, participant selection and measurement. UKB-only FI is part of the meta-analysis, so the two FI estimates are correlated; the paired calculation explicitly uses their shared delete-block variation. This is not independent replication of FI or FFS, and genetic correlation is not an individual risk score or causal effect.

## Methods and diagnostics

The project’s previously harmonized, common-support FinnGen R13 clinical-insomnia file was intersected with the existing UKB-only FI/HFRS sensitivity set. The latter already required FI, FFS, HFRS and UKB-only FI availability; therefore this is an HFRS-restricted subset of the earlier R13 analysis, not an independently selected replication. Matched variants retained identical SNP and A1/A2 order. LDSC used the same European reference and weights and 200 blocks in all three runs. The ordered SNP hash and all block boundaries match, and the insomnia h² estimate and all delete-block estimates are identical across fits. The three reported correlation SEs were independently reproduced from their delete-block values; paired SEs were computed from differences of the same delete blocks. The FI, FFS and insomnia h²/SE ratios exceed four, the project's precision gate, but this does not prove unbiased correlations.

FinnGen insomnia's univariate intercept was 1.077 (SE 0.009). The reference is a general European LD panel rather than a Finnish-specific panel; phenotype, ancestry, imputation and confounding differences can affect results. The R13 total N is used in the existing input rather than per-variant N, and the UKB-only FI uses N=164,610 as an approximation because variant-specific sample sizes were not verified. Exact participant overlap between source cohorts has not been established. Identical sample sets within a cohort should not be assumed.

`ukb_fi_finngen_insomnia_sensitivity.json` records point estimates, uncertainty, SNP hash and diagnostics. The three `*_R13_ukbcommon.exact.json` files retain native LDSC output and delete-block results. The prior common-support receipt and UKB FI source handling are documented in `ukb_fi_hfrs_sensitivity.md` and `hfrs_fi_ffs_paired_contrast.md`. No individual-level data, prediction or new GWAS was generated.
