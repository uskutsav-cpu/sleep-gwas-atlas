# UK Biobank FI sensitivity: correlation with hospital frailty

The UK Biobank-only FI–HFRS estimate is very similar to the FI meta-analysis–HFRS estimate. On the same 1,134,953 variants, the genetic correlations (LDSC rg ± standard error) are:

| Frailty measure | Correlation with Finnish HFRS | SE |
|---|---:|---:|
| Atkins UK Biobank-only FI | 0.6421 | 0.0244 |
| Atkins FI meta-analysis (UK Biobank and TwinGene) | 0.6399 | 0.0257 |
| Fried Frailty Score (FFS) | 0.5492 | 0.0230 |

The paired difference, UK Biobank-only FI minus meta-analysis FI, is **0.0022** (SE 0.0054; 95% CI −0.0083 to 0.0128; two-sided p = 0.677). FFS minus UK Biobank-only FI is **−0.0929** (SE 0.0262; 95% CI −0.1441 to −0.0416; p = 0.000385). The latter is a post hoc, unadjusted comparison. Both comparisons use a paired 200-block jackknife, accounting for covariance between estimates. This sensitivity suggests the Swedish TwinGene contribution does not explain the observed difference between FI–HFRS and FFS–HFRS correlations; it does not establish the cause of that difference.

## Method and checks

We intersected the previously harmonized FI, FFS, and Finnish HFRS set with Atkins cohort-specific UKB1 summary statistics. UKB1 Z scores were calculated as BETA_UKB1 / SE_UKB1 with N = 164,610; per-variant sample sizes were not available or verified. The source effect allele was aligned to the existing common set after checking rsID, GRCh37 chromosome and position, and alleles; ambiguous duplicated source rsIDs were excluded. All three native LDSC runs used the same ordered SNPs, EUR LD reference and weights, and 200 blocks. The ordered SNP hashes, block boundaries, and HFRS heritability and delete-block estimates matched across runs; the reported genetic correlation standard errors were independently reproduced from the delete-block values. Exact fit outputs and the input-hash receipt accompany this note.

## Interpretation and limitations

FI includes sleep/tiredness-related deficits and FFS includes exhaustion; this analysis does not remove that content overlap. HFRS is a hospital-record frailty proxy, not a clean sleep-free frailty measure. Sample overlap and phenotype definitions may differ between datasets. HFRS has an elevated LDSC intercept (1.188, SE 0.014), and a general European LD reference may not optimally model the Finnish sample. Constant cohort N approximates variable per-SNP N. These are exploratory population-level genetic correlations, not causal effects, clinical predictions, evidence of association for individual candidate variants, or external validation of a risk score. UK Biobank-only FI is part of the meta-analysis, so the two estimates are not independent; the paired analysis accounts for their shared block-level variation.

Source: [Atkins et al. cohort-specific GWAS files](https://figshare.com/articles/dataset/Genome-Wide_Association_Study_of_the_Frailty_Index_-_Atkins_et_al_2021/9204998). Detailed calculations are in `ukb_fi_hfrs_sensitivity.json`; the three `*_HFRS_ukbcommon.exact.json` files contain the native fit and block-level estimates.
