# Checking 22 candidates in another FI and physical-frailty GWAS pair

All 22 variants in the exploratory insomnia–hospital-frailty candidate set have unique matching records in [Mak et al.'s UK Biobank FI and modified Fried phenotype GWAS](https://doi.org/10.6084/m9.figshare.24591573.v1), by rsID, assumed GRCh37 chromosome/position and allele pair. After aligning to the same counted allele, **all 22 FI point estimates have the same sign as the original Atkins FI estimates, and all 22 physical-phenotype point estimates have the same sign as the original Ye FFS estimates**. All 22 have Mak FI p < 0.05 and 11 have both Mak FI and physical-phenotype p < 0.05.

This is a **within-UK-Biobank consistency check**, not independent replication: Mak, Atkins and Ye use UK Biobank participants, and the candidates were selected using related prior information. The overlap in people and definitions can contribute to the strong agreement. All p-values are unadjusted, the analyses were explored after selection, and some variants are correlated. Matching signs for imprecise estimates does not establish that they affect both traits. A non-significant result is not evidence of no effect.

## The six candidates with nominal evidence across the earlier six sources

| Variant | Positional label | Mak FI p | Mak physical-phenotype p | Original Ye FFS p |
|---|---|---:|---:|---:|
| rs12536335 | FOXP2 | 1.21×10⁻¹⁴ | 0.00279 | 0.0017 |
| rs817509 | RP11-384F7.1 | 3.36×10⁻⁶ | 0.00534 | 0.0099 |
| rs2838788 | ADARB1 | 2.62×10⁻⁵ | 0.00813 | 0.0040 |
| rs34290943 | LAMB2 | 1.10×10⁻⁸ | 0.000198 | 0.00011 |
| rs12671330 | FOXP2 | 1.53×10⁻¹⁴ | 0.00152 | 0.0017 |
| rs7734045 | AC009014.3 | 0.00134 | 0.0517 | 0.015 |

For rs7734045 the Mak physical-phenotype p-value is just above 0.05 while the prior FFS p-value is below it. This difference in threshold status does **not** demonstrate different underlying effects; comparing effects requires a suitably specified model with their covariance, and the phenotypes are not identical. Both point-estimate signs agree after alignment. The two FOXP2-labeled variants lie in the same region and have nonzero reference-panel LD, so should not be counted as independent loci without further analysis. Positional labels do not identify causal genes. Source-specific multiallelic coding and the assumed build still need confirmation. Original FI cohort estimates for rs12536335 differ in sign under the inferred TwinGene adapter; that is an additional reason for caution despite its strong UKB-based associations.

## Traceability

We scanned all 8,983,921 paired source rows in the previously checksum-verified FI and FP files and aligned the released fastGWA A1-effect betas to each candidate's counted allele. The [lookup table](mak2023_22_candidate_lookup.csv) retains effects, standard errors, Wald intervals, p-values, sample sizes and source alleles, along with the original FI/FFS estimates. FI and FP effect units differ, so coefficients must not be added or directly compared by magnitude. Source documentation, checksums and phenotype descriptions are in [the Mak LDSC review](mak2023_same_cohort_fi_fp_ldsc_review.md); [counts and output hash](mak2023_22_candidate_lookup.json) and `work/mak2023_22_candidate_lookup.py` document the extraction. No individual outcomes or predictive model were assessed.
