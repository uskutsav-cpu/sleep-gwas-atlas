# SMAD5 naïve CD4 QTL colocalization with eligible credible-set indels

The [earlier Stage 4 analysis](stage4_cd4_smad5_colocalization_review.md) excluded indels. We added molecular credible-set indels with a unique exact GWAS match, reciprocal GRCh38/GRCh37 mapping and reference/normalization checks where needed, QTL MAF ≥ 1%, GWAS alternate-allele frequency between 1% and 99%, and QTL imputation r² ≥ 0.8. These additions are limited to the 11 indels in the reported CD4 SMAD5 credible set; other previously excluded variants remain excluded. Multiple source rsIDs for the same allele were deduplicated only after confirming identical QTL values. Two indels had no exact Mak representation and three had no exact FinnGen record; [exclusions](stage5_cd4_smad5_indel_exclusions.csv) are not null effects.

| Trait | Added indels | Credible-set variants before → after | Credible-set PIP sum covered before → after | H1: QTL only | H3: distinct signals | H4: shared signal |
|---|---:|---:|---:|---:|---:|---:|
| Mak FI | 9 | 53 → 62 of 70 | 0.777 → 0.797 | 0.0027 | 0.9964 | 0.0009 |
| Mak physical frailty | 9 | 53 → 62 of 70 | 0.777 → 0.797 | 0.9759 | 0.0139 | 0.0102 |
| FinnGen R13 insomnia | 8 | 51 → 59 of 70 | 0.669 → 0.693 | 0.9135 | 0.0312 | 0.0553 |

The reported credible-set PIP sum is 0.9504, so the updated input still misses part of that molecular set. These baseline posteriors use the same single-signal approximate-Bayes-factor model and priors as Stage 4 (QTL SD 0.15; p1 = p2 = 10⁻⁴; p12 = 10⁻⁵; original trait-specific scales). H0 and H2 account for the remainder. Adding these indels leaves the earlier qualitative pattern: Mak FI favors distinct signals, while Mak physical frailty and FinnGen insomnia are dominated by QTL-only evidence. These are conditional results under the model, not evidence that sharing is impossible or a conclusion about all FI/FFS studies. The full [prior sensitivity table](stage5_cd4_smad5_indel_coloc_sensitivity.csv), [coverage](stage5_cd4_smad5_indel_coloc_coverage.csv), and source-specific updated aligned CSVs are included.

Before additions, the original baseline was recomputed successfully. The new baseline calculations agreed with `coloc` v5.2.3 formula calculations to less than 2 × 10⁻¹⁶ ([formula check](stage5_cd4_smad5_indel_coloc_formula_check.csv)). This is an arithmetic check, not independent validation of data preprocessing, imputation, prior choices or the single-signal assumption. Missing molecular-set variants, unexamined noncredible-set indels and palindromic SNPs, sample overlap, ancestry, phenotype scaling and multiple independent signals limit the interpretation. The effect allele and source build for Mak remain inherited assumptions. Molecular QTLs and chromatin assays may overlap AlphaGenome training; no mediation or disease causality is established.

Full-precision AlphaGenome expression/splice/chromatin predictions and allele- and cell-matched independent functional validation are still outstanding.
