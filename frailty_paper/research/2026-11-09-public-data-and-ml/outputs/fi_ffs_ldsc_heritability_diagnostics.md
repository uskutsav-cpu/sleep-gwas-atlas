# FI–FFS LDSC: heritability and overlap diagnostics

The FI–FFS genetic correlation of approximately **0.757** is supported by relatively precise SNP-heritability estimates for both input traits in the same LDSC fit. On the 1,137,927-SNP common-support analysis:

| Diagnostic | FI | FFS |
|---|---:|---:|
| Observed-scale SNP heritability estimate | 0.1089 (SE 0.0052) | 0.0621 (SE 0.0027) |
| Estimate / SE | 20.8 | 23.1 |
| Univariate LDSC intercept | 1.0228 (SE 0.0090) | 1.0270 (SE 0.0101) |

The genetic covariance estimate was 0.06225 (SE 0.00301), yielding r_g = 0.7570 (SE 0.0212) after dividing by the square root of the product of the two heritability estimates. These numbers are LDSC estimates tied to these studies, preprocessing, sample-size assumptions, LD reference and phenotype definitions. SNP heritability is not the proportion of an individual's risk caused by genetics or an absolute prediction, and the two measures have different distributions and constructs. The table should not be read as showing that FI is intrinsically more genetic than FFS. Precision does not guarantee unbiasedness or validate the model.

The **cross-trait intercept is 0.2224** (SE 0.0073), estimated freely in this fit. FI and FFS both use UK Biobank data, so sample overlap is plausible and may contribute; correlated confounding and other model deviations can also affect this intercept. It cannot be interpreted directly as the number of overlapping participants, a phenotypic correlation, or genetic correlation. Exact shared participants are unknown. A [methodological assessment of LDSC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6758917/) explains why overlap and confounding can affect the bivariate intercept. Constraining it to zero would need a defensible no-overlap and no-correlated-confounding assumption; it was not done here.

## Stability across available SNP supports

| SNP support | SNPs | FI h² (SE) | FFS h² (SE) | FI–FFS r_g (SE) |
|---|---:|---:|---:|---:|
| FinnGen-insomnia matched | 1,138,262 | 0.1089 (0.0052) | 0.0621 (0.0027) | 0.7570 (0.0212) |
| MVP-insomnia matched | 1,151,928 | 0.1087 (0.0054) | 0.0623 (0.0027) | 0.7548 (0.0200) |
| FI–FFS–HFRS common | 1,137,927 | 0.1089 (0.0052) | 0.0621 (0.0027) | 0.7570 (0.0212) |

These analyses reuse the same FI and FFS GWAS and heavily overlapping variant sets. They are sensitivity checks of SNP support, **not independent replications**. The HFRS-restricted support contains only FI and FFS in the two-trait fit displayed here, although HFRS availability determined which SNPs were retained. The study estimates differ in participant selection, age and measurement. UKB FI includes sleep and tiredness; FFS includes exhaustion. No sleep-excluded FI was analyzed.

The [diagnostic table](fi_ffs_ldsc_heritability_diagnostics.csv) contains full precision, intercepts, covariance and source hashes. The values were extracted from the project's native exact LDSC outputs; each r_g was recomputed from the two heritabilities and covariance, matching the stored result to numerical precision. The [direct-correlation audit](fi_ffs_risk_followup.md) independently recalculated the 200-block jackknife SEs. No new GWAS or individual-level data were used. Script: `work/fi_ffs_ldsc_diagnostics.py`.
