# rs323509-region testis splicing and trait colocalization

For the testis splice junction `5:104739412:104739497:clu_21568_-` in `ENSG00000251574`, an exploratory single-signal model gives **H4 = 0.663 for FinnGen insomnia** at the baseline priors below. This is a lead for further study, not a confirmed shared variant or a disease mechanism in testis. The same baseline model favors different molecular and GWAS signals for the three frailty traits. No three-trait colocalization was performed.

| GWAS | Matched variants | Molecular credible-set variants represented | H1: QTL only | H3: distinct signals | H4: shared signal |
|---|---:|---:|---:|---:|---:|
| Mak FI | 2,286 | 4/4 | <0.001 | 0.929 | 0.071 |
| Mak modified physical frailty | 2,286 | 4/4 | <0.001 | 0.997 | 0.003 |
| Ye FFS | 2,213 | 4/4 | <0.001 | 0.913 | 0.087 |
| FinnGen R13 insomnia | 2,257 | 4/4 | 0.158 | 0.179 | 0.663 |

H0/H2 are negligible. Posterior probabilities are conditional on a single causal variant per trait, the included region and specified priors. The default illustration sets p1 = p2 = 10⁻⁴, p12 = 10⁻⁵, QTL effect prior SD = 0.15 on an assumed normalized splice scale, and trait prior SD = 1.1004 for Mak FI, 0.12325 for Mak physical frailty, 0.12082 for Ye FFS, or 0.20 log-odds for insomnia. The Mak scales are working estimates from another region; the Ye scale was estimated from SE, total sample size 386,565, and QTL MAF as a proxy for unavailable GWAS allele frequencies. With GWAS prior SD multiplied by 0.5–2 and p12 varied from 10⁻⁶ to 5×10⁻⁵, **insomnia H4 ranges 0.119–0.928**; FI H4 ranges 0.0076–0.277 and Ye FFS H4 0.0095–0.324. [Complete sensitivity results](stage3_rs323509_testis_coloc_sensitivity.csv).

All four reported molecular credible-set variants are represented in all four matched inputs, but this does not prove that all trait variants are present or that each trait has a single signal. The QTL was fine-mapped in a 2-Mb source window; this comparison covers only chr5:104246478–105246478, a 1-Mb candidate-centered window. The p12 prior and scaling assumptions are influential. rs323509 and the other molecular credible-set variants are correlated; the analysis cannot reliably distinguish them and should not be interpreted as assigning a causal allele. High H4 for one insomnia GWAS and lower H4 for FI/FFS do not establish shared insomnia–frailty genetics. Testis is not established as a disease-relevant tissue; other tissues may share regulatory effects.

## Inputs and validation

QTL: eQTL Catalogue r7 QTD000340, GTEx v8 testis leafcutter splicing, 322 samples per metadata. The `cc` source contains selected molecular events; it is not a comprehensive negative screen. [Regional receipt](stage3_rs323509_testis_extract_receipt.json), [molecular fine mapping](stage3_rs323509_testis_finemapping.md), [harmonization receipt](stage3_rs323509_testis_harmonization_receipt.json), and [coverage](stage3_rs323509_testis_coloc_coverage.json). The source sample overlap and population matching with GWAS require further verification. Mak FI and FP share a cohort; Ye FFS also uses UK Biobank, so they are not independent studies.

Matching required a unique rsID, coordinate and allele pair, positive-strand unique GRCh37-to-GRCh38 liftover for the UK Biobank GWAS, non-palindromic biallelic SNPs and QTL MAF ≥1%; Mak/FinnGen GWAS MAF ≥1% was also required. Ye's source filename indicates MAF ≥1%, but variant-specific GWAS frequencies were unavailable to check. The Catalogue QTL imputation r² is missing for the candidate and was not filtered; GWAS imputation quality was not independently checked. Allele coding and the multiallelic status of rs323509 warrant further source review. The Ye allele direction is inferred; ABFs use squared Z and do not depend on sign.

The custom calculations agree on identical inputs with the `approx.bf.estimates` and `combine.abf` reference functions from [coloc v5.2.3](https://github.com/chr1swallace/coloc/blob/v5.2.3/R/claudia.R), with maximum posterior difference below 7×10⁻¹⁵ ([check](stage3_rs323509_testis_coloc_formula_check.csv)). This checks formulas, not source preprocessing or biological assumptions. Follow-up requires a larger matched region, ancestry-matched LD, multi-signal analysis, sensitivity to priors and an independent insomnia dataset; no outcome-level risk estimate is derived here.

A subsequent [comparison using the original Atkins FI and Jansen insomnia GWAS](stage3_original_gwas_coloc_review.md) reports separate pairwise results and explains phenotype overlap, prior sensitivity and credible-set coverage.
