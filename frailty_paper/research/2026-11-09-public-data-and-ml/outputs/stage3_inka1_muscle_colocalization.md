# INKA1 muscle eQTL: exploratory regional comparisons

For the 2-Mb INKA1 skeletal-muscle region, single-signal approximate-Bayes-factor comparisons favor **different regional QTL and GWAS signals** for Mak FI, Mak modified physical frailty, and Ye FFS under the illustrative model below. FinnGen insomnia has weaker regional evidence: the QTL-only hypothesis dominates. These are exploratory results with incomplete molecular credible-set coverage and a complex, gene-rich region; they are not definitive evidence that every underlying signal is different, and they do not test a three-trait shared variant.

| GWAS | Matched variants | Muscle molecular credible-set variants covered | H1: QTL only | H3: different signals | H4: shared signal |
|---|---:|---:|---:|---:|---:|
| Mak FI | 2,085 | 73/88 | <0.0001 | 0.9999999 | 0.00000009 |
| Mak modified physical frailty | 2,085 | 73/88 | <0.0001 | 0.999947 | 0.000034 |
| Ye FFS | 2,097 | 74/88 | 0.000129 | 0.999824 | 0.000048 |
| FinnGen R13 insomnia | 2,048 | 70/88 | 0.8386 | 0.0724 | 0.0889 |

H0/H2 are negligible in the displayed model. These are posterior probabilities conditional on the tested model and priors, not p-values. For p1 = p2 = 10⁻⁴ and p12 = 10⁻⁵, the expression-QTL prior SD is 0.15 on an assumed standardized scale; GWAS prior SDs are 1.1004 (Mak FI), 0.12325 (Mak physical frailty), 0.12051 (Ye FFS) and 0.20 log-odds (FinnGen insomnia). The Mak scales were working estimates from another regional analysis. The Ye scale was estimated from SE, published total N = 386,565 and **QTL-sample MAF as a proxy** for unavailable FFS GWAS allele frequencies; this is an approximation, not a verified phenotype SD. Across tested GWAS effect-prior multipliers 0.5–2 and p12 = 10⁻⁶–5×10⁻⁵, insomnia H4 ranges 0.0051–0.4743. [Sensitivity results](stage3_inka1_muscle_coloc_sensitivity.csv).

One reported muscle QTL credible set contains 88 distinct variants, but 14–18 are missing from the matched analyses. Its rs34290943 exact-marker association and fine mapping are described in the [molecular review](stage3_inka1_muscle_finemapping.md). There may be more than one trait or QTL association signal in this region even though only one credible set was reported for this molecular phenotype. A single-signal model can produce misleading H3/H4 if variants are missing or multiple signals exist. Improved variant coverage, ancestry-matched LD, independently confirmed phenotype scales, and multi-signal fine-mapping are needed before strong interpretation.

## Inputs and checks

QTL input: eQTL Catalogue r7 QTD000281 GTEx v8 skeletal-muscle INKA1 expression, GRCh38 chr3:48803261–50803261; 702 QTL samples in source metadata ([receipt](stage3_catalogue_inka1_2mb_extract_receipt.json)). The GWAS are Mak FI and modified physical frailty, Ye FFS, and FinnGen R13 insomnia. Mak FI and physical frailty share a cohort; Ye FFS also uses UK Biobank, so the frailty results are correlated. This is not comparison with the Atkins FI and Jansen insomnia effects displayed in the main candidate table.

We required unique rsID, coordinate and allele-pair agreement, with unique positive-strand GRCh37-to-GRCh38 liftover for the three UK Biobank GWAS. The analysis retains non-palindromic biallelic SNPs with QTL MAF ≥1% and finite effect/SE; Mak and FinnGen GWAS allele frequency was required to be between 1% and 99%. Ye source filename indicates MAF ≥1%, but its per-variant frequency was unavailable for independent checking. The Catalogue muscle QTL imputation r² is `NA` for the candidate, so an r² filter was not applied; GWAS imputation quality was not independently verified. Missing alleles, indels, palindromic SNPs and nonunique identifiers are excluded. Effect allele was aligned to QTL ALT, but Ye A1 coding is inferred from the project source; ABF calculations use squared Z and do not depend on its sign. [Harmonization receipt](stage3_inka1_harmonization_receipt.json) and [coverage](stage3_inka1_muscle_coloc_coverage.json).

On identical aligned inputs and baseline priors, the calculations agree with `approx.bf.estimates` and `combine.abf` from [coloc v5.2.3](https://github.com/chr1swallace/coloc/blob/v5.2.3/R/claudia.R) to < 3×10⁻¹⁶ ([formula check](stage3_inka1_coloc_formula_check.csv)). This is a numerical check, not independent validation of preprocessing, scales, single-signal assumptions, source data or causality.

A subsequent [comparison using the original Atkins FI and Jansen insomnia GWAS](stage3_original_gwas_coloc_review.md) reports separate pairwise results and explains phenotype overlap, prior sensitivity and credible-set coverage.
