# ADARB1 muscle splicing: exploratory colocalization

For the selected ADARB1 skeletal-muscle splice junction, the single-signal model at the baseline priors below is dominated by **QTL-only H1** for Mak FI, Mak modified physical frailty and Ye FFS. The FinnGen insomnia comparison gives H1 = 0.339, H3 = 0.556 and H4 = 0.105. These results do not establish shared insomnia–frailty genetics or exclude a shared signal; H4 depends on priors, and two molecular credible-set variants are missing.

| GWAS | Matched SNPs | Molecular credible-set coverage | H1: QTL only | H3: distinct signals | H4: shared signal |
|---|---:|---:|---:|---:|---:|
| Mak FI | 3,230 | 13/15 | 0.897 | 0.083 | 0.020 |
| Mak modified physical frailty | 3,230 | 13/15 | 0.882 | 0.072 | 0.046 |
| Ye FFS | 3,181 | 13/15 | 0.843 | 0.079 | 0.079 |
| FinnGen R13 insomnia | 3,109 | 13/15 | 0.339 | 0.556 | 0.105 |

H0 and H2 are negligible under these assumptions. Separately, the [original Atkins FI and Jansen insomnia comparisons](stage3_original_gwas_coloc_review.md) also have high H1 (0.938 and 0.949, respectively). No direct three-trait colocalization or causal mediation test was performed.

The baseline sets p1 = p2 = 10⁻⁴, p12 = 10⁻⁵, a molecular prior SD of 0.15 and trait prior SDs 1.1004 (Mak FI), 0.12325 (Mak FP), 0.12121 (Ye FFS) and 0.20 log-odds (FinnGen insomnia). The Mak values are working estimates from another region; Ye's estimate uses source SE, total N 386,565 and QTL MAF as a proxy for unavailable GWAS MAF. The molecular scale is assumed approximately standardized. When trait-prior SD is multiplied by 0.5–2 and p12 ranges 10⁻⁶–5×10⁻⁵, H4 ranges 0.0011–0.159 for Mak FI, 0.0025–0.310 for Mak FP, 0.0044–0.440 for Ye FFS and 0.0084–0.429 for FinnGen insomnia ([full table](stage3_adarb1_muscle_splice_coloc_sensitivity.csv)).

The source molecular credible set has 15 variants; each comparison covers 13. The two excluded variants are C/G palindromic SNPs with reported molecular PIPs 0.0226 and 0.0481; the 15 PIPs sum to approximately 0.962. They are excluded because strand resolution was not independently established, and missing variants and multiple trait or molecular signals can distort the single-signal H3/H4 comparison. The source molecular fine mapping used a 2-Mb region; this test used the 1-Mb candidate-centered window `chr21:44622459–45622459` GRCh38. The [fine-mapping review](stage3_adarb1_muscle_splice_finemapping.md) describes the specific junction and its relation to rs2838788.

## Inputs and checks

The molecular source is eQTL Catalogue r7 QTD000285, GTEx v8 reprocessed skeletal-muscle splicing, metadata N = 702. Its `cc` file contains selected molecular events and does not establish negative results for other junctions. [Extraction receipt](stage3_catalogue_adarb1_splicing_extract_receipt.json), [GWAS matching receipt](stage3_adarb1_muscle_splice_harmonization_receipt.json), and [coverage](stage3_adarb1_muscle_splice_coloc_coverage.json). Sample overlap and ancestry matching need independent review; the Mak traits share a cohort and the Ye study also uses UK Biobank participants.

Matching required unique rsID, chromosome, coordinates and allele pair, and positive-strand unique GRCh37-to-GRCh38 liftover for the UK Biobank studies. Non-palindromic biallelic SNPs with QTL MAF ≥1% and valid beta/SE were retained; Mak and FinnGen GWAS MAF ≥1% was checked. Ye per-variant GWAS MAF was unavailable, although its source filename specifies a 1% MAF filter. The candidate QTL imputation r² is missing, and no QTL r² filter was applied; GWAS imputation quality was not independently verified. Ye A1 direction is inferred, while the ABF uses squared Z and is not affected by sign.

On identical inputs and baseline priors, the custom formulas agree with [coloc v5.2.3](https://github.com/chr1swallace/coloc/blob/v5.2.3/R/claudia.R) reference functions to <4×10⁻¹⁵ ([check](stage3_adarb1_muscle_splice_coloc_formula_check.csv)). This checks numerical implementation only; it does not validate input processing, prior choices, signal multiplicity or tissue mechanism. More complete allele-resolved coverage and ancestry-matched multi-signal analysis remain necessary.
