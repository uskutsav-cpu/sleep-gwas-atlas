# FOXP2 pituitary expression: exploratory regional comparisons

A 2-Mb single-signal comparison of the pituitary FOXP2 eQTL with Mak FI, Mak modified physical frailty, Ye FFS and FinnGen insomnia favors different molecular and trait signals under the baseline assumptions. **This cannot rule out a shared signal:** the analysis excludes the insertion in the molecular credible set, which has molecular PIP 0.312, and may miss other variants or multiple signals.

| GWAS | Matched SNPs | Molecular credible-set variants represented | H1: QTL only | H3: distinct signals | H4: shared signal |
|---|---:|---:|---:|---:|---:|
| Mak FI | 3,230 | 3/4 | <0.001 | 0.9999 | 0.00010 |
| Mak modified physical frailty | 3,230 | 3/4 | 0.0042 | 0.9955 | 0.00033 |
| Ye FFS | 3,183 | 3/4 | 0.0343 | 0.9642 | 0.00146 |
| FinnGen R13 insomnia | 3,179 | 3/4 | 0.3282 | 0.6637 | 0.00806 |

H0/H2 are small under this illustration. The three retained SNPs account for approximately 0.661 of the reported 0.972 summed PIPs in this reported molecular set; the missing insertion accounts for 0.312. This is molecular-set coverage, not GWAS fine-mapping or proof that the set is exhaustive. The single-signal model, missing indel, multiple signals, unverified source overlap and scale assumptions materially limit inference. The three project markers themselves are not in the reported molecular set ([fine-mapping review](stage3_foxp2_pituitary_finemapping.md)).

The baseline uses p1 = p2 = 10⁻⁴, p12 = 10⁻⁵, QTL prior SD 0.15 and trait prior SD 1.1004 (Mak FI), 0.12325 (Mak FP), 0.12225 (Ye FFS) or 0.20 log-odds (insomnia). The Mak SDs are working estimates from a prior region, and Ye's SD uses SE, total N = 386,565 and QTL MAF as a proxy for missing GWAS frequencies. Varying trait SD by 0.5–2 and p12 from 10⁻⁶ to 5×10⁻⁵ gives H4 ranges 0.000010–0.000509 for Mak FI, 0.000033–0.00165 for Mak FP, 0.00014–0.00741 for Ye FFS and 0.00061–0.0473 for insomnia ([full sensitivity](stage3_foxp2_pituitary_coloc_sensitivity.csv)). Those ranges do not address the omitted indel or multiple signals.

The molecular source is eQTL Catalogue r7 QTD000301, reprocessed GTEx v8 pituitary expression, metadata N = 237, GRCh38 chr7:113086327–115086327 ([extraction receipt](stage3_catalogue_foxp2_pituitary_extract_receipt.json)). Matching used unique rsID, position and allele pair; UK Biobank data were uniquely lifted from GRCh37 on the positive strand. Non-palindromic biallelic SNPs and QTL MAF ≥1% were retained; Mak/FinnGen allele frequency ≥1% was also required. Ye frequency was unavailable for variant-level checking; its source filename indicates a 1% MAF filter. Molecular and GWAS imputation quality and sample overlap need further review. The Mak and Ye frailty studies reuse UK Biobank participants, so they are not independent. [Matching](stage3_foxp2_pituitary_harmonization_receipt.json), [coverage](stage3_foxp2_pituitary_coloc_coverage.json).

The custom H0–H4 formulas agree with [coloc v5.2.3](https://github.com/chr1swallace/coloc/blob/v5.2.3/R/claudia.R) functions on the same inputs to <4×10⁻¹⁵ ([check](stage3_foxp2_pituitary_coloc_formula_check.csv)). That validates the calculation, not input preparation, prior choice or biological inference. Indel-inclusive, ancestry-matched multi-signal fine-mapping is a necessary follow-up.

A subsequent [comparison using the original Atkins FI and Jansen insomnia GWAS](stage3_original_gwas_coloc_review.md) reports separate pairwise results and explains phenotype overlap, prior sensitivity and credible-set coverage.
