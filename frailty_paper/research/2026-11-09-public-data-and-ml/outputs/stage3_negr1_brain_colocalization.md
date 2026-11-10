# rs2815764: NEGR1 caudate expression fallback

The initial worklist selected a cerebellar RPL31P12 eQTL for rs2815764, but the queried Catalogue r7/v8 dataset returned no rows for that gene. As a separate fallback, I examined **NEGR1** (ENSG00000172260) expression in brain caudate: GTEx v10 reports the candidate eQTL p = 1.84 × 10⁻⁴, and the related Catalogue r7 QTD000156 GTEx v8 reprocessing reports p = 5.63 × 10⁻⁵ for the exact `chr1:72288440 A/G` allele. These are related source data and not independent replications. The fallback is a different gene/tissue choice, selected after the original query returned no rows.

The Catalogue's caudate molecular fine-mapping file reported **no NEGR1 credible set**. The following single-signal comparisons therefore have unknown molecular-signal completeness and are screening results, not evidence that NEGR1 mediates a trait or that a shared causal variant exists.

| GWAS | Matched SNPs | H1: QTL only | H3: distinct signals | H4: shared signal |
|---|---:|---:|---:|---:|
| Atkins FI | 3,170 | 0.669 | 0.219 | 0.085 |
| Jansen insomnia | 3,239 | 0.778 | 0.175 | 0.018 |
| Mak FI | 3,324 | <0.001 | 0.511 | 0.474 |
| Mak modified physical frailty | 3,324 | 0.053 | 0.713 | 0.210 |
| Ye FFS | 3,307 | 0.095 | 0.728 | 0.152 |
| FinnGen R13 insomnia | 3,183 | 0.763 | 0.078 | 0.131 |

H0/H2 are omitted from the table and are available in the [full model output](stage3_negr1_coloc_sensitivity.csv). For Mak FI, H4 ranges 0.083–0.818 across the tested trait-prior multipliers and p12 settings; for Atkins FI it ranges 0.005–0.424. Different phenotypes, samples, variant support and priors can contribute to differing results; no direct difference test was performed. No claim of shared insomnia–frailty risk follows. One reported set being absent does not prove that there is no molecular signal; it makes credible-set coverage unavailable for this evaluation.

The QTL input is Catalogue r7 QTD000156, GTEx v8 caudate expression (194 samples in metadata), in a 2-Mb candidate-centered GRCh38 window ([extraction receipt](stage3_catalogue_negr1_caudate_extract_receipt.json)). Matching uses the public [regional GWAS caches](stage3_gwas_regional_cache_receipt.json), unique rsIDs, positive-strand unique liftover for GRCh37 sources, exact position/allele pairs, non-palindromic biallelic SNPs and QTL MAF ≥1%. Variant-specific GWAS MAF was checked where available; Ye's was not available. QTL and GWAS imputation quality, ancestry matching and sample overlap remain unverified. Mak/Ye/Atkins/Jansen use overlapping or potentially overlapping UK Biobank participants. [Matching receipt](stage3_negr1_qtl_harmonization_receipt.json), [model coverage](stage3_negr1_coloc_coverage.json).

The baseline prior is p1 = p2 = 10⁻⁴, p12 = 10⁻⁵, QTL effect SD 0.15 and a phenotype-specific trait effect SD. The sensitivity grid varies the trait SD by 0.5–2 and p12 from 10⁻⁶ to 5×10⁻⁵; it does not resolve missing molecular fine mapping, prior misspecification or multiple signals. The custom ABFs agree with [coloc v5.2.3](https://github.com/chr1swallace/coloc/blob/v5.2.3/R/claudia.R) reference functions for all six baseline tests to <9×10⁻¹⁵ ([check](stage3_negr1_coloc_formula_check.csv)). That is an arithmetic check, not validation of the biological assumptions.
