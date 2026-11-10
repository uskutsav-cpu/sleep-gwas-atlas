# FOXP2 pituitary colocalization with the insertion included

The earlier [FOXP2 analysis](stage3_foxp2_pituitary_colocalization.md) excluded the molecular credible-set insertion chr7:114534306 T>TA. The Stage 5 reciprocal mapping and reference checks identified an exact representation in Mak FI, Mak modified physical frailty, Ye FFS and FinnGen R13 insomnia. This sensitivity analysis adds that one insertion to each earlier matched regional table; all other inclusion and exclusion decisions are unchanged. The insertion's two eQTL Catalogue rsIDs were deduplicated by exact variant and identical molecular estimates. The [coverage table](stage5_foxp2_indel_coloc_coverage.csv) and four `stage5_foxp2_pituitary_*_with_indel_aligned.csv` files preserve the inputs. All four unique variants in the reported molecular credible set are now included, compared with three before; this does not establish coverage of every causal variant or signal in the locus.

Under the same single-signal approximate Bayes-factor model and baseline priors as the earlier analysis (p1 = p2 = 10⁻⁴, p12 = 10⁻⁵; QTL prior SD 0.15; original trait-specific scales), the posteriors are:

| Trait | Matched variants | H3: distinct QTL and GWAS signals | H4: shared signal |
|---|---:|---:|---:|
| Mak FI | 3,231 | 0.99978 | 0.00021 |
| Mak physical frailty | 3,231 | 0.99541 | 0.00038 |
| Ye FFS | 3,184 | 0.96421 | 0.00149 |
| FinnGen R13 insomnia | 3,180 | 0.66380 | 0.00795 |

H3 and H4 do not sum to one: the remaining posterior mass includes no association or association with only one trait. In these baseline models, adding the insertion does not materially change the earlier finding favoring distinct signals over a shared signal. Across the tested prior scales and shared priors, H4 ranges from 0.00002–0.00108 (Mak FI), 0.00004–0.00191 (Mak physical frailty), 0.00014–0.00758 (Ye FFS), and 0.00060–0.04670 (FinnGen). The full [sensitivity table](stage5_foxp2_indel_coloc_sensitivity.csv) is exploratory, not a proof that sharing is absent.

The insertion's pituitary QTL is based on 237 samples; GTEx molecular data can overlap model training, and some frailty GWAS share participants. The single-causal-variant assumption, excluded noncredible-set variants and other indels, variant coverage, missing imputation-quality values, phenotype scales, chosen priors and possible multiple signals remain important limitations. Mak sources' GRCh37 build and effect allele are inherited assumptions. Ye's A1 direction remains inferred. Since ABFs here depend on squared Z, a sign flip alone does not change these H0–H4 values, but incorrect variant matching or effect scales could. The result does not establish that the insertion affects FOXP2 or causes insomnia/frailty.

The original baseline was recomputed before adding the insertion and agreed with the earlier output. The new baseline calculations were checked against `coloc` v5.2.3 formulas to maximum difference below 2 × 10⁻¹⁵ ([formula checks](stage5_foxp2_indel_coloc_formula_check.csv)). This checks arithmetic, not source correctness or model assumptions. Separate full-precision AlphaGenome indel scoring and allele-specific experimental validation are still missing.
