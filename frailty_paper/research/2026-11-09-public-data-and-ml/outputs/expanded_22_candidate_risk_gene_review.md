# Twenty-two candidate variants: shared-risk follow-up and nearby genes

This is a research shortlist for discussing insomnia–frailty risk, **not a validated risk score**. It contains all 22 of the 109 eligible Song et al. candidates with nominal p<0.05 in both the MVP insomnia and Finnish hospital-frailty (HFRS) GWAS. These markers were selected after prior cross-trait results and follow-up; they are not a prespecified or independent discovery set. All have the same estimated direction in the two follow-up outcomes, but same-sign point estimates do not establish a shared causal variant or mechanism.

**Allele orientation:** the counted alleles in this original evidence table match the earlier trait-weighted prototype, not necessarily the later [shared-direction count](shared22_score_research_guide.md). Eight of 22 are reversed between those two files. Both orientations are valid only with their corresponding effect signs: for example, rs12536335 is G in this table but A in the shared-direction count. Use the [allele orientation crosswalk](shared22_allele_orientation_crosswalk.csv) when comparing files; do not apply a coefficient to the other allele without reversing its sign under the verified biallelic diploid model.

| Variant | Counted allele | Nearest gene-body label | FI p | FFS p | Larger of MVP insomnia and HFRS p-values |
|---|---|---|---:|---:|---:|
| rs12536335 | G | FOXP2 | 2.63e-06 | 0.0017 | 4.57e-08 |
| rs28528997 | G | MAD1L1 | 8.04e-05 | 2.3e-12 | 2.16e-06 |
| rs2815764 | G | NEGR1 | 0.000136 | 1.2e-05 | 2.3e-06 |
| rs2155292 | G | NCAM1 | 2.75e-08 | 0.54 | 3.89e-06 |
| rs817509 | G | RP11-384F7.1 | 0.0223 | 0.0099 | 2.86e-05 |
| rs1909650 | A | PCDH8P1 | 6.32e-07 | 0.097 | 0.000386 |
| rs2838788 | C | ADARB1 | 0.000791 | 0.004 | 0.00188 |
| rs34290943 | C | LAMB2 | 2.6e-05 | 0.00011 | 0.00401 |
| rs12671330 | C | FOXP2 | 2.02e-09 | 0.0017 | 0.00718 |
| rs9536347 | T | OLFM4 | 0.00779 | 0.0038 | 0.00901 |
| rs6452792 | A | CTC-498M16.4 | 0.00166 | 0.36 | 0.0121 |
| rs7691487 | C | MMRN1 | 0.000841 | 0.0021 | 0.0139 |
| rs4886140 | G | RPP40P2 | 7.89e-05 | 0.45 | 0.0154 |
| rs2604551 | G | RP11-665G4.1 | 0.000363 | 0.38 | 0.016 |
| rs7734045 | C | AC009014.3 | 0.00545 | 0.015 | 0.0184 |
| rs6058750 | C | RP11-410N8.1;RP11-410N8.4 | 8.94e-06 | 0.77 | 0.026 |
| rs9845387 | A | LSAMP | 0.00922 | 0.21 | 0.0263 |
| rs12524473 | G | CD83 | 0.000302 | 0.36 | 0.0272 |
| rs17312861 | G | FOXP2 | 0.0032 | 0.39 | 0.0355 |
| rs6938026 | G | RP3-330M21.5;CUL9 | 0.000955 | 0.021 | 0.0413 |
| rs7731706 | G | TRPC7 | 0.000333 | 0.064 | 0.0478 |
| rs11636381 | C | RP11-24D15.1 | 0.06 | 7.9e-08 | 0.0493 |

The labels identify the nearest annotated gene body within ±100 kb of the position assumed to be GRCh37. They do not establish an effector gene. Some variants overlap noncoding RNAs or pseudogenes; a neighboring protein-coding gene or a distant regulated gene could be relevant instead. For example, our nearest-body lookup labels rs34290943 as LAMB2, while Song et al. discussed DALRD3 in their broader gene analyses; these represent different mapping approaches, not a confirmed assignment of either gene to this variant.

## What the list can and cannot support

- Six rows meet the numerical max(p) × 124 < 0.05 benchmark in the earlier exploratory screen. The benchmark is not a valid post-selection confirmation without additional assumptions and an independent prespecified analysis.
- Eleven variants have matching point-estimate signs and nominal p<0.05 across the reviewed FI, FFS, UKB insomnia, MVP insomnia and HFRS inputs. These observations are affected by the selection, multiple testing, overlapping or unknown participants and phenotype differences. The FFS effect-allele convention is inferred. Failure of a nominal threshold is not evidence of no effect.
- All 22 specific two-allele contrasts matched biallelic records in the 503-person European 1000 Genomes reference panel. Fourteen rsIDs have more than two alleles in broader Ensembl GRCh37 records, so this check does not settle original study-specific coding or establish that every source used an equivalent contrast.
- Five same-chromosome pairs within 1 Mb were evaluated in the reference panel. For the tested contrasts, their haplotype r² values range from 0.031 to 0.160; the highest was rs1909650–rs9536347. These values are specific to this panel, and low pairwise r² does not establish conditional independence or shared causality. Long-range LD, other variants, and population differences were not exhaustively evaluated.
- The originally prepared 15-marker FI, FFS and insomnia weighted sums are separate technical prototypes, not absolute risk probabilities or evaluated predictors. This 22-marker set was selected using MVP and HFRS results and has **not** been fitted, calibrated, or evaluated against independent outcomes. Do not add FI, FFS, HFRS and insomnia coefficients, which use different scales.

## Recommended next analysis

Freeze a study question and candidate-selection rule before looking at a new validation cohort. Verify alleles and ancestry-matched LD, consider locus-level fine-mapping and colocalization, and link variants to genes with regulatory evidence. For a 15–30-variant score, choose one outcome and time horizon, select and estimate weights in a development set, then evaluate it in independent genotyped participants against clinical baseline, FI and FFS scores. Assess calibration and incremental prediction; perform a sleep/fatigue-excluded FI sensitivity if data permit. No such outcome-based evaluation was possible with the reference genotypes used here.

## Files and provenance

`expanded_22_candidate_risk_gene_table.csv` retains coordinates, the counted/other alleles, signed source effects, p-values, gene annotations, and reference flags. Its SHA-256 is `d500617af1a76f8d7408c26a2b4ba277a3f736335aab7efd3f364459872fc874`. `song124_22_nominal_positional_genes.csv` preserves Ensembl annotation URLs and response hashes. `song124_22_nominal_reference_ld.json` records reference-panel matching and the five local pairwise LD estimates. The input screen and assumptions are documented in `song124_mvp_hfrs_followup.md`; the presumed GRCh37 build of the Song table remains not independently authenticated. No individual-level genotypes or scores were saved.
