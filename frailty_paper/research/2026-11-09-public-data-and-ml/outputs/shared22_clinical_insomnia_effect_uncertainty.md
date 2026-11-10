# Insomnia odds ratios and uncertainty for the 22 candidates

The [candidate table](shared22_clinical_insomnia_effect_uncertainty.csv) puts the two clinical insomnia sources on the same counted-allele orientation. It reports the source-rounded MVP odds ratio and the FinnGen odds ratio with an approximate 95% Wald interval. Each is a **per-copy relative odds estimate**, conditional on the source GWAS model. It is not an absolute probability, a lifetime risk or a risk estimate for the combined 22-marker count.

| Candidate and counted allele | MVP OR, approximately; p | FinnGen OR (approx. 95% Wald CI); source p |
|---|---|---|
| rs12536335 A (FOXP2-labeled) | 1.032; 4.57 × 10⁻⁸ | 1.0185 (1.0052–1.0319); 0.00625 |
| rs28528997 A (MAD1L1-labeled) | 1.036; 3.24 × 10⁻⁷ | 1.0052 (0.9908–1.0198); 0.484 |
| rs2815764 G (NEGR1-labeled) | 1.041; 1.13 × 10⁻⁷ | 1.0042 (0.9854–1.0233); 0.667 |
| rs2155292 G (NCAM1-labeled) | 1.026; 3.89 × 10⁻⁶ | 1.0229 (1.0090–1.0369); 0.00116 |

The FinnGen intervals for rs28528997 and rs2815764 include both a small decrease and an increase in odds. They do not prove that the corresponding MVP estimates are wrong or that the variant has no effect. A comparison between sources would require a formal, appropriately designed analysis of differences, overlapping participants, allele identity, phenotype ascertainment and uncertainty. The beta estimates and source p-values can arise from tests that are not exactly Wald tests, so the reconstructed Wald interval and reported p-value need not agree perfectly.

All 22 candidates were already selected partly on their MVP insomnia associations; all 22 FinnGen point estimates happen to be positive on the same counted-allele orientation, and 9 have unadjusted FinnGen p < 0.05. This is a post-selection description, not 22 independent validations or a calibrated test of directional sharing. The variants include correlated markers and the full pattern can be affected by prior selection, overlap, population differences, and allele harmonization. Nearby genes are positional labels, not causal assignments.

## Calculation and availability

The inspected MVP GWAS Catalog release `GCST90475826` records `#NA` in the standard-error field for all 22 matched candidates. Its odds ratios are source-rounded; reversing the reported effect allele uses the reciprocal OR. We did not infer confidence intervals from a rounded OR and a p-value, because that would require an unverified Wald-test equivalence and could misstate uncertainty. The source [MVP study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12857194/) describes SAIGE-based association analysis; that fact alone does not establish the exact test for each row in this extracted file. A source supplying the actual log-OR standard errors would be needed for reliable conventional intervals.

FinnGen log-odds beta and SE come from the aligned project matrix. The displayed OR and bounds are `exp(beta)` and `exp(beta ± z₀.₉₇₅ × SE)`; this is an ordinary approximate Wald interval, not a multiple-testing- or selection-adjusted interval. The input source beta and SE, exact alleles and p-values are retained in the CSV. The arithmetic was checked for all 22 rows, and the MVP reciprocal orientation was checked against the prior source extraction. Assumed build conversion, nonpalindromic allele matching, possible multiallelic representations and sample-overlap limitations still apply; [screen provenance](song124_mvp_hfrs_followup.md) and [eight-source review](shared22_eight_source_evidence_review.md).
