# Brain6 insomnia and ADHD versus frailty: genome-wide LDSC follow-up

A new post-hoc LDSC analysis of the original Brain6 insomnia input against the native frailty-index (FI) input estimated **r_g=0.6405 (SE 0.0235)** and a **cross-trait covariance intercept of 0.1593 (SE 0.0063)**. The intercept suggests nonzero shared test-statistic structure under the LDSC model. Both studies include UK Biobank and FI includes an insomnia item, but the intercept does not identify the number of shared participants, establish their phenotypic correlation, or uniquely distinguish sample overlap from other causes.

| Pair | Valid merged SNPs | r_g (SE) | p | Bonferroni p, these 3 post-hoc pairs | Cross-trait intercept (SE) |
|---|---:|---:|---:|---:|---:|
| Original Brain6 insomnia–FI | 1,130,007 | 0.6405 (0.0235) | 4.32×10⁻¹⁶³ | 1.30×10⁻¹⁶² | 0.1593 (0.0063) |
| Original Brain6 insomnia–FFS | 1,130,007 | 0.4839 (0.0266) | 3.23×10⁻⁷⁴ | 9.70×10⁻⁷⁴ | 0.1383 (0.0065) |
| Brain6 ADHD–FI | 1,087,104 | 0.5467 (0.0285) | 3.15×10⁻⁸² | 9.45×10⁻⁸² | 0.0072 (0.0069) |

The p-values are those printed in the detailed LDSC sections; the summary table rounds them to zero. The correction is limited to this three-pair follow-up, not all tests in the atlas. The near-zero ADHD–FI intercept does not prove sample independence. SNP support differs across pairs and the estimates have correlated uncertainty; do not infer a statistically significant difference from their point estimates. They reuse the original GWAS and are not independent replication.

These results concern the original Brain6 insomnia phenotype; the earlier [three-frailty comparison](../insomnia_frailty_same_support/RESULTS.md) used FinnGen registry insomnia and common support. Differences between them cannot be attributed to a single cause, since phenotypes, cohorts and SNP support differ. FI includes sleep and tiredness items and FFS includes exhaustion, so these results do not separate sleep-independent physical frailty.

## Relation to the chromosome 7 pilot

The positive insomnia–FI intercept strengthens the reason to question the pilot's independent-statistics assumption. It **cannot be inserted directly** as the HyPrColoc correlation parameter: LDSC's genome-wide regression intercept is not a validated local residual-correlation estimate, and trait scale, overlap weighting, ancestry, LD and model mismatch need evaluation. The existing [hypothetical correlation stress test](../../frailty_coloc/chr7_three_trait/overlap_sensitivity/RESULTS.md) does not become an overlap-corrected analysis. Among 1,055 selected pilot markers, 225 have valid uniquely matched FI records in the native munged file; all 225 have concordant signs after allele alignment. This limited check supports consistency of those values, but does not prove full source identity or agreement for missing variants. ADHD–FI's near-zero fitted intercept likewise does not certify independence.

The [protocol](PROTOCOL.md), [run script](run.sh), complete logs, [table](results.tsv), [verification script](verify.py), and [receipt](receipt.json) document the calculations. Input hashes match the earlier same-data rechecks, and the LDSC script and 44 LD-score/reference files match their earlier receipt. This is computational provenance, not validation of all underlying phenotype, allele, reference or participant assumptions. These estimates do not rescue failed LAVA QC, validate PLACO+, establish an effector gene or demonstrate a causal mechanism.

A later [common-support paired jackknife](paired_common_support/RESULTS.md) quantifies the FI–FFS contrast for original Brain6 insomnia; its limitations include sleep-item overlap and post-hoc selection.
