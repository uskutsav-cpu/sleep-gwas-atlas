# The 22 shared-direction candidates across eight GWAS results

The [combined evidence matrix](shared22_eight_source_evidence_matrix.csv) aligns the 22 alleles in the unweighted research count to eight available result sets: Atkins FI, Ye FFS, Mak FI, Mak modified physical phenotype, FinnGen hospital frailty (HFRS), UKB insomnia, MVP insomnia and FinnGen R13 insomnia. It retains each source effect and p-value and, where available, the source standard error and approximate Wald interval. All effects are on their **original phenotype scales**. The MVP log-odds values are reconstructed from rounded source odds ratios and are approximate; coefficients from different traits cannot be added or compared as if they had the same units.

Five candidates have a positive point estimate and unadjusted p < 0.05 in all eight rows after alignment to the score's counted allele:

| Variant | Counted allele | Nearest gene-body label | FinnGen insomnia p | Mak physical-frailty p |
|---|---|---|---:|---:|
| rs12536335 | A | FOXP2 | 0.00625 | 0.00279 |
| rs817509 | G | RP11-384F7.1 | 0.00272 | 0.00534 |
| rs2838788 | T | ADARB1 | 0.000345 | 0.00813 |
| rs34290943 | T | LAMB2 | 0.00110 | 0.000198 |
| rs12671330 | C | FOXP2 | 0.000384 | 0.00152 |

This is a **descriptive shortlist**, not a set of confirmed shared causal alleles. The 22 candidates were selected using previous insomnia–FI analyses and the MVP-insomnia/HFRS results; using those same data in the table does not produce independent replication. UKB contributes to several FI, physical-frailty and insomnia GWAS, and FinnGen releases may overlap. Statistical tests and selection were repeated, so the nominal p-values and counts are not adjusted for the discovery process. The two FOXP2-labeled markers are in the same region and have nonzero reference LD, and should not be counted as independent signals. Positional gene labels do not identify an effector gene; rs34290943 illustrates that alternative mapping approaches can nominate different genes.

An [allele-specific transcript annotation](shared5_variant_annotation_review.md) is available for these five markers. It adds transcript context but does not establish gene function, a shared causal variant or an effector gene.

Variants not meeting every nominal threshold have not been shown to lack an effect. A positive point estimate alone is not a statistically established association, and comparing p-value thresholds across studies does not establish effect heterogeneity. Multiallelic representation, assumed genome build, the FI TwinGene sign convention and source-specific allele coding remain material issues. The FFS effect-allele convention is supported by the [source-table check](ffs_effect_allele_source_consistency_audit.md), not verified for every candidate. The FI versus FFS study definitions overlap with sleep/fatigue measures, and HFRS is a different, diagnosis-based construct. None of these results validates the [unweighted allele-count prototype](shared22_score_research_guide.md) as an individual risk score.

The matrix combines the previously documented [22-candidate screen](expanded_22_candidate_risk_gene_review.md), [FinnGen-insomnia follow-up](22_candidate_six_trait_followup.md), and [Mak FI/FP lookup](mak2023_candidate_consistency_review.md). The [receipt](shared22_eight_source_evidence_matrix.json) documents source names, units, counts and the output hash. Script: `work/build_shared22_evidence_matrix.py`. No individual-level genotype or outcome data were analyzed.

An [arithmetic audit](shared22_numerical_consistency_audit.json) recalculated two-sided normal-Wald p-values and stored Wald intervals for all 22 rows in each of Mak FI, Mak FP, HFRS and FinnGen insomnia, the four sources for which the matrix retains SEs. All 88 p-values agreed with the recalculation to within approximately 0.02% relative difference, and the intervals matched within floating-point rounding. This is a consistency check, not independent verification of source alleles, phenotype units, study design, selection or causal interpretation; p-values from other test types need not generally equal Wald approximations. It does not cover the four matrix sources without SEs.
