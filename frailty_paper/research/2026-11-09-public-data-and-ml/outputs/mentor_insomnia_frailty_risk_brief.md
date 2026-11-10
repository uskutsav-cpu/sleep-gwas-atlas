> **Updated synthesis:** For the current findings, validation limits and AI plan, see [the 28 October mentor synthesis](mentor_insomnia_frailty_synthesis_2026-10-28.md).

# Insomnia and frailty: risk-focused mentor brief

Prepared 27 October 2026. This brief summarizes work using existing project results and free public data. The estimates are research findings, not an individual's risk of disease.

## The question and why it matters

Insomnia and frailty are associated in population studies, but frailty can be defined in different ways. The **frailty index (FI)** aggregates health deficits across domains; the **Fried frailty score (FFS)** counts five physical domains: weakness, slow walking, exhaustion, low activity and weight loss. The project's UK Biobank FI includes a sleeplessness question and tiredness, whereas FFS includes exhaustion. This creates a direct measurement-overlap concern when studying sleep and frailty. Removing an item can also change the phenotype itself, so the overlap does not invalidate the question, but it must be tested rather than assumed away. [Definitions and source-specific caveats](fi_ffs_component_crosswalk.md).

## FI versus FFS

The project's linkage disequilibrium score regression (LDSC) analysis of the Atkins FI and Ye FFS GWAS estimated **genome-wide genetic correlation r_g = 0.757** (SE 0.021; approximate 95% CI 0.715–0.799). A second analysis using Mak's FI and modified Fried physical phenotype estimated **r_g = 0.796** (SE 0.012). The studies overlap or may overlap through UK Biobank and use different definitions. These are correlations of genetic effects across variants, not correlations of people's clinical measurements, proportions of genes shared, or causal effects. [Methods and matched comparison](mak2023_same_cohort_fi_fp_ldsc_review.md).

In a separate clinical cohort, a 27-item FI and Fried-criteria count had a reported Spearman correlation of **0.542**. This is an individual-level correlation in a different sample; it cannot be compared numerically with genome-wide r_g as if the two measured the same thing. [Clinical evidence](fi_ffs_same_cohort_comparison.md).

## Insomnia and the two frailty measures

| Public insomnia GWAS | Insomnia–FI r_g | Insomnia–FFS r_g | Paired difference, FFS minus FI (95% CI) |
|---|---:|---:|---:|
| UK Biobank symptom measure | 0.640 | 0.483 | −0.156 (−0.213 to −0.100) |
| FinnGen clinical insomnia | 0.426 | 0.245 | −0.181 (−0.260 to −0.102) |
| MVP clinical insomnia, European ancestry | 0.549 | 0.398 | −0.151 (−0.220 to −0.082) |

FI shows the higher estimate in all three comparisons. Paired uncertainty uses matched SNP support and block-jackknife covariance. These analyses reuse the same FI and FFS GWAS and have differing phenotype definitions, ascertainment, LD references and possible sample overlap. The pattern motivates a study of overlap and shared genetic architecture; it cannot tell us which mechanism explains the difference or establish a causal direction. [Analysis and diagnostics](insomnia_frailty_risk_briefing.md).

## Shared markers and a research score

Song et al.'s prior cross-trait study supplied 124 candidate variants. Among 109 eligible matches in our subsequent screen, 22 showed nominal p < 0.05 in both MVP insomnia and Finnish hospital frailty (HFRS), with estimated directions aligned after allele orientation. This screen was post-selected; traits, participants and studies overlap or may overlap. HFRS is not FI or FFS, and joint marginal association does not prove a shared causal variant. [Selection details](song124_mvp_hfrs_followup.md).

The [candidate table](shared22_mentor_risk_table.md) gives each allele, nearest gene-body label and separate evidence across traits. Examples include markers near **FOXP2, MAD1L1, NEGR1 and ADARB1**. These are positional labels rather than established effector genes. Allele/build and multiallelic coding require further verification. Five of the 22 were nominal with positive estimates in all eight inspected GWAS result sets, but the sets are not eight independent replications and nearby markers can be correlated.

The [22-marker prototype](shared22_score_research_guide.md) is an equal-weight count from 0 to 44 for complete diploid calls. It was not trained or tested against observed outcomes and does not produce an absolute probability. Removing the three FOXP2-labeled markers changes the reference-panel count (correlation 0.898 with the full count in 503 EUR-labeled people); this is a genotype-only sensitivity check, not prediction or evidence that FOXP2 is causal. [Reference sensitivity](shared22_leave_region_reference_review.md).

## Next steps with free public data

1. **Test the measurement explanation.** Search for full, publicly accessible GWAS statistics from an FI excluding insomnia or fatigue. Mak et al. report a modified FI without eight mental-wellbeing items, including those two, but its full genome-wide file has not been located in the inspected release. The already analyzed diagnosis-based HFRS offers different-phenotype triangulation, though it includes some sleep-adjacent diagnostic categories. [Available options](public_frailty_sleep_overlap_sensitivity_options.md).
2. **Resolve variant identity and shared signals.** Verify build, allele and multiallelic representation and use appropriate ancestry-matched LD. When suitable full regional statistics and LD are available, apply quality-controlled fine-mapping and colocalization; an association at the same marker alone does not establish the same causal variant. The project's local genetic-correlation attempts failed their quality gates and should not be reported as negative findings.
3. **Assess score performance only if eligible public outcome data are found.** Freeze candidates and alleles, establish participant and relatedness independence from selection and GWAS, and compare association, calibration and incremental prediction with clinical predictors and published genome-wide scores. Public summary statistics alone cannot validate the prototype for individual prediction. No suitable fully public participant-level validation cohort has yet been established.
4. **After the risk analysis, follow function.** Investigate candidate-gene and cell-type evidence, and use sequence models such as AlphaGenome to generate testable variant-effect hypotheses. Gene proximity or model predictions do not establish mechanism.

The [detailed risk assessment](mentor_risk_assessment_updated.md) contains source-specific caveats, additional sensitivity analyses and links to the underlying result files.
