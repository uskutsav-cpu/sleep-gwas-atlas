# Insomnia and frailty: an updated assessment of shared genetic risk

Prepared for research discussion • 27 October 2026

## What can be said now

FI and FFS are substantially genetically correlated in the available GWAS, but they measure different constructs. Across three insomnia GWAS, the genetic correlation with FI is higher than that with FFS. This consistent pattern motivates a comparison of definitions and genetic architectures; it does not establish that insomnia causes frailty. FI contains sleep-related items, which makes phenotype overlap an important alternative explanation.

The project has an exploratory set of 22 variants with evidence in both an insomnia GWAS and a hospital-frailty GWAS, with nearby-gene labels and separate FI-, FFS- and insomnia-weighted sums. **Neither the variants nor the sums constitute a validated shared-risk predictor.** Use them to formulate testable questions, not to give anyone an absolute risk or diagnose a condition.

For a single consistently oriented view of the complete candidate set, see the [22-marker mentor table](shared22_mentor_risk_table.md).

## 1. Compare what the frailty measures capture

The frailty index (FI) operationalizes accumulated deficits across multiple health domains. The Fried frailty score (FFS) counts five physical domains: weakness, slow walking, exhaustion, low activity and weight loss. The original Fried definition refers to unintentional weight loss; the UK Biobank question in the Hanlon/Mak adaptation does not distinguish intentional from unintentional loss. Ye cites the Hanlon adaptation, but the complete Ye implementation needs verification. In the project's analyses FI is based on the [Atkins et al. GWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/), including UK Biobank and TwinGene; FFS is based on the [Ye et al. GWAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/). Their samples partly reuse UK Biobank, while exact participant intersection is not established. [Source crosswalk](mak_hanlon_physical_frailty_crosswalk.md).

There are two different comparisons:

- **Genetic correlation:** On common SNP support, LDSC estimates FI–FFS r_g = **0.757**, SE 0.021, approximately 95% CI **0.715–0.799**. It summarizes genome-wide correlation of estimated genetic effects for these GWAS phenotypes. It is not the clinical correlation between an individual's FI and FFS, a fraction of variants shared, or a causal estimate. [Detailed matrix](frailty_three_way_common_support_ldsc.md). The contributing LDSC heritability estimates are reasonably precise in this fit, and the freely estimated cross-trait intercept is nonzero; shared UK Biobank participants may contribute, but the intercept cannot determine overlap by itself. [Diagnostics](fi_ffs_ldsc_heritability_diagnostics.md).
- **Measured clinical correlation:** A separate small Spanish study reported Spearman **0.542** for its 27-item FI and Fried-criteria count. A different Canadian study reported phi **0.38** for binary frailty classifications. Those cohorts, scales and statistics differ from each other and from the project GWAS, so these numbers should not be directly contrasted with r_g. [Clinical comparison and sources](fi_ffs_same_cohort_comparison.md).

A second analysis using Mak et al.’s released UK Biobank FI and modified Fried physical-phenotype GWAS estimated **r_g = 0.796** (SE 0.012). On a common SNP set, it was 0.798 versus 0.757 for the original pair; the paired post hoc difference was 0.041 (95% CI 0.007–0.075). The studies reuse UK Biobank and differ in scoring, so this is not independent replication and does not identify why they differ. [Full analysis](mak2023_same_cohort_fi_fp_ldsc_review.md).

The UK Biobank FI derivation includes sleeplessness and tiredness; FFS includes exhaustion. The tiredness question in a published FI derivation matches the exhaustion question in Mak's physical phenotype, although the precise coding used in each GWAS requires confirmation. Direct item overlap can contribute to the observed correlations, but the current summary statistics cannot quantify its contribution or distinguish it from other shared influences. [Detailed question and coding comparison](ukb_fi_fp_exhaustion_overlap.md).

## 2. What the insomnia correlations show

| Insomnia GWAS | Insomnia–FI r_g | Insomnia–FFS r_g | FFS minus FI, paired estimate (95% CI) |
|---|---:|---:|---:|
| UK Biobank symptom measure | 0.640 | 0.483 | −0.156 (−0.213 to −0.100) |
| FinnGen clinical insomnia | 0.426 | 0.245 | −0.181 (−0.260 to −0.102) |
| MVP European clinical insomnia | 0.549 | 0.398 | −0.151 (−0.220 to −0.082) |

Paired differences use matched SNP support and jackknife covariance. The same FI and FFS GWAS are reused across the rows; these are comparisons across insomnia sources, not three independent replications of the frailty studies. Different ascertainment, populations, age, sample overlap and phenotype construction remain plausible explanations. Refer to [the source-level analysis](insomnia_frailty_risk_briefing.md).

An independent arithmetic audit of the saved delete-block values reproduces the paired SEs for the later four-frailty-source comparisons; see [paired uncertainty audit](insomnia_frailty_paired_uncertainty_audit.md). That audit uses a different SNP intersection from the table above, so small numerical differences are expected.

A sensitivity analysis using the **UK Biobank-only FI** instead of the UK Biobank-plus-TwinGene meta-analysis gave r_g **0.419** with FinnGen insomnia versus **0.426** for meta-FI on the same SNP support. The paired difference was −0.0074 (95% CI −0.0215 to +0.0066). This suggests that adding TwinGene does not explain much of that specific FI–FFS contrast; it does not test removal of sleep or fatigue items. [Sensitivity details](ukb_fi_finngen_insomnia_sensitivity.md).

In the Mak study’s FI and modified physical-frailty GWAS, the FinnGen insomnia correlations were **0.401 and 0.237**, respectively, on the same 1,129,146 SNPs. The paired FP-minus-FI difference was −0.164 (95% CI −0.216 to −0.112); the analogous original-pair difference on that support was −0.185. These overlapping studies provide a second phenotype-construction sensitivity, not independent replications or proof of the reason for the difference. [Analysis](mak2023_insomnia_fi_fp_comparison.md).

The pattern also appears with MVP insomnia: Mak FI r_g = **0.536** versus modified physical frailty r_g = **0.404**, with a paired difference of −0.131 (95% CI −0.183 to −0.080). This additional matched-support analysis is exploratory; the MVP LDSC intercept is elevated and the studies do not provide independent frailty replications. [Analysis](mak2023_mvp_insomnia_fi_fp_comparison.md).

The UK Biobank insomnia source produces the same direction for the Mak pair: FI r_g = **0.641** versus modified physical frailty r_g = **0.493**; the paired difference is −0.148 (95% CI −0.188 to −0.108). This result has substantial potential participant and phenotype-item overlap. [UKB source analysis and allele QC](mak2023_ukb_insomnia_fi_fp_comparison.md).

The hospital frailty risk score (HFRS) is a separate ICD-code-based phenotype. Its common-support r_g values are approximately **0.640 with FI** and **0.550 with FFS**; the paired FFS–HFRS minus FI–HFRS difference is −0.090 (95% CI −0.142 to −0.038), post hoc. The Finnish HFRS analysis has an elevated univariate LDSC intercept (~1.188), a general European rather than Finnish-specific LD reference and approximate sample-size assumptions, so it deserves particular caution. HFRS is not demonstrably free of sleep or exhaustion overlap. [Methods and diagnostics](frailty_three_way_common_support_ldsc.md).

One observational cohort removed insomnia symptoms from its FI and still reported an association between baseline insomnia and later worsening among prefrail participants; an actigraphy study found an association between sleep fragmentation and later physical frailty. These are complementary clinical observations, not genetic or causal validation. [Risk context](sleep_frailty_observational_risk_context.md).

Longitudinal results also depend on the sleep and frailty definitions: a Japanese cross-lagged study reported associations in both directions, while a three-wave UK Biobank study reported reciprocal associations for its social-frailty measure but not its physical measure. These observational designs do not prove causality or show that the studies directly contradict each other. [Study details and limitations](sleep_frailty_observational_risk_context.md).

In the second frailty pair, HFRS was also more genetically correlated with Mak FI (0.644) than with the modified physical phenotype (0.537). This broad-phenotype context does not explain the insomnia difference; the elevated HFRS intercept and Finnish LD-reference mismatch are material limitations. [Matched comparison](mak2023_hospital_frailty_comparison.md).

## 3. Shared-variant candidates and nearby genes

[Song et al.](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2024.1358996/full) reported 124 insomnia–FI cross-trait candidate variants. In our follow-up, 109 had eligible position and allele-pair matches in both MVP insomnia and FinnGen HFRS inputs under the documented build assumption; 22 had nominal p < 0.05 in both. The six rows below have the strongest values under a max(p) × 124 screening benchmark. Selection occurred after earlier results, so the benchmark is not a substitute for an independent prespecified test. HFRS and MVP insomnia do not establish shared FFS–insomnia causality; LD can tag different underlying signals.

| Candidate | Higher-MVP-insomnia allele / comparator | Nearest gene-body label | Approx. MVP insomnia OR per allele | FFS p |
|---|---|---|---:|---:|
| rs12536335 | A / G | FOXP2 | 1.032 | 0.0017 |
| rs28528997 | A / G | MAD1L1 | 1.036 | 2.3 × 10⁻¹² |
| rs2815764 | G / A | NEGR1 | 1.041 | 1.2 × 10⁻⁵ |
| rs2155292 | G / A | NCAM1 | 1.026 | 0.54 |
| rs817509 | G / A | RP11-384F7.1 | 1.038 | 0.0099 |
| rs1909650 | G / A | PCDH8P1 | 1.024 | 0.097 |

ORs are approximate transformations of rounded source estimates and are relative odds for the indicated allele contrast, not absolute risks. Weak or nonsignificant FFS results do not prove an absence of effect. Gene labels are proximity annotations, not established effector genes. Five of these six rsIDs have multiple alleles in the broader reference; source-specific multiallelic coding and the assumed build require verification. A source-table check supports the FFS A1 effect-allele convention at 37 published lead loci, though it does not individually verify every candidate ([audit](ffs_effect_allele_source_consistency_audit.md)). The FI component estimates for rs12536335, rs28528997 and rs2815764 differ in sign under the inferred TwinGene adapter, which is strongly supported by numerical reconstruction but not author-confirmed. [Full selection and provenance](song124_mvp_hfrs_followup.md); [cohort inconsistency](six_leads_fi_cohort_comparison.md); [all 22 directions and gene labels](22_candidate_risk_direction_review.md).

The wider 22-candidate screen has now also been checked in FinnGen clinical insomnia. Six variants have matching point-estimate directions and nominal p < 0.05 across FI, FFS, three insomnia sources and HFRS; this is post-selected and the sources are not all independent. The two FOXP2-labeled variants also have measurable reference-panel LD. [Full six-source follow-up](22_candidate_six_trait_followup.md).

The same 22 variants were also checked in Mak et al.’s FI and modified physical-frailty GWAS. All estimated directions agreed with the respective original frailty sources; 11 were nominal in both Mak outcomes. Five candidates had positive estimates and nominal p < 0.05 in all eight available FI, physical-frailty, insomnia and HFRS result sets. These are post-selected descriptive findings, not eight independent replications; the two FOXP2-labeled markers are in the same region. [Eight-source matrix and limitations](shared22_eight_source_evidence_review.md).

One hypothesis-generating gene lead is **ADARB1**, whose transcripts contain an annotation for candidate rs2838788. Separate mouse experiments implicate Adarb1 in circadian RNA editing and rhythms. Those experiments did not test the human variant, insomnia or frailty; neither the causal variant nor the effector gene is established in this project. [Gene background and evidence boundaries](shared_candidate_gene_background.md).

Public genomic-SEM frailty factor results offer further context but are not sleep-excluded FI measurements. An exploratory indexed lookup of the 22 markers has uneven coverage across general, multimorbidity and disability factors, whose releases filter heterogeneous SNPs; missing rows are not null results. [Factor definitions and lookup](foote_latent_frailty_gwas_availability.md).

## 4. What the sparse scores do and do not tell us

The [15-variant](unvalidated_15_variant_score_prototype.csv) and [22-variant](unvalidated_22_variant_score_prototype.csv) files contain **three separate weighted sums**, using FI, FFS or UK Biobank insomnia effects. Each sum is a prototype based on marginal GWAS coefficients and selected markers. It is not calibrated to an outcome, trained as a joint score, or evaluated in a cohort with observed outcomes. Coefficients from the three phenotype scales must not be added together.

The two marker lists share only four variants. In 503 European-labeled reference samples without phenotype outcomes, correlations between the 15- and 22-marker sums for the same outcome were 0.346 (FI), 0.198 (FFS) and 0.293 (insomnia). This shows sensitivity to this particular choice of markers and weights; it does not compare predictive accuracy. In the 22-marker FFS sum, rs28528997 accounts for approximately 28% of the sum's sample variance under a covariance allocation. That fraction is not variance in FFS or disease risk explained. [Score comparison](reference_22_vs_15_score_diagnostic.md); [variance diagnostic](reference_22_score_variance_review.md).

Published genome-wide FI, FFS and HFRS scores are useful benchmarks, but their reported results cannot be assigned to our sparse prototypes. For example, [Foote et al.](https://www.nature.com/articles/s41588-025-02269-0) evaluated FI- and Fried-derived scores against observed FI in several cohorts; the score comparison, phenotype overlap and a discrepancy between reported source tables are reviewed [here](multiscore_frailty_prediction_review.md). A credible new model should be judged by whether it improves performance beyond appropriate baseline and existing-score models, not simply by association p-values.

A further cohort-weight sensitivity exposes a major uncertainty: for the same 22 variants in 503 reference individuals, the UKB- and TwinGene-weighted FI sums correlate at **−0.455**, while the UKB- and meta-weighted sums correlate at **0.996**. This relies on an inferred, not author-confirmed, TwinGene sign convention and selected, noisy marginal effects; it is not evidence of opposite biological effects or predictive performance. [Cohort-weight sensitivity](22_fi_cohort_weight_sensitivity.md). A separate frequency audit strongly supports an opposite-allele orientation for the TwinGene column, although it does not confirm beta coding. [Frequency audit](twingene_frequency_orientation_review.md). Removing any one marker leaves the correlation negative; removing four post hoc selected markers changes it to −0.018, showing sensitivity to that subset rather than establishing a rule for exclusion. [Follow-up](22_fi_score_cohort_sensitivity_followup.md).

A concrete source for a sleep-overlap sensitivity has now been identified: Mak et al. report an FI GWAS excluding eight mental-wellbeing items, including insomnia and fatigue. Full summary-statistics availability for this particular modified FI remains unverified, and removing the full category changes more than sleep. [Source audit](mak2023_sleep_excluded_fi_lead.md).

A second lead is Tseng et al.'s diagnosis-and-laboratory FI GWAS. Its published 52-item list does not explicitly name insomnia or fatigue, but broader code groups and indirect overlap require checking; this is a different construct, and a complete summary file was not located in targeted searches. [Item and availability review](tseng_diagnosis_lab_fi_gwas_lead.md).

An exploratory lookup of Tseng's selected FI leads found two nominal insomnia associations in both clinical sources, neither meeting a 31-marker multiplicity benchmark. At APOE rs429358 the insomnia effects point opposite to the published FI effect under an unverified Tseng allele convention; the other marker lies in the MHC. Neither is confirmed shared risk. [Details](tseng_diagnosis_lab_fi_gwas_lead.md).

A separate [22-marker shared-direction count](shared22_score_research_guide.md) is now available as an unweighted research prototype. Its alleles were oriented to matching MVP insomnia and HFRS point-estimate directions after selection; no outcome-based validation, calibration or clinical interpretation is available.

## 5. Next analyses, in order

1. **Separate construct overlap from shared risk using public data.** Seek a freely available genome-wide FI GWAS without the sleep item and compare it against the full FI using matched SNP support and covariance-aware uncertainty. The published Mak category-stripped FI is the closest identified lead, but its full summary statistics have not been located in the inspected release and six additional mental-wellbeing items were removed. The public diagnosis-based HFRS results offer limited triangulation, not a sleep-item deletion. A same-participant deletion GWAS would need participant-level data that the team has said it does not have; it is a future design, not the current public-data plan. [Sensitivity options](public_frailty_sleep_overlap_sensitivity_options.md).
2. **Resolve variant identity and follow up loci.** Confirm genome builds, effect alleles, multiallelic coding and source overlap; use ancestry-appropriate LD, conditional analysis and fine-mapping, then colocalization and functional gene mapping. Nearby-gene labels are starting hypotheses only. The existing FI-by-sleep local genetic-correlation analysis failed its project gates and cannot establish local sharing.
3. **Evaluate prediction if suitable free, public outcome data become available.** Freeze outcomes, marker selection and weights before testing. Use genotype QC and confirm non-overlap and relatedness with discovery and candidate-selection cohorts. Compare baseline clinical models, FI and FFS genome-wide scores, insomnia scores and sparse prototypes using appropriate discrimination, calibration and incremental performance. Assess both FI and FFS outcomes and sleep-excluded FI when measured. [Cohort feasibility and detailed design](independent_outcome_validation_cohort_feasibility.md).
4. **Only after risk is established, explore function.** Use single-cell expression/regulatory datasets to prioritize plausible tissues and cell types, and use AlphaGenome or other models to generate variant-effect hypotheses; evaluate these against fine-mapping and external molecular evidence. Predicted function is not proof of mechanism, and these downstream analyses are future work.

All analyses above are research analyses, not clinical risk estimates. The underlying outputs include exact counts, file hashes, methods, and data limitations; no participant-level outcome validation was performed.

## Public latent-factor follow-up

Foote et al.’s public, filtered latent-frailty summary statistics provide an exploratory seven-factor look at the candidate markers, with unequal marker availability across factors. These factors are not FI or FFS, and the general model includes insomnia. The [factor definitions and seven-factor results](foote_latent_frailty_gwas_availability.md) explain why this is a hypothesis-generating phenotype comparison rather than independent validation or a sleep-item removal test.

A [public-data sensitivity decision table](public_frailty_sleep_overlap_sensitivity_options.md) compares the modified FI leads with the available HFRS and latent-factor results, and identifies what each can and cannot answer about sleep-item overlap.

A [leave-one-region reference analysis](shared22_leave_region_reference_review.md) shows how the proposed unweighted count changes when candidate groups are omitted. It assesses genotype-sum sensitivity only, not predictive value.

For allele-aligned FI and physical-frailty estimates and uncertainty for every candidate, see the [22-marker side-by-side comparison](shared22_mak_fi_fp_comparison.md). The two source outcomes have different scales, so their raw coefficients cannot be compared directly.

For the candidate insomnia odds ratios and the available FinnGen uncertainty intervals, see the [clinical insomnia comparison](shared22_clinical_insomnia_effect_uncertainty.md). The inspected MVP file lacks standard errors for these rows, so MVP intervals were not reconstructed.

The [chromosome 3 gene-mapping audit](chr3_song_candidate_gene_mapping_audit.md) distinguishes Song’s DALRD3 label for rs34290943 from this project’s LAMB2 positional and transcript annotations, and separates cross-trait p-values from marginal results.
