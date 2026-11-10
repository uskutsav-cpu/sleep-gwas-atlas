> Update: The original Atkins supplementary FI item tables were located; see the [source crosswalk](fi_ffs_component_crosswalk.md) for the stated coding and remaining code-verification caveat.

# Insomnia–frailty risk: FI versus Fried score and an evidence inventory

Research note, 27 October 2026. This is a read-only synthesis and a check of existing project outputs, not a newly run GWAS, a systematic review, or a validated prediction model.

## Main finding and correction

An October project package contains a covariance-aware comparison that was absent from the older September snapshot I consulted earlier. My earlier statement that no sleep–FFS analysis had been run was therefore incorrect. The three comparisons below reuse the same FI and FFS frailty GWAS and use three different insomnia sources. In all three, insomnia's estimated genetic correlation with the particular FI is higher than with the particular FFS. This is a comparison of genome-wide correlation estimates, not a comparison of clinical incidence or individual risks.

| Insomnia source | FI r_g | FFS r_g | FFS − FI | SE of difference | 95% CI | Raw p | p × 396 (capped at 1) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Jansen UKB symptoms | 0.640 | 0.483 | -0.156 | 0.029 | [-0.213, -0.100] | 6.2e-08 | 2.46e-05 |
| FinnGen R13 clinical insomnia | 0.426 | 0.245 | -0.181 | 0.040 | [-0.260, -0.102] | 7.2e-06 | 0.00286 |
| MVP European clinical insomnia | 0.549 | 0.398 | -0.151 | 0.035 | [-0.220, -0.082] | 1.6e-05 | 0.0065 |

The paired difference is computed with covariance from matched 200-block jackknife estimates, so its SE is not calculated by treating the correlations as independent. An October 10 audit independently recomputed the estimates and uncertainty from saved native outputs; that is not an end-to-end rerun from the raw inputs. The original symptom comparison used the 396-test family; the clinical-source families were 16 and 4 respectively, and the last column is a conservative 396-test overlay. Different insomnia sources share the same frailty sources: they are sleep-source sensitivities, not independent frailty replication, and the project has not tested whether their differences from one another are significant. The first comparison used its historical paired-uncertainty audit; the two clinical comparisons also passed point and 200-delete joint three-trait matrix checks. Do not infer a novel biological mechanism from the contrasts.

Source artifacts: `.../2026-10-10/new-chat/outputs/strengthening/tables/THREE_SOURCE_INSOMNIA_CONTRASTS.tsv`; `.../2026-10-10/sleep-frailty-v8-isolated/frailty_paper/analysis/v9_phase1_recheck_2026-10-10.md`; and `.../strengthening/reviews/statistical_genetics/INDEPENDENT_STATISTICAL_REVIEW.md`. Paths are identifiers of the existing local evidence, not new results.

## What FI and FFS actually measure

The FI is a proportion of accumulated health deficits. The specific Atkins source meta-analyzed 164,610 UK Biobank participants (49 deficits) and 10,616 TwinGene participants (44 deficits), with a quantile-normalized FI; the study samples were of European descent. The Fried score in the Ye GWAS is a 0–5 score for five physical domains: weight loss, exhaustion, low activity, slow gait, and weakness. It is a different measurement, not a clean replacement for FI. [Atkins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/); [Ye et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/).

The UK Biobank FI's items include sleeplessness/insomnia (field 1200) and tiredness/lethargy (field 2080). The Atkins supplementary TwinGene table lists 44 items with no explicit sleep or fatigue question; the executed code has not been independently audited. The full meta-analysis is therefore not a uniform 49-item sleep-containing measure. The public coding sources also differ for one fatigue response, and the exact executed original scoring and Ye's full field/cutoff mapping remain unverified. FFS exhaustion and FI fatigue also overlap conceptually. Participant intersections, including those shared by the frailty GWAS and original UKB insomnia study, remain unknown.

**Why this matters:** part of the FI–FFS contrast may reflect phenotype content or other differences in sample, measurement, GWAS model, ascertainment or precision. Its causes cannot be isolated by comparing the existing correlations. A direct same-participant FI analysis with the sleep item removed, then with sleep and fatigue removed, compared with a matched FFS, remains unperformed. Marginal GWAS results cannot reconstruct this exact deletion, particularly because the primary FI was quantile-normalized. The October protocol treats that as a future authorized study, not an existing result.

An illustrative *direct observed-score* comparison in 909 community-dwelling adults aged 65 or older reported a correlation of 0.76 between its FI and Fried phenotype, but frailty classification agreement was only kappa 0.38. Those are cohort-specific clinical measurements, not the present project's genetic correlation, and the specific score constructions differ. [Thompson et al.](https://digital.library.adelaide.edu.au/items/d8c5819d-7317-47e2-a4e6-d57933328854)

## Interpret the word “risk” carefully

* **FI and FFS scores** measure observed frailty; neither is an individual's genetic risk score.
* **Genetic correlation** describes the normalized genome-wide genetic covariance of two traits. It is not an odds ratio, absolute risk, correlation of two measured clinical scores, or probability of future frailty. The FI-vs-FFS contrast above is not the genetic or phenotypic correlation between FI and FFS themselves.
* **An associated variant** is a marker with statistical association. A shared genomic region is not proof that the same causal variant affects both traits; the effect allele is not necessarily a risk-increasing allele for both.
* **A polygenic score** sums allele dosages weighted by specified effect estimates. It requires defined alleles/build, genotype data, a prespecified model, external validation and calibration. A later technical check used public reference genotypes but no frailty or insomnia outcomes; no personal risk or validated 15–30-variant score can be calculated here.

## Nineteen quantitative risk-related signals already in the project

The companion CSV contains 19 distinct *genetic correlation comparisons*: twelve sleep traits with FI and insomnia with seven latent-frailty factors. They are a reproducible evidence inventory, not 19 risk loci, causal factors, predictors or items of a risk score. The 396-family and 84-family q values come from different frozen analysis families and are not interchangeable. Direction follows original phenotype coding; a negative sign is not automatically a protective effect of changing the trait. No correction is recomputed on this displayed subset.

### Sleep trait versus FI — inherited 396-test BH family

| Sleep trait | Frailty measure | r_g | SE | q |
|---|---|---:|---:|---:|
| insomnia | FI | +0.6405 | 0.0235 | 1.71e-160 |
| sleepdur | FI | -0.2220 | 0.0275 | 7.71e-15 |
| shortsleep | FI | +0.4543 | 0.0272 | 1.57e-60 |
| longsleep | FI | +0.2925 | 0.0375 | 6.38e-14 |
| chronotype | FI | -0.0378 | 0.0253 | 0.243 |
| sleepiness | FI | +0.2779 | 0.0269 | 1.32e-23 |
| napping | FI | +0.2506 | 0.0228 | 1.37e-26 |
| snoring | FI | +0.2516 | 0.0256 | 1.72e-21 |
| sleep_apnea | FI | +0.4687 | 0.0321 | 2.32e-46 |
| sleep_efficiency | FI | -0.1520 | 0.0365 | 0.00012 |
| accel_sleep_duration | FI | -0.0479 | 0.0320 | 0.242 |
| sleep_timing | FI | +0.0668 | 0.0411 | 0.197 |

Nine of the 12 FI rows pass q≤0.05. The original insomnia value 0.6405 is the atlas pair estimate; the comparison table uses a matched SNP support across FI and FFS and has 0.6396. They answer closely related but not identical computational questions.

### Insomnia versus latent frailty — separate 84-test BH family, sensitivity only

| Sleep trait | Latent measure | r_g | SE | q |
|---|---|---:|---:|---:|
| insomnia | frailty_general | +0.6658 | 0.0256 | 1.73e-147 |
| insomnia | frailty_factor_6 | +0.2101 | 0.0364 | 3.45e-08 |
| insomnia | frailty_factor_3 | +0.1276 | 0.0279 | 1.47e-05 |
| insomnia | frailty_factor_2 | +0.1326 | 0.0306 | 4.32e-05 |
| insomnia | frailty_factor_5 | +0.0728 | 0.0275 | 0.0153 |
| insomnia | frailty_factor_1 | +0.0644 | 0.0389 | 0.15 |
| insomnia | frailty_factor_4 | +0.0331 | 0.0412 | 0.506 |

The general latent factor includes insomnia itself among its underlying indicators, so its high correlation is part–whole dependent and is not independent replication. All seven comparisons are sensitivity analyses; overlap is unresolved. The source authors' tentative descriptions are Factor 1 limited social support; Factor 2 unhealthy lifestyle; Factor 3 multimorbidity; Factor 4 metabolic problems; Factor 5 poorer cognition; Factor 6 disability. These are labels for residual statistical factors, not definitive mechanisms. [Foote et al.](https://www.nature.com/articles/s41588-025-02269-0) Two do not pass the 84-test q≤0.05 threshold. [Project manuscript](https://github.com/uskutsav-cpu/sleep-gwas-atlas/blob/main/frailty_paper/paper/manuscript.md).

## Published candidate shared region and gene links

Song et al. reported a 3p21.31 local association and two variants in their cross-trait analysis: rs34290943, mapped in their paper to **DALRD3**, and rs10865954, mapped to the **MIR191** locus. Their published multi-trait p-values should not be reinterpreted as independent validation or as a calibrated risk score; the paper reports 124 cross-trait SNPs, but we have not verified all underlying supplement rows, LD independence, effect directions, builds and gene mappings here. Their unconstrained global r_g estimate was 1.14, outside the theoretical range, whereas a constrained-intercept result was 0.61; this is an important model-sensitivity warning. Use these two as published follow-up candidates, not genes established by this project's failed local analysis. [Song et al. 2024](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2024.1358996/full).

The project's FI×sleep LAVA run completed all 29,940 slots, but its post-remount audit recorded failure of all 12 predeclared 1% trait gates, including 5,125 negative-variance failures, 252 no-reference-SNP failures and seven other failures. These count process-failure events in the frozen analysis, not independent causal loci. No local findings were promoted; this does not prove absence of local sharing. Audit: `.../sleep-frailty-v8-isolated/frailty_paper/analysis/lava_intervention_recheck_2026-09-29.md`.

A later paper by Ye et al. reports an insomnia–FFS correlation of 0.824, but the inspected methods and supplementary audit identify **GNOVA** for that correlation, not the paired LDSC procedure above; the exact insomnia input/configuration was not fully reconciled. Do not use the 0.824 as if it were a direct replication or a like-for-like contradiction of the matched LDSC FFS estimate. [Ye et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/)

## What a genuine risk-score study would require

A published comparator exists: [PGS005229](https://www.pgscatalog.org/score/PGS005229/) is an FFS-derived SBayesR score with 2,761,868 variants based on the European-ancestry Ye UKB GWAS. Foote et al. evaluated it against FI in three separate European-ancestry cohorts; the catalog reports cohort- and model-specific measures. That makes it a useful *frailty-score comparator*, not an insomnia–frailty combined score. R² and effect metrics depend on model, cohort and outcome, and do not yield a person's absolute risk without an appropriate calibrated model.

To develop a 15–30-variant research score, first define the predicted outcome and time horizon; select and freeze variants and weights in a development set using harmonized independent discovery GWAS, LD-aware selection and the chosen endpoint; evaluate it in genuinely independent genotyped cohorts against baseline and genome-wide score comparators; and report ancestry, sample overlap, calibration, uncertainty, discrimination and incremental performance. A fixed 15–30-variant size is a research constraint, not evidence of predictive adequacy. Selection from published lead variants without validated weights, or from genes nearest to them, does not create a valid score. Aggregate GWAS statistics alone cannot calculate an individual's dosage-weighted score or validate its prediction performance.

## Priority work and evidence gaps

1. Verify and preserve the existing paired contrast and its source provenance; reconcile the exact Ye GWAS conventions and FI item coding. The reported comparison is already available; do not rerun merely to manufacture novelty.
2. Obtain appropriate authorization and same-person phenotype/genotype data for the sleep-item and fatigue-item deletion study; require a genuinely independent matched frailty cohort for validation. No such access or sufficient power is established here.
3. For locus discovery, define a separately justified analysis and frozen multiplicity family before results; use harmonized dense summaries and ancestry-matched LD, assess overlap, failure rates, local h², and fine-mapping diagnostics. Do not interpret the failed frozen LAVA family as a passing screen.
4. Build a variant evidence table with genomic build, effect alleles and directions for each trait, association and LD context, fine-mapping uncertainty, and explicit locus/variant-level evidence. Add nearest genes as distance annotations only, then QTL/colocalization and regulatory links if source and QC gates pass.
5. For cell types, use a prespecified, source-licensed single-cell expression/chromatin atlas and appropriate stratified LDSC or related analyses with baseline annotations, multiple-testing control and sensitivity checks. Enrichment in both single traits is not proof of the same shared causal variant or mediation.
6. Compare a published frailty PGS and a separate insomnia PGS in eligible independent data, then test any proposed combined or sparse score against them. No combined score or target-cohort validation has been run.

This note does not infer causal direction, a protective intervention, target gene, cell type, or individualized clinical risk from genetic correlations.

## Subsequent source-audit update

A direct FI–FFS correlation was subsequently located in saved auxiliary LDSC outputs: r_g=0.7570 (SE 0.0212) and r_g=0.7548 (SE 0.0200) on two matched SNP sets, using the same frailty GWAS. The point estimates and block-jackknife SEs were independently recalculated from the saved native files. The full explanation, published FI/FFS PGS comparison and locus-candidate source audit are in the accompanying `fi_ffs_risk_followup.md`. The Song et al. supplement was also recovered after this note was first written: its 124 cross-trait candidates and a 30-row ranking are now supplied as separate CSVs. These are not proven shared causal variants or a risk score.
