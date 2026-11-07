# Brain6 and sleep–frailty: research overview

A concise [research report](RESEARCH_REPORT.md) brings together the genetic results, mechanism hypotheses, ML findings and unresolved work.

A [public physical-function composite source](component_exclusion/schoeler_physical_function/SOURCE.md) adds a phenotype without explicit sleep or exhaustion items. Its [genome-wide comparison](component_exclusion/schoeler_physical_function/RESULTS.md) is complete, with important intercept and overlap limitations; it is not a matched item-excluded frailty score.

A new [explicit genetic physical-function index](component_exclusion/physical_factor_workaround/INDEX_RESULTS.md) combines walking, grip and activity genetic components and compares the defined index with registry and questionnaire insomnia. It avoids explicit sleep or exhaustion inputs and the residual-covariance identification problem found in the [latent-factor sensitivity](component_exclusion/physical_factor_workaround/RESIDUAL_RESULTS.md). It remains a post-hoc physical-function comparison with source and overlap limitations, rather than a person-level frailty score or matched item-deletion GWAS.

A [summary-statistics reconstruction audit](component_exclusion/identifiability/LINEAR_ROUTE.md) identifies the exact conditions under which raw additive score effects can be reconstructed without participant data, and demonstrates why the available transformed and binary outcomes do not meet those conditions.

The [6 November retained-input LDSC replay](component_exclusion/physical_factor_workaround/VERIFICATION.md) reproduced all nine physical-index input correlations and SEs and all 27 block-delete vectors exactly after the SSD reconnected. This extends the saved-summary replay to the LDSC stage using the same inputs and reference.

A [working manuscript draft](MANUSCRIPT_DRAFT.md) consolidates the research questions, methods, results, limitations and selected references.

For concise, writing-ready drafts, see the [two introductions, prior-study comparison and hypotheses](INTRODUCTIONS_AND_COMPARISON.md).

The [LDSC input audit](data_manifest/README.md) verifies the external files for the two main genome-wide rechecks and describes how to rerun them without overwriting the saved results. A separate [sampled FinnGen source-concordance audit](replication/finngen_adhd_archived_munged/source_concordance/RESULTS.md) strengthens the clinical ADHD association-statistic provenance while leaving the original ingest and cohort overlap unresolved.

This guide summarizes the current evidence and links to the complete [research atlas](RESEARCH_ATLAS.md), which contains the two introduction drafts, prior-study comparison, hypotheses and novelty ledger. The user asked us to disregard the unavailable proposal PDF; the [literature map](LITERATURE_MAP.md) is an independently assembled substitute, not a reconstruction of that document. The original data and larger reference resources remain external; the [selected data manifest](manifest.tsv), [bundle inventory](BUNDLE_INTEGRITY.md) and linked analysis protocols document what is included and what requires reacquisition.

## Brain6: what the evidence supports

Published 2026 work already reported shared insomnia–ADHD genetics and the rs62264764 region. Brain6's re-execution gives a genome-wide genetic correlation **r_g = 0.4033, SE = 0.0318** on the same saved inputs. This is a consistency check, not independent cohort replication, a causal effect, or evidence that individual loci share a causal variant. See the [recheck](global_audit/insomnia_adhd_rerun/RESULTS.md) and [literature comparison](RESEARCH_ATLAS.md#8-prior-study-comparison).

At **rs62264764**, AlphaGenome predicts an allele-dependent frontal DNase signal in a 501-bp window; it is a model prediction, not a measured accessibility change. LSAMP and LINC03051 are overlapping candidate genes, not established mediators. Existing brain-cell molecular associations involve LINC03051, and the published LaDDR study found brain RNA-feature–insomnia associations; distinct tissues, donors and RNA phenotypes limit direct comparison. One-signal molecular colocalization pilots have not established a shared disease–gene mechanism. The [gene and pathway review](pathway_hypotheses/GENE_PATHWAY_REVIEW.md), [brain xQTL review](variant_context/rs62264764_brain_xqtl/RESULTS.md) and [LaDDR review](variant_context/linc03051_laddr_2026/RESULTS.md) keep these evidence types separate. A [same-marker FinnGen sensitivity](variant_context/rs62264764_brain_xqtl/excitatory_coloc/finngen_sensitivity/RESULTS.md) is dominated by the expression-only hypothesis and does not resolve whether the disease and molecular signals share a causal variant.

The chromosome 5 shared region has encouraging **single-signal** insomnia–ADHD colocalization, PP4 0.991 at p12=10⁻⁵ and 0.917 at p12=10⁻⁶. It is previously reported, model-dependent, and does not establish the gene or repair the failed local-correlation gates. The [six-region locus atlas](locus_atlas/LOCUS_ATLAS.md) also distinguishes the separate chromosome 5 ADHD-prioritized region from this shared region. No independently verified two-trait locus replication is established; the [follow-up comparison](replication/six_lead_joint/RESULTS.md) explains the promising but post-selected related-phenotype evidence.

A [public clinical ADHD follow-up](replication/finngen_adhd_raw_lookup/RESULTS.md) finds a same-direction association for the chromosome 5 lead rs2431108 in FinnGen R13. It is a retrospective, selected-locus comparison with cohort-overlap and source-verification limits; it does not establish a shared causal variant or repair local QC. A [separate MVP-insomnia and FinnGen-ADHD comparison](replication/mvp_finngen_clinical/RESULTS.md) finds no lead crossing its descriptive six-lead benchmark; the comparison is post-selected and does not establish independent replication.

A [sleep-side effect sensitivity](replication/mvp_finngen_clinical/sleep_effect_sensitivity/RESULTS.md) compares MVP and FinnGen insomnia estimates at the six leads; reconstructed uncertainty and unknown cohort dependence prevent a calibrated heterogeneity conclusion. A [source diagnostic](replication/mvp_finngen_clinical/sleep_effect_sensitivity/source_ci_audit/RESULTS.md) found the CI–p mismatch in a broader, non-random file segment; its cause is unverified. A [check of the original MVP release](replication/mvp_finngen_clinical/sleep_effect_sensitivity/source_ci_audit/original_release/RESULTS.md) confirms the selected values and that the original file has no standard-error column.

An [Estonian clinical-insomnia preprint lookup](replication/estonian_insomnia_preprint/RESULTS.md) finds two of the 425 candidate variants in a selected follow-up table, neither with nominal insomnia support. None of the six original leads occurs in that table, which is not complete regional summary data and cannot settle their replication status.

Potential **neuronal adhesion, circuit development and excitatory/inhibitory signaling** pathways are hypotheses informed by prior work, not pathways demonstrated to be mediated by these variants. No clinical prediction or causal gene assignment follows from the current evidence.

## Sleep and frailty: what differs

A [sleep–ADHD–frailty mechanism crosswalk](biological_mechanisms/three_trait_bridge/RESULTS.md) joins the existing reporter candidates and genetic follow-ups, identifying hypotheses and separating frailty-item overlap from independent evidence.

An [exploratory NHANES phenotype pilot](component_exclusion/nhanes_public/phenotype_pilot/RESULTS.md) estimates same-person correlations using approximate FI and adapted Fried scores, including versions with sleep or the walking-based exhaustion proxy removed. A [same-sample mobility-skip sensitivity](component_exclusion/nhanes_public/phenotype_pilot/visit_code_audit/mobility_same_sample/RESULTS.md) shows that the working FI–adapted-Fried correlation changes substantially under two hypothetical scores for skipped walking items. The scores and missing-data assumptions need validation; this is not an item-excluded genetic analysis. A [later public NHANES availability check](component_exclusion/nhanes_public/prepandemic_availability/RESULTS.md) identifies a two-item physical sensitivity without exhaustion and a weight-loss intentionality skip that needs explicit handling. A [later-cycle observational pilot](component_exclusion/nhanes_public/prepandemic_availability/physical_sleep_pilot/RESULTS.md) compares this partial score with and without exhaustion; both sleep-duration correlations are imprecise.

A [clinical-insomnia sensitivity](replication/mvp_frailty_genomewide/RESULTS.md) also finds positive genome-wide correlations with FI, FFS and HFRS, with elevated intercepts and phenotype overlap limiting interpretation.

For individual physical functions, a [clinical-insomnia sensitivity](component_exclusion/mvp_physical_components/RESULTS.md) found negative genome-wide correlations with walking pace and grip strength. These are not modified frailty scores; phenotype coding, intercepts and unknown cohort overlap limit interpretation.

FI counts a broad collection of health deficits; the five-item Fried score (FFS) focuses on weight loss, exhaustion, weakness, slowness and low activity; hospital frailty risk score (HFRS) uses hospital diagnosis codes. They are different phenotypes measured in different settings. The FI used here includes insomnia and tiredness, while FFS includes exhaustion, making direct item overlap relevant to sleep comparisons.

The same-data FI–FFS LDSC rerun gives **r_g = 0.7534, SE = 0.0199**. This is a genetic-effect correlation, not a correlation between individual scores. In a separate published UK Biobank 3×3 category table, the reconstructed within-person ordinal Spearman correlation is **0.2845**; it is not the continuous-score correlation for the GWAS participants. Definitions, samples and estimands differ, so these numbers should not be compared as if they measured the same quantity. See the [genetic recheck](global_audit/fi_ffs_rerun/RESULTS.md) and [phenotypic comparison](phenotype_comparison/PHENOTYPIC_VS_GENETIC.md).

On common SNP support, the original insomnia source gives **r_g 0.6404 with FI, 0.4840 with FFS and 0.4355 with HFRS**. Post-hoc paired contrasts and their limitations are in the [three-frailty comparison](global_audit/brain6_frailty_covariance/three_frailty_common_support/RESULTS.md). Differences cannot be attributed to physical frailty alone: phenotypes, cohorts and sleep-related content differ. A [published-context review](component_exclusion/zenodo_published_context/RESULTS.md) clarifies that the public component GWAS, an observational exhaustion score and the published composite GWAS are distinct resources. A [public sarcopenia-factor screen](component_exclusion/sarcopenia_public_candidate/RESULTS.md) adds selected-marker evidence for a related physical construct, but its input traits and sample-size semantics need documentation. A [four-trait feasibility check](component_exclusion/sarcopenia_public_candidate/ldsc_feasibility/RESULTS.md) obtained numerical LDSC estimates, but these remain unvalidated and do not constitute an item-excluded frailty GWAS. The requested matched FI-minus-sleep and FFS-minus-exhaustion GWAS have **not** been estimated; the [construction protocol](component_exclusion/CONSTRUCTION_PROTOCOL.md) and [observational precedent](component_exclusion/zhu_observational/REVIEW.md) explain what is and is not available.

## Local results, multiple testing and ML

A [candidate-level experimental and AlphaGenome crosswalk](biological_mechanisms/cross_assay_prediction/RESULTS.md) joins 27 published neural-progenitor reporter results to saved predictions and highlights window sensitivity and unresolved gene targets.

The frozen local-correlation quality gates failed. A targeted chromosome 3 replay identified sensitivity to SNP support, LD components and reference choice, plus numerical defects in the extracted approximate UKB LD matrix. An alternative small 1000 Genomes reference does not validate a replacement. See the [local audit](local_evidence/LOCAL_FRAILTY_AUDIT.md) and [reference comparison](local_evidence/song_lava_receipt/cross_trait/marker_overlap/replay/reference_comparison/RESULTS.md). Missing or failed local estimates are unknown, not zero.

A [one-region HDL-L pilot](local_evidence/hdll_pilot/RESULTS.md) returned positive insomnia–frailty estimates under hypothetical sample-overlap settings. They are uncorrected, post-selected and based on the same data; they do not repair the failed LAVA gates or establish replication.

A [public neural-progenitor MPRA](variant_context/crossdisorder_mpra/published_table/RESULTS.md) includes 42 of 425 ranked variants in its count matrices, but not rs62264764. Its recovered published workbook matches 27 of those 42 and labels rs325485 and rs6874138 as allele-sensitive active elements. Their moderate reference LD with the shared-region lead does not establish a causal disease variant. An ABC link from rs6874138 to NUDT12 motivates a circadian-related hypothesis, but the [available astrocyte expression associations](variant_context/crossdisorder_mpra/published_table/nudt12_brain_qtl/RESULTS.md) are weak and no target mechanism is established.

In the separate ADHD-prioritized chromosome 5 region, [rs62369151](variant_context/crossdisorder_mpra/published_table/CROSS_ASSAY_RESULTS.md) has allele-effect signals in two selected reporter assays and reference r²≈0.83 with rs4916723. Both linked variants show [brain-cell expression associations](variant_context/crossdisorder_mpra/published_table/rs62369151_brain_qtl/RESULTS.md) with LINC02060 and TMEM161B-AS1, while a prior single-signal analysis favored distinct dominant ADHD and expression signals. These assays and associations do not establish a shared insomnia–ADHD mechanism.

An [older-adult low-grip-strength sensitivity](component_exclusion/older_weakness/RESULTS.md) found positive genetic correlations with questionnaire and registry insomnia (0.169 and 0.171 respectively). Grip strength is one physical component, not an exhaustion-excluded frailty score; cohort overlap, approximate sample size, reference fit and the broader testing burden limit interpretation.

For a physical-function comparison that does not itself include a sleep item, a [published walking-pace GWAS](component_exclusion/walking_pace_published/RESULTS.md) reported genetic correlations with insomnia and ADHD. It does not replace the still-unavailable matched sleep-item-excluded frailty analysis. A [falls-GWAS sensitivity](component_exclusion/falls_sleep_ldsc/RESULTS.md) also finds positive genome-wide correlations with two insomnia sources; falls are a distinct outcome, not an exhaustion-excluded frailty score.

A [new exploratory registry-insomnia comparison](component_exclusion/physical_function_ldsc/RESULTS.md) with UK Biobank walking pace and left grip strength found only grip significant after the two-test correction. Both are individual physical-function traits, not modified frailty composites.

On shared markers, a [two-source sensitivity](component_exclusion/physical_function_ldsc/insomnia_source_sensitivity/RESULTS.md) found a materially more negative walking-pace correlation for questionnaire insomnia than registry insomnia. Differences in cohorts, measurement and reference fit prevent attributing this to one cause.

An [expanded five-outcome comparison](component_exclusion/physical_function_ldsc/insomnia_source_sensitivity/frailty_extension/RESULTS.md) shows that FI and FFS correlations are higher for questionnaire insomnia, while HFRS is higher for registry insomnia. These are source-sensitive, post-hoc associations with different cohort overlaps and phenotype definitions.

The authors' newly located public UK Biobank HFRS data enabled a [two-cohort hospital-frailty comparison](global_audit/ukb_hfrs_public/RESULTS.md). The results vary across both HFRS and insomnia sources; an unresolved sample-count discrepancy in the UKB release limits interpretation.

A [dementia-weight exclusion sensitivity](global_audit/ukb_hfrs_public/dementia_sensitivity/RESULTS.md) gives similar original and modified UKB HFRS–insomnia genetic correlations. It does not answer whether associations persist after removing sleep or exhaustion from FI or FFS.

An [HFRS screen of the 425 Brain6 candidates](variant_context/brain6_ukb_hfrs/RESULTS.md) finds a strong chromosome 7 signal including rs1476535. It is a selected same-cohort screen, not independent sleep–ADHD–frailty replication.

A [MVP-insomnia and Finnish-HFRS candidate comparison](variant_context/brain6_mvp_finngen_hfrs/RESULTS.md) adds matching sleep–frailty evidence for rs62264764 and rs1476535, but their separate childhood ADHD-symptom associations are weak. The selected markers do not establish three-trait replication or a common causal variant.

The [expanded childhood-symptom check](variant_context/brain6_three_trait_candidate/RESULTS.md) found no crossing of the exploratory three-phenotype benchmark among the 425 selected variants. This does not rule out shared biology or an effect below the available power.

A [chromosome 7 HFRS–insomnia colocalization pilot](frailty_coloc/ukb_hfrs_chr7/RESULTS.md) favors distinct dominant signals under its single-signal model and tested priors. It cannot exclude secondary sharing or identify a causal gene.

An [insomnia-source sensitivity on the same variants](frailty_coloc/ukb_hfrs_chr7/insomnia_source_sensitivity/RESULTS.md) gives a different and prior-sensitive questionnaire result. Questionnaire insomnia and HFRS both use UK Biobank, while the coloc model assumes independent statistics.

A [correlated-statistic stress test](frailty_coloc/ukb_hfrs_chr7/insomnia_source_sensitivity/overlap_stress/RESULTS.md) shows that a different method's questionnaire cluster disappears at tested positive hypothetical correlations. The actual overlap correlation is not known.

At the original chromosome 5 region, a [three-trait HFRS pilot](frailty_coloc/chr5_hfrs_three_trait/RESULTS.md) retains strong model-based insomnia–ADHD sharing but does not return HFRS in the cluster. Prior, overlap and single-signal limitations remain, and the failed local gate is unchanged.

Near rs62264764, a [chromosome 3 three-trait pilot](frailty_coloc/chr3_hfrs_three_trait/RESULTS.md) returns a tentative cluster under independence, but HFRS drops out at tested positive hypothetical dependence. The true dependence is unknown and no causal gene or three-trait replication is established.

The conservative [three-frailty SNP conjunction screen](conjunction_screen/three_frailty_union/RESULTS.md) produced no corrected marker–pair rejections on its restricted supports. This does not prove an absence of shared genetics. The exploratory [PLACO+ diagnostic](placo_diagnostic/completed_audit/RESULTS.md) has 82 numerical crossings, but null-tail stress testing raises a calibration concern. A separate [conjunction-FDR test](conditional_enrichment/conjfdr_readiness/r_implementation/synthetic_ld_overlap/RESULTS.md) also raises calibration concerns; no real-data conjFDR finding is claimed. No analysis here resolves a shared causal variant solely by crossing a marginal threshold.

For context, [PLACO+ random-effect null calculations](placo_diagnostic/tail_stress/model_null/RESULTS.md) are conservative under the specific distributions stipulated by the method and the fitted parameters, while the earlier fixed-effect stress cases can be anti-conservative. The real-data null distribution remains unverified, so this model check does not promote the crossings.

The [observed marginal Z-score tails](placo_diagnostic/tail_stress/empirical_margins/RESULTS.md) are heavier than a single fitted normal reference. Since they include true signals, LD and selected rows, that descriptive comparison does not measure error rates or test PLACO+'s component null directly.

The [trait embedding](multiscale_embedding/REPORT.md) is a research prototype: its tested masked-cell performance did not beat simple baselines, and it is not an individual risk model. A better next experiment would use appropriately licensed genetic and molecular features, hold out whole LD blocks and independent source studies, and compare retrieval against association and fine-mapping baselines. AlphaGenome output terms must be checked before using any model outputs for training.

An [MVP-specific feature-transfer experiment](multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/RESULTS.md) supplies six clinical-insomnia–disease genetic correlations to the existing model. Its Fried-score prediction was 0.607 versus the observed correlation 0.400; reused GWAS sources and already-inspected targets make this an exploratory stress test, not independent validation. A [fixed-model input sensitivity](multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/bounded_model/mvp_feature_transfer/input_stress/RESULTS.md) shows that six simultaneous one-standard-error input changes leave the prediction above that target point estimate; the range is not a confidence interval.

In a later [held-out-group ridge model](multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/RESULTS.md), [feature ablations](multiscale_embedding/hfrs_extension/heldout_sleep/group_holdout/feature_ablation/RESULTS.md) show that the small evaluation changes with input choice. This does not establish feature importance or independent generalization; the model requires existing sleep–brain correlation estimates as inputs.

## Highest-priority remaining work

1. Obtain or generate authorized same-person, sleep/exhaustion-excluded FI and FFS GWAS, with matched sample and model; re-estimate paired genetic comparisons.
2. Obtain study-appropriate LD and overlap information, freeze a local-analysis and calibration plan, and rerun quality gates before reporting local correlations.
3. Validate a two-trait locus in independently ascertained sleep and ADHD or frailty data; follow with multi-signal molecular colocalization and targeted perturbation.
4. Calibrate cross-trait methods under relevant composite-null, LD and sample-overlap conditions before interpreting PLACO+ or conjFDR results.
5. Test proposed gene/pathway and ML hypotheses against independent data and simple baselines, respecting source-specific terms.

The [scope audit](SCOPE_AUDIT.md) is the detailed checklist of completed work, unresolved evidence and links.

For a concise mechanistic next-step plan, see [rs62264764 evidence and follow-up](variant_context/RS62264764_FOLLOWUP_PLAN.md).

The [external dependency inventory](data_manifest/EXTERNAL_DEPENDENCIES.md) lists which key source files are mounted and hash-matched, and what remains unavailable or unverified.
