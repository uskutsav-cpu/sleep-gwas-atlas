# Insomnia–frailty: public experimental follow-up

## What the data show

The 488-allele queue has a checksum-verified AlphaGenome response for every variant. Its variant-level index, `stage5_variant_phenotype_experimental_index.csv`, joins the six-study disease association matrix to response coverage and the public reporter-assay matches. Detailed score rows and candidate tissues are in `stage5_live_tissue_evidence_atlas.csv.gz`; all raw outputs and metadata are in `stage5_live_raw_scores.zip`. QTL records and event information remain in `stage5_variant_gene_tissue_phenotype_integrated.csv`. These linked tables are a research atlas, not a clinical risk predictor.

Two public allele-specific reporter sources provide experimental context:

1. Ho et al. 2025's primary CD4 T-cell and Jurkat results each exactly match 88 queue alleles. One primary-cell allele, rs9814873 (chr3:49416679:A>G), is classified by the source as an expression-modulating variant; the same locus has an INKA1 muscle eQTL association and frailty associations in this project. This experiment does not establish INKA1 as the reporter's target. Insomnia association differs between the two examined sources, so a shared disease effect is unproven.
2. Siraj et al. 2026's selected-library analysis matches 42 reference-checked variants across available cell lines; five meet the authors' activity-and-skew rules in SK-N-SH. Those five may also recur in other lines or source subsets. The SMAD5 alternative chr5:136207208:G>A has a strong raw SK-N-SH reporter result but is absent from the source's selected analysis subset and must remain exploratory. rs323509 was tested in SK-N-SH and did not pass the activity-and-skew screen; this does not assess its testis sQTL.

See `stage5_independent_reporter_evidence_review.md`, `stage5_siraj2026_neural_reporter_review.md` and `stage5_siraj2026_allcell_reporter_review.md` for exact figures, sources, genome-build matching, prior colocalization findings and reproducibility. The source-specific coverage and track tables retain missing results and experiment details rather than converting untested variants into negatives.

## What remains unresolved

These reporter studies are separate experimental measurements of engineered constructs. Their complete independence from AlphaGenome training and benchmarking has not been established; the `ALL_FOLDS` model is not a held-out genomic-region test. Reporter RNA activity and predicted DNase or other model outputs have different scales and biology, and a count-based cross-check supports ALT-versus-REF orientation for five selected SK-N-SH candidates, but not globally; see `stage5_siraj2026_reporter_direction_review.md`. Some source variants are repeated across libraries, tissues or traits; this does not create independent replications. SK-N-SH is a neuroblastoma cell line, not normal brain.

There are no exact reporter matches for the original FOXP2 insertion or original SMAD5 marker in the two queried reporter sources. Three splice modalities were not returned for 221 variants, and no event/tissue match passed the project's narrow splice-junction screen. Public nonmatches do not establish zero effects, and the selected comparisons do not establish gene mediation or shared causal variants. Previous FOXP2 and SMAD5 colocalization sensitivity analyses carry incomplete coverage, reference-LD, sample-overlap and model-assumption limitations.

The next decisive evidence would be a publicly available or newly generated endogenous, allele-specific assay in a relevant cell state, with a verified variant–gene link, together with well-powered colocalization using ancestry- and study-matched LD and a prespecified, independent phenotype cohort. Free public summary data can support further triangulation, but cannot substitute for an assay or cohort that does not exist or is inaccessible. No independent clinical prediction, causal mechanism or comprehensive validation is claimed.

A provenance audit found a cell-state label discrepancy for two model CAGE strands: the API labels them naïve CD4, while the matching published FANTOM source-file description says memory conventional CD4. Their equivalence to the current API model version needs confirmation. The affected 176 primary-T-cell reporter/model comparison rows are flagged; see `stage5_alphagenome_track_independence_audit.md`. This affects CAGE interpretation, not the separate SK-N-SH DNase direction check.

The count-based reporter direction check was extended to all 284 selected cell-line result groups: nine have a different sign from the study's fitted skew. Because the simple ratio and fitted GLM are different statistics, this is an unresolved consistency issue, not proof of an allele reversal. See `stage5_siraj2026_allcell_orientation_review.md`; signed interpretation should not be generalized from the five highlighted candidates.

Siraj et al.'s published methods define their fitted allelic effect as ALT versus REF. The resulting descriptive comparison with matched-label AlphaGenome DNase tracks is in `stage5_siraj2026_fitted_direction_review.md`; it is not an independent model-accuracy benchmark. The separate Ho-study reporter direction remains unresolved.

A later public FinnGen immune multiome query added portal-returned molecular records for 460 of the 488 variants and position-level no-result responses for 28. For the original SMAD5 marker and a fine-mapped alternative, the top naïve-CD4 eQTL feature is SMAD5 and the top caQTL is a distant peak; neither establishes a shared causal effect. This provides a separately generated molecular dataset, but FinnGen donor overlap with phenotype GWAS and the current model's training provenance remain unresolved. See `stage5_cascade_multiome_readme.md` and `stage5_cascade_smad5_priority_review.md`; the augmented variant index is `stage5_variant_phenotype_experimental_multiome_index.csv`.

The public CASCADE colocalization portal was additionally checked in ±250 kb windows around all 23 exploratory markers. No returned phenotype code explicitly contained “insomn”, “sleep” or “frail”, but the portal displays only a filtered subset and the phenotype-code directory was unavailable, so this is not a negative colocalization result. The SMAD5-region query returned TGFBI–vascular phenotype pairs, not SMAD5–sleep/frailty results. See `stage5_cascade_regional_coloc_review.md` for scope, raw responses and limitations.

## Additional brain-cell molecular evidence

The public SingleBrain top-association files independently add brain expression context: four queued candidate alleles are top SMAD5 associations in excitatory neurons, inhibitory neurons, oligodendrocytes and their precursors, and three are top TMEM161B-DT associations. These are distinct variants and are not evidence of shared disease causality or proof of cell-type specificity. Only one top variant per gene and cell type is available in this screen. See `stage5_singlebrain_top_association_review.md` and `stage5_singlebrain_gene_top_associations.csv` for allele-level details and limitations.

The SingleBrain browser was subsequently queried at variant level for four priority genes, yielding 584 matched variant–gene–cell rows for 104 distinct queued alleles. For the top SMAD5 excitatory-neuron eQTL chr5:136159634:A>G, the aligned disease GWAS results are nominal and postselected; no shared causal effect is established. Query scope, reproducibility and mixed-statistic limitations are recorded in `stage5_singlebrain_priority_variant_review.md`.

For the SMAD5 excitatory-neuron data, the underlying public full-statistics file was accessed by byte range and its effect-allele and standard-error fields were audited. All 35 portal/queue matches agree in beta and have ALT as their source effect allele, but the retrieved source has a large positional gap encompassing the original disease marker. It does not support a complete regional colocalization. One other source row has an effect allele inconsistent with its REF/ALT. See `stage5_singlebrain_ext_smad5_source_review.md`.

An underlying SingleBrain component study, Bryois et al., provides separately released excitatory-neuron SMAD5 statistics that include the original rs7734045 marker and 1,621 SMAD5 records in the apparent SingleBrain release gap. The marker has nominal p=0.038 in that selected molecular test; this is not an independent replication of SingleBrain and does not establish shared disease causality. The checksum-verified lookup and exact-allele matches are documented in `stage5_bryois_exc_smad5_review.md`.

A feasibility check of two public PsychAD brain vascular expression-prediction databases found two of 97 previously nominated genes meet the paper’s model-quality thresholds in mural cells and none in endothelial cells. This is model eligibility, not a disease association, clinical validation, or a screen of neuronal models. See `stage6_psychad_small_model_feasibility.md`.

A subsequent exploratory S-PrediXcan analysis used the eligible PsychAD RBM6 brain mural-cell model with seven public GWAS. It yielded stronger associations with UK Biobank frailty measures than with FinnGen clinical insomnia or hospital frailty; uncalibrated p values, post-selection, LD and cohort overlap prevent claims of shared causality or independent validation. The full results and independent formula check are in `stage6_psychad_rbm6_twas_review.md`.

The exploratory RBM6 gene-model associations were sensitive to removal of rs2353586, which weakens the Jansen insomnia association especially. The removal changes model variance and is not a causal SNP test; see `stage6_psychad_rbm6_sensitivity_review.md`.

A direct five-GWAS comparison in a two-megabase RBM6-centered window found 2,256 exactly aligned common SNPs. The Jansen insomnia and FinnGen hospital-frailty peak markers in this restricted set have low pairwise LD in the public European reference, but this is not a colocalization test and does not prove distinct causal signals. See `stage6_rbm6_regional_shared_risk_review.md`.

The hospital-frailty common-set peak rs2282749 is also an RBM6 excitatory-neuron eQTL in SingleBrain; the ALT allele points toward lower expression and higher hospital-frailty score, but marginal association and directions do not establish colocalization or mediation. The influential model SNP rs2353586 has unusable zero-SE SingleBrain excitatory-neuron rows and a partial CASCADE immune response with multiple top genes. See `stage6_rbm6_molecular_followup.md`.

CASCADE reports precomputed RBM6 immune eQTL colocalization signal-pair results for pain-related FinnGen endpoints, including CD4 results with PP.H4 around 0.88. The portal results are selected, phenotypes and cells overlap, and no literal sleep/frailty code appeared among the queried rows; these are not sleep–frailty colocalizations or independent validation. See `stage6_rbm6_cascade_coloc_review.md`.

A direct exploratory SuSiE-RSS feasibility attempt on the RBM6 region produced two insomnia candidate sets but no reported hospital-frailty credible set using a 503-person European reference. The common-SNP subset, rank-deficient out-of-sample LD and Z/LD inconsistency constrain interpretation; no sleep–frailty SuSiE-coloc result is claimed. See `stage6_rbm6_multisignal_feasibility_review.md`.

A new two-cohort FOXP2-region SuSiE-coloc analysis of MVP insomnia and FinnGen hospital frailty is highly sensitive to reference LD: an apparently shared EUR-reference signal was not reproduced by the much smaller Finnish subset, on the same SNP support. Neither reference settles the causal question. See `stage6_foxp2_independent_cohort_coloc_review.md`.

A larger public Finnish SISu 4.2 pairwise LD lookup supports strong linkage of the FOXP2 candidate with the EUR-fit lead and identifies differences from FIN99 for a selected credible set. Pairwise r² is not a signed full LD matrix and does not settle the colocalization discrepancy. See `stage6_foxp2_sisu_ld_reference_review.md`.

A larger public UKB European signed-LD sensitivity yielded a high modeled shared-signal posterior on 553 shared HapMap3+ variants, with similar results across the three references on this smaller subset. Both original lead variants were absent, the UKB matrix needed a PSD correction, and the full 3,019-variant reference discrepancy remains unresolved. See `stage6_foxp2_ukb_ld_sensitivity_review.md`.

A narrower FOXP2-region MVP insomnia–FinnGen HFRS analysis using public gnomAD NFE and FIN signed LD gave PP.H4=0.979 for one signal pair under a selected prior (0.826–0.996 over the originally tested shared priors 1e-6 to 5e-5; an expanded post hoc range 1e-8 to 5e-5 gives 0.045–0.996). The result is postselected, excludes a distal signal, relies on external LD and does not establish a causal gene or validate the risk score. See `stage6_foxp2_gnomad_matched_coloc_review.md`.

FOXP2 expression follow-up using the same larger references favors separate pituitary eQTL versus MVP insomnia and HFRS signal pairs; excitatory-neuron FOXP2 versus HFRS similarly favors distinct signals at the illustrated prior, while the neuron/MVP comparison yielded no reported credible set. Missing variants, including a pituitary credible-set indel, and uncertainty in brain effective N limit these findings. See `stage6_foxp2_expression_coloc_review.md`.

The previously omitted FOXP2 pituitary eQTL credible-set insertion was matched in both disease GWAS and the gnomAD FIN/NFE signed LD matrices. Adding it covers all four original unique credible-set alleles and leaves the modeled preference for distinct expression versus disease signals largely unchanged; see `stage6_foxp2_pituitary_indel_coloc_review.md`.

### Expanded brain-expression versus insomnia comparison

A direct FOXP2 excitatory-neuron expression–MVP insomnia analysis retained 857 exact reference-matched variants, including all 720 from the earlier three-way analysis. Both SuSiE fits converged, but the insomnia fit still returned no reported credible set, so a multi-signal colocalization posterior could not be computed. This does not establish either shared or distinct causal signals. The SingleBrain insertion record has zero standard error and cannot be used as an expression estimate. Details and diagnostics: `stage6_foxp2_brain_mvp_pair_review.md`.

### Exploratory 15-marker weighted association with HFRS

A newly calculated reference-LD-adjusted weighted-Z test of the three existing 15-marker weight sets against FinnGen HFRS was positive for all three and remained positive after omitting rs12536335. This is a post hoc summary-statistics association check, with approximate null calibration and possible selection/sample overlap; it is not independent validation of the sums or a comparison of their predictive performance. See `stage6_hfrs_sparse_score_summary_review.md` for the assumptions and complete sensitivity results.

For the exploratory 15-marker HFRS weighted-Z analysis, replacing FIN99 reference LD with the overlapping EUR503 panel and bootstrapping the FIN99 reference gave similar statistics. This sensitivity check does not remedy post-selection, phenotype differences, study overlap or lack of target individual data; see `stage6_hfrs_sparse_score_summary_review.md`.

### Device-measured sleep context

Public Portas et al. (2026) sex-combined UK Biobank REM/NREM statistics have exact allele matches for the 15 project candidates. The MEIS1-region rs113851554 T allele is associated with shorter inferred REM and longer inferred NREM, while the FOXP2-region rs12536335 G allele has a postselected NREM association. A 1,225-SNP FOXP2 SuSiE analysis returned no NREM credible set, so no sleep-stage–disease colocalization was computed. Possible UK Biobank overlap, sleep-stage inference error, selection and allele representation prevent treating these findings as independent validation. The paper's female-only NREM rs1799252 result is a separate analysis; its pairwise LD with the project leads could not be established from the queried reference. See `stage6_device_sleep_followup_review.md`.

### External sleep-phenotype evidence at MEIS1

Published polysomnography and sleep-clinic studies link the project's rs113851554 to periodic limb movements and raise the possibility of a restless-legs contribution to insomnia associations. Their cohorts and phenotypes differ, and two movement publications reuse Wisconsin and MrOS participants. This supports follow-up of phenotype specificity, not shared sleep–frailty causality or validation of the project score. Details, cohort overlap and a contrasting clinical comparison are in `stage6_meis1_independent_sleep_phenotype_review.md`.

A separate published HypnoLaus/MESA replication analysis reported the MEIS1 rs113851554 movement association, although source text and table p-values differ; this is replication for periodic movements, not the project score or frailty. A fresh GWAS Catalog v2 availability audit found no accessible full movement summary statistics in its two entries and no study records for four targeted frailty publications. These are scoped Catalog results, not proof of global unavailability. See `stage6_meis1_independent_sleep_phenotype_review.md` and `stage6_gwas_catalog_v2_access_review.md`.

### Musculoskeletal factor context

A complete stream of Gong et al.'s public degenerative-musculoskeletal factor GWAS matched 12 of the 15 selected candidates. The FOXP2-region rs12536335 and chromosome 3 rs10865954 had strong factor associations. The factor uses UK Biobank/FinnGen or potentially overlapping sources and is not FI, FFS or independent score validation. The paper also reports RBM6 expression colocalizations in selected immune-cell labels, but its 89-donor CD4 single-cell description and broader cell labels leave source provenance unresolved. See `stage6_gong2026_musculoskeletal_followup.md`.

### Musculoskeletal/insomnia/frailty locus comparison

An exploratory 790-SNP FOXP2 analysis gave high modeled pairwise shared-signal posteriors for the musculoskeletal factor versus MVP insomnia and versus FinnGen HFRS. However, it retained only 6 of 18 SNPs (47.8% of original component mass within the credible set) from the original MVP insomnia credible set, and the factor uses overlapping/potentially overlapping studies and a Genomic SEM effective N. These outputs are highly provisional, not joint three-trait colocalization or independent validation. See `stage6_gong_foxp2_signal_sharing_review.md`.

### FOXP2 FI cohort discrepancy

An exploratory 1,187-variant comparison of the original FI cohort-level statistics found that the UKB/TwinGene opposite-direction pattern at rs12536335 extends across 12 correlated high-LD markers, conditional on the project's unconfirmed TwinGene sign adapter. The Swedish SuSiE fit returned no credible set, so no FI cohort colocalization was calculated. This does not provide independent score validation or establish a biological direction difference. See `stage6_foxp2_fi_cohort_regional_review.md`.

### HRS–ELSA FI locus and external sleep lookup

The published HRS–ELSA FI study reports a nominally replicated GRIN2B-region rs7134291 signal. An exact-allele lookup found imprecise associations in MVP insomnia (OR 1.0052, p=0.465) and FinnGen HFRS (beta −0.0196, p=0.145); neither supports a null conclusion or score validation. Its 32 selected published FI rows contain no exact rsID from either project's prototype, which does not imply the candidates were absent or null in the cohort. See `stage6_hrs_elsa_grin2b_cross_trait_review.md`.

### RBM6 larger-reference follow-up

In a narrowed 1-Mb, 1,034-SNP RBM6-region analysis, one insomnia–FI signal pair had high modeled PP.H4 at the illustrated prior, but both source GWAS use overlapping UK Biobank data and the FI includes insomnia. The narrower window retained only 4.6% of one earlier insomnia credible set's component mass and omitted the earlier insomnia peak. HFRS yielded no credible set, so no insomnia–HFRS posterior was calculated. These are not independent confirmations or gene mediation evidence. See `stage6_rbm6_gnomad_coloc_review.md`.

### RBM6 excitatory-neuron eQTL follow-up

An exploratory 541-SNP brain eQTL–disease fit reported candidate signal pairs with insomnia and FI, but retained only about 23–25% of the probability mass within the prior insomnia credible sets and 0–71% across FI sets. Brain credible-set counts changed with assumed sample size; FFS and HFRS returned no sets. The apparent pairwise overlap does not establish mediation or independent validation. See `stage6_rbm6_brain_molecular_coloc_review.md`.


A GTEx bulk-cortex RBM6 follow-up retained the two selected earlier insomnia and FI credible sets, but omitted the published cortex eQTL set’s highest-PIP variant, rs10049477 (PIP 0.503). Its pairwise colocalization posteriors are fragile and do not validate the risk scores or establish a cell-specific mechanism. See `stage6_rbm6_gtex_cortex_coloc_review.md`.


An expanded GTEx cortex–frailty analysis includes the previously omitted rs10049477 and 44 of 49 unique published RBM6 eQTL-set variants. It favors sharing with one FI component under the selected model and distinct signals with FFS, with FFS allele frequency unavailable for reference QC. This does not establish independent score validation or resolve insomnia colocalization on the expanded variant set. See `stage6_rbm6_gtex_direct_coloc_review.md`.


A direct insomnia–cortex comparison expanded variant support to 985 but still omitted the leading published RBM6 expression variant; its similar conditional posterior cannot establish independent confirmation. See `stage6_rbm6_gtex_jansen_coloc_review.md`.


An expanded 1,000-variant cortex–MVP insomnia analysis includes the published leading RBM6 expression variant and favors distinct signals at the illustrated prior; the result is prior-sensitive and does not test the scores. See `stage6_rbm6_gtex_mvp_coloc_review.md`.


On a shared 960-variant set, Jansen and MVP insomnia show different pairwise RBM6 cortex colocalization posteriors, but the set misses the leading published expression variant and is not a formal heterogeneity or individual-score validation test. See `stage6_rbm6_sleep_matched_review.md`.


A post hoc reference-LD-adjusted 15-marker weighted-Z check against MVP and FinnGen R13 insomnia is positive under all three existing weight sets. This is summary-data association under strong assumptions, not a prespecified replication or individual-score validation; see `stage6_insomnia_sparse_score_summary_review.md`.
