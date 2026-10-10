#!/usr/bin/env python3
"""Assemble bounded reviewer judgments from independently inspected small artifacts.

These are machine-readable evidence/figure descriptions, not manuscript legends.
No estimates, data transformations, plotted figures, or source bodies are produced.
"""
import csv
import datetime
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
V4=ROOT/'sleep_unified_research_v4'
H='discovery_extension/sleep_submission_evidence_v1/'
V1='sleep_unified_research_v1/'
V2='sleep_unified_research_v2/'
V3='sleep_unified_research_v3/'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(65536),b''):h.update(b)
    return h.hexdigest()

def write_json(path,value):
    path.write_text(json.dumps(value,indent=2)+'\n')

def write_tsv(path,records):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader()
        for row in records:
            w.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in row.items()})

def main():
    now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    probe=json.loads((V4/'qc/independent_figure_export_probe_v4.json').read_text())
    inventory=json.loads((V4/'statistical_validation/independent_table_inventory_v4.json').read_text())
    lookup={x['path']:x for x in inventory['tables']}
    descriptions=[
      {'id':'01_unified_study_design','scope':'Three historical testing families and execution status',
       'factual_values':{'core_tests':396,'core_primary_fdr_positive':153,'core_all_fdr_positive':161,'extension_tests':1200,'extension_fdr_positive':603,'validation_locked':217,'validation_estimated':41,'validation_outcome_side_positive':23,'validation_concordant':18,'validation_qc_ineligible':17,'validation_absent':159,'processed_input_pilot_pairs':1,'both_trait_independent_validations':0},
       'encoding':['Three text boxes; counts and execution status','Color emphasizes incomplete execution'],
       'uncertainty':['Counts have no sampling interval','One pilot does not mean one full raw-to-estimator pair'],
       'issues':['Preserve the v1 snapshot date; current raw-source/preprocessing progress is not shown','Any v4 replacement should explicitly say processed-input pilot; native replay alone can be read too broadly'],
       'use':'research status graphic; update snapshot wording before final scientific use'},
      {'id':'02_global_core_atlas','scope':'12 sleep traits × 33 outcomes, original 396-test family',
       'factual_values':{'cells':396,'bh_positive':161,'primary_positive':153,'sensitivity_positive':8,'sensitivity_rows':['T2D','melanoma'],'color_scale':[-1,1],'data_rg_range':[-.2976,.6405]},
       'encoding':['Cell color is historical global rg','Dot is within-396-family BH q < 0.05','Hatching is historical QC-failed sensitivity status'],
       'uncertainty':['SE exists in source table; no per-cell interval is drawn','Threshold dots do not show precision or covariance','No equal-power comparison of objective and self-report traits'],
       'issues':['Source values are historical serialized estimates; one pilot has separately preserved full precision','Mixed vector/raster PDF/SVG; flat cell raster is 12×33, not a photographic resolution requirement','Keep liability-scale and lipid QC qualifications linked to the table'],
       'use':'qualified historical atlas; final scientific interpretation pending native/QC review'},
      {'id':'03_global_extension_atlas','scope':'12 sleep traits × 100 extension outcomes, original 1,200-test family',
       'factual_values':{'cells':1200,'bh_positive':603,'color_scale':[-1,1],'data_rg_range':[-.4098,.5881],'duplicate_label_groups':4,'rows_with_duplicate_display_labels':10},
       'encoding':['Cell color is historical global rg','Dot is within-1,200-family BH q < 0.05'],
       'uncertainty':['SE exists in source table; no per-cell interval is drawn','Selection, cohort overlap and power differences are not solved by the heatmap'],
       'issues':['Single PDF is approximately 12.12×22.98 inches; 7-point row labels become roughly 3.35 points if height is scaled to 11 inches','Four repeated label groups: Operation code (3), Pain type(s) experienced in last month (2), Treatment/medication code (2), Mouth/teeth dental problems (3)','Category code/source ID is needed to distinguish those rows','Mixed vector/raster PDF/SVG; flat cell raster 12×100 and raster colorbar'],
       'use':'zoomable historical supplementary atlas; pagination or defined panel split and explicit codes required for print'},
      {'id':'04_validation_evidence','scope':'217 locked historical validation candidates',
       'factual_values':{'outcome_side_positive':23,'concordant':18,'qc_ineligible':17,'no_eligible_external_source':159,'estimated_pairs':41,'both_trait_independent':0,'conditional_MVP_candidates':10},
       'encoding':['Bar lengths show mutually exclusive historical classes','Large zero shows both-trait independence count'],
       'uncertainty':['Class counts do not measure replication precision','All 41 estimated rows reuse discovery sleep inputs','QC category combines 9 low-h2-Z and 8 high-intercept exclusions; underpowered is an overbroad historical label'],
       'issues':['MVP has since completed result-free preprocessing; its h2/rg and all ten conditional candidate gates remain pending','Do not describe the 23 as fully independent replication'],
       'use':'qualified outcome-side validation design graphic'},
      {'id':'05_question_admission_audit','scope':'Five prospective question classes; administrative ordinal scores',
       'factual_values':{'questions':5,'score_range':[0,4],'admitted_primary_question':0},
       'encoding':['Numbers and color encode reviewer ordinal scores; no biological effect sizes'],
       'uncertainty':['Scores have no validated quantitative scale or sampling distribution'],
       'issues':['Question labels intentionally truncate with ellipses; complete IDs/questions must remain available in the linked feasibility table','Score provenance is an automated review, not human scientific approval','Mixed vector/raster heatmap'],
       'use':'technical review packet or supplement; no scientific-result figure'},
      {'id':'06_covariance_heterogeneity_diagnostic','scope':'41 historical estimates, including 23 outcome-side positives',
       'factual_values':{'scatter_points':41,'positive_subset':23,'nominal_heterogeneity_all41':16,'nominal_heterogeneity_positive23':7,'histogram_z_range':[-5,5],'omitted_histogram_points':0},
       'encoding':['Scatter x/y are discovery/outcome-side rg','Horizontal/vertical errors are ±1.96 marginal SE','Dashed scatter line is equality','Histogram is effect-difference Z assuming Cov = 0; red lines approximately ±1.96'],
       'uncertainty':['Marginal CIs are not a covariance-aware interval for the difference','Shared sleep inputs, unknown participant intersections and discovery selection make the heterogeneity test uncalibrated','Nominal 16/41 and 7/23 flags are not confirmed biological differences'],
       'issues':['Keep Cov = 0 and uncalibrated wording in any replacement','No sampled block-coordinate alignment proof or joint covariance is supplied by this graphic'],
       'use':'supplementary diagnostic only; no confirmatory heterogeneity claim'},
      {'id':'07_core_QC_intercept_diagnostic','scope':'33 outcomes; standalone intercept and range across 12 sleep pairings',
       'factual_values':{'outcomes':33,'lipid_pairs_above_1p2':36,'lipid_ranges':{'LDL':[1.222,1.247],'HDL':[1.408,1.423],'TG':[1.205,1.214]},'threshold':1.2},
       'encoding':['Open circle is standalone outcome LDSC intercept','Horizontal segment is pairwise min/max range','Red highlights three lipid outcomes; dashed line is QC threshold'],
       'uncertainty':['Pairwise min/max is not a confidence interval','Standalone and pairwise SNP sets/two-step estimator conventions differ','Native rsID prefilter reproduction does not resolve intercept inflation'],
       'issues':['Historical tiers are preserved, not newly cleared','Intercept SE is available in source results but not drawn'],
       'use':'supplementary unresolved-QC diagnostic'},
      {'id':'MVP_source_preprocessing_QC_v3','scope':'Result-free MVP GIA clinical-insomnia preprocessing',
       'factual_values':{'source_rows':19703815,'retained_rows':826027,'exclusion_categories':6,'P_CI_discrepancy_buckets':7,'h2_estimated':False,'rg_estimated':False},
       'encoding':['Panel A horizontal bars use a log row-count axis and exclusive exclusions','Panel B is the percentage of retained rows by relative discrepancy from source-reported P'],
       'uncertainty':['CI-derived Z assumes a 95% log-scale Wald CI; source SE is absent','Discrepancy buckets are descriptive and do not filter retained rows','826,027 retained variants do not demonstrate estimator power, phenotype equivalence or independence'],
       'issues':['No h2/intercept/power gate has passed','Worker/reference/source seals are carried from independently checked v3; this pass rehashes only small receipt/table/figure artifacts'],
       'use':'supplementary preprocessing QC graphic; no biological result'}]
    for d,p in zip(descriptions,probe['figures']):
        assert d['id']==p['id']
        d['input_tables']=[{'path':s,'sha256':lookup[s]['sha256'] if s in lookup else None} for s in p['sources']]
        d['exports']=p['exports']
        d['visual_QA']={'PNG_native_inspected':True,'PDF_render_inspected':True,'clipping_or_overlap_observed':False,'legible_at_native_export':True,'extension_print_legibility_exception':d['id']=='03_global_extension_atlas','review_type':'automated visual inspection; human author review pending'}
        d['accessibility']={'svg_has_text':True,'svg_title_or_desc_present':False,'pdf_tagged':False,'formal_color_vision_test_performed':False,'numerical_table_alternative_available':True}
        plot_script=V3+'scripts/24_make_mvp_qc_figure.py' if d['id']=='MVP_source_preprocessing_QC_v3' else V1+'scripts/08_make_evidence_figures.py'
        figure_manifest=V3+'manifests/mvp_QC_figure_receipt_v3.json' if d['id']=='MVP_source_preprocessing_QC_v3' else V1+'manifests/figure_manifest_v1.json'
        d['generation_provenance']={'script_path':plot_script,'script_sha256':sha(ROOT/plot_script),'manifest_path':figure_manifest,'manifest_sha256':sha(ROOT/figure_manifest),'matplotlib':'3.10.7','numpy_recorded':'2.5.1' if d['id']!='MVP_source_preprocessing_QC_v3' else 'not separately recorded in v3 plot manifest','font_family':'DejaVu Sans','pdf_fonttype':42,'svg_fonttype':'none','complete_environment_or_font_binary_lock':False}
        d['publication_readiness']='PARTIAL'
        d['description_kind']='factual metadata; not manuscript caption or alt-text prose'
    write_json(V4/'qc/independent_figure_readiness_v4.json',{'baseline_commit':'2ef91a5e','reviewed_utc':now,'figure_count':8,'original_exports_preserved':True,'figures_regenerated':False,'plots_are_scientific_validation':False,'official_guidelines_url':'https://academic.oup.com/sleep/pages/author-guidelines','guideline_access_note':'Official-domain search retrieved current author-guidelines text on 2026-10-09; direct page opens timed out. Parent is independently checking live requirements.','format_resolution_findings':['All eight have PDF, SVG and PNG exports','PDF is a listed SLEEP main-figure format; PNG/SVG are useful research exports','All v1 PNGs are approximately 220 dpi; v3 PNG is 240 dpi, below the stated 300dpi raster minimum','Atlas and question-audit PDFs/SVGs contain rasters; the remaining five PDFs have no embedded images listed by Poppler','No source-body hashes were recomputed for this inspection'], 'figures':descriptions})
    write_tsv(V4/'qc/independent_figure_readiness_v4.tsv',[{'figure_id':d['id'],'scope':d['scope'],'readiness':d['publication_readiness'],'use':d['use'],'issues':d['issues'],'uncertainty_limits':d['uncertainty'],'factual_values':d['factual_values'],'source_tables':d['input_tables']} for d in descriptions])

    evidence=[]
    def item(domain,status,scope,coverage,paths,missing,dependencies,limit):
        evidence.append({'evidence_domain':domain,'status':status,'expected_scope':scope,'observed_coverage':coverage,'artifacts':[{'path':p,'sha256':lookup[p]['sha256'] if p in lookup else sha(ROOT/p) if (ROOT/p).is_file() else None,'rows':lookup[p]['rows'] if p in lookup else None,'columns':lookup[p]['columns'] if p in lookup else None} for p in paths], 'missing_columns_or_values':missing,'dependencies_to_complete':dependencies,'interpretation_limit':limit})
    item('GWAS source identities','PARTIAL','45 core trait roles + 100 extension + 13 historical validation outcomes + MVP candidate (159 trait/source-role records; not 159 unique GWAS)',
         'Core registry: 43 containers; v3 42 original containers matched and 45 raw projections available. Extension 100 and historical validation 13 processed inputs recovered. Rich metadata125 = 12sleep+100extension+13validation; thin unified metadata145 =45core+100extension.',
         ['config/public_gwas_sources.tsv',V1+'tables/phenotype_and_source_metadata.tsv',H+'tables/phenotype_and_source_metadata.tsv',V3+'tables/CURRENT_EVIDENCE_STATUS_V3.json'],
         ['Single current joined table lacks container SHA, raw projection SHA, processed SHA, source/release/accession, code/receipt and temporal status together for all159 roles','Historical rich metadata current_raw_source_integrity has112NA; do not read it as current v3 availability','Original prostate container hash is unresolved despite retained raw projection'],
         ['Use exact v3 current source ledgers and provenance reviewer updates; retain original fields separately','Resolve original prostate bytes or preserve a specifically bounded qualified projection chain'], 'HEAD/public access and processed input recovery are not full raw-source authentication or rights approval.')
    item('Phenotype definitions','PARTIAL','All traits and candidate source with construct, direction, ascertainment and measurement',
         'Core45 definitions; richer125 metadata definitions/effect conventions; sleep12 measurement table; v2/v3 MVP source contract.',
         ['config/analysis_panel.tsv','config/sleep_measurement_modes.tsv',H+'tables/phenotype_and_source_metadata.tsv',H+'tables/sleep_trait_metadata.tsv',V2+'source_definition_review/provenance_v2_definition_contract.json'],
         ['Unified145 table lacks cases, controls, effective_n and effect_convention','Richer125 excludes33 non-sleep core outcomes and new MVP candidate','Uniform fields for collection window, respondent, diagnostic coding version, source covariates and cross-source estimand equivalence are not complete','Frozen core mdd build is UNRESOLVED; current qualified source determination must be joined rather than silently changed'],
         ['Current source/phenotype review, source primary documentation and sleep expert adjudication'], 'Name/ICD resemblance is not identical ascertainment or estimand; actigraphy is not PSG and morningness is not measured biological phase.')
    item('Cohort overlap','PARTIAL','Four discovery/validation intersections plus estimator covariance',
         'v1 table records all217 four-intersection statuses; all41 estimated use same discovery sleep input;58 external-source mappings are cohort-source distinct, individual intersections unknown.',
         [V1+'REPLICATION_RESULTS.tsv',V1+'REPLICATION_SOURCE_LEDGER.tsv',H+'tables/replication_source_comparison.tsv'],
         ['Numeric intersection counts, participant-level or study-specific evidence and cross-estimator covariance remain unknown','Legacy comparison table participant_overlap_status=NON_OVERLAPPING_CONFIRMED is superseded by current cohort-source-distinct qualification'],
         ['Cohort composition records; genuinely independent sources or calibrated covariance-aware design'], 'Distinct outcome cohort does not make both-trait replication independent.')
    item('Original trait heritability','COMPLETE_WITH_QUALIFICATIONS','45 core +100 extension +13 validation historical h2 rows',
         '45 and100 dedicated tables;13 validation h2 source rows; frozen core43pass/2drop and extension100pass.',
         [V1+'sources/recovered/results/tables/h2_summary.tsv',V1+'tables/original_extension_h2_100.tsv','discovery_extension/results/replication/replication_source_h2.tsv'],
         ['Historical full-precision h2 absent for most rows','Single combined historical/current native table with N convention, liability assumptions, allele/reference code and execution receipt is not complete'],
         ['Full native h2 receipts for158 planned h2 jobs and scale/QC adjudication'], 'Historical numerical/log concordance is not native rerun; liability h2>1 is not a physical proportion interpretation.')
    item('All396 core correlations','COMPLETE_WITH_QUALIFICATIONS','396 original pairs, original BH396 and preserved tiers',
         '396unique, allrg/SE/Pfinite, BH161 =153primary+8sensitivity;372primary rows24sensitivity.',
         [V1+'tables/original_core_396.tsv',V1+'tables/original_vs_native_pilot.tsv'],
         ['Validated native full-precision columns absent for395pairs','Per-pair execution receipt, before/after raw/QC/processed/reference/code hashes absent for whole campaign'],
         ['57 core jobs and raw-to-estimator chain; independent precision-aware comparison'], 'This is complete historical table coverage; not completed full native396 reproduction.')
    item('All1200 extension correlations','COMPLETE_WITH_QUALIFICATIONS','1,200 original pairs and original BH1200',
         '1200unique allrg/SE/Pfinite;603BHpositive; richer derived table contains CIs/printed precision bounds.',
         [V1+'tables/original_extension_1200.tsv',H+'tables/global_1200.tsv'],
         ['Validated native full-precision estimates/SE/P and per-pair execution receipt are absent for1200pairs','Raw-source body authentication is not implied by100HEAD matches'],
         ['112 extension native jobs plus original/raw-to-processed evidence'], 'Historical CI is based on serialized SE; no new native or primary interpretation.')
    item('All41 estimated validation correlations','COMPLETE_WITH_QUALIFICATIONS','All41 estimated rows, including18 below validation threshold',
         '41-rowHETEROGENEITY_MASTER includes all41 estimates and historical CIs/difference/Cov0 sensitivity grids.',
         [H+'tables/HETEROGENEITY_MASTER.tsv',V1+'REPLICATION_RESULTS.tsv'],
         ['Native full-precision values/receipts absent41','Validated joint covariance and covariance-aware differenceCI/SE/P absent'],
         ['21 historical validation jobs and joint/overlap design review'], '23positive+18concordant are outcome-side validation with sleep reuse.')
    item('Original217 candidate classifications','COMPLETE_WITH_QUALIFICATIONS','Original locked217 retained regardless availability/QC/result',
         '217unique:23positive18concordant17QC159absent;17QC =9h2Zfail+8interceptfail.',
         [V1+'REPLICATION_RESULTS.tsv',V1+'REPLICATION_SOURCE_LEDGER.tsv',H+'tables/replication_family_217.tsv'],
         ['Current joined classifications must retain source snapshot date and current QC reasons','Historical UNDERPOWERED label is overbroad for8high-intercept exclusions'],
         ['Do not shrink217 denominator; preserve original and current qualified interpretation separately'], 'Classification coverage is complete, independence and scientific eligibility are not.')
    item('Native reproduction comparisons','PARTIAL','190 planned job receipts:158h2+32rg batches;396/1200/41pair estimates',
         '0 campaign execution receipts in sealedv1-v3;1processed-input pilot10metriccomparisonrows;6successful preprocessing units (3lipids,2cancer,1MVP).',
         [V1+'manifests/native_reproduction_jobs_v1.json',V1+'tables/original_vs_native_pilot.tsv',V3+'tables/NATIVE_CANCER_RAW_CONTENT_COMPARISON_V3.tsv'],
         ['For every planned job: execution status, command, input/reference/code before+after SHA, original serialized/native full precision estimate/SE/P, precision convention, QC verdict and discrepancy class','Six preprocessing units must not populate estimator-completion columns'],
         ['3GiB internal launch floor; original data and dependencies; scoped native launch review; independent per-pair numerical audit'], 'No full raw-to-estimator core/extension/validation pair completed. Pilot is separately qualified.')
    item('Lipid statistical QC','PARTIAL','36core lipid pairings and proposed lipid-conditioned model',
         'All36pairwise outcome intercepts>1.2; standalone LDL1.141HDL1.196TG1.165; exact3native rsID prefilters verified.',
         [V1+'tables/core_intercept_stage_diagnostic.tsv',V1+'tables/original_core_396.tsv',V2+'tables/NATIVE_LIPID_PREFILTER_COMPARISON.tsv',V1+'reviews/statistical_validity_v1.md'],
         ['Resolved estimator/SNP-set explanation, native matched-rule sensitivity estimates/SE/CI and final scientific inclusion verdict absent','No calibrated lipid-conditioned primary outcome table'],
         ['Full native/source-harmonization evidence; same-SNP-set/two-step diagnostics; methods and statistical review'], 'Matching rsID prefilter bytes does not clear intercept inflation or admit conditioning.')
    item('Liability-scale QC','PARTIAL','29binary core traits and all relevant prevalence/sample conventions',
         'Coreconfig includes29case/control/population prevalence rows and citations; frozen MS liabilityh2=5.982 and melanoma3.731.',
         ['config/analysis_panel.tsv',V1+'sources/recovered/results/tables/h2_summary.tsv',V1+'reviews/statistical_validity_v1.md'],
         ['Unified table of observed versus liability h2, correct effectiveN, ascertainment, prevalence sensitivity grid and adjudicated assumptions not completed','Current primary-eligible scale interpretation remains unresolved'],
         ['Verified case/control/effect/N convention; native h2 and feasible prevalence sensitivity'], 'Positive liability rescaling leaves rg invariant but does not authenticate h2 magnitude or biological variance claims.')
    item('Covariance validity','PARTIAL','Covariance-aware joint estimator contrasts with aligned genome blocks',
         'One pilot200delete arrays independently checked; historical exportedGenomicSEM matrices and Cov0/rho grids audited.',
         [V1+'reviews/independent_numerical_v1.json',V1+'reviews/statistical_validity_v1.json',H+'tables/HETEROGENEITY_MASTER.tsv'],
         ['Genome coordinate/order/SNP-set/block alignment proof and calibrated joint sampling covariance for targeted contrasts absent','Validated rg uncertainty propagation absent for standardized covariance export'],
         ['Joint native estimators with exact aligned block contract or another validated covariance method; matrix numerical and statistical checks'], 'Equal-length arrays or matching delete count do not prove genomic alignment; Cov0 grid is a diagnostic.')
    item('New primary analysis','NOT_ADMISSIBLE','One genuinely eligible prospectively supported sleep-focused primary analysis',
         'No question admitted; no new h2/rg or pair outcomes; result-free MVP preprocessing only.',
         [V1+'HYPOTHESIS_FEASIBILITY_MATRIX.tsv',V1+'NEW_SCIENTIFIC_RESULTS.tsv',V3+'tables/CURRENT_EVIDENCE_STATUS_V3.json'],
         ['Admittedquestion, estimand, prespecifiedfamily, valid covariance, power/QC gates, estimate/SE/CI/P, native receipt, sensitivity and corroboration fields are unpopulated'],
         ['Choose only scientifically justified question after source/statistical gates; human adjudication'], 'NO_GO/placeholder rows are not null biological findings or completed analyses.')
    item('Independent replication','BLOCKED_EXTERNAL_DATA','Genuine independence on both trait inputs with phenotype comparability',
         '0completed both-trait independent; historical41reuse sleep; newMVP candidate preprocessed but all estimator/power/eligibility gates pending.',
         [V1+'REPLICATION_RESULTS.tsv',V3+'logs/mvp_preprocessing_execution_receipt_v3.json'],
         ['Eligible independent outcome/source pairing, both-trait intersections, match adjudication, source h2/intercept/power, corrected result and native comparison absent'],
         ['Appropriate independent data plus3GiB execution floor and prespecified gate'], 'Source recovery is not successful scientific replication.')
    item('Measurement-specific contrasts','PARTIAL','Subjective, objective, morningness and disease-like sleep comparisons with uncertainty',
         '12sleep modes and source/effect directions documented; automated phenotyping reviews available.',
         ['config/sleep_measurement_modes.tsv',V1+'reviews/sleep_phenotyping_v1_core.tsv',H+'tables/sleep_trait_metadata.tsv'],
         ['Per-target contrast estimand, rg1/rg2, aligned covariance, differenceSE/CI/P, effectiveN/h2-aware power evaluation, clinical interpretation adjudication absent'],
         ['Phenotype expert review and validated shared-estimator covariance'], 'Counting significant rows cannot compare measurement specificity or power.')
    item('Sensitivity results','PARTIAL','Each major claim: original/sensitivity estimate, uncertainty, QC, result interpretation',
         '281historical audit-summaryrows and41Cov0rho-grid rows; precision/correction and subset audits.',
         [H+'tables/sensitivity_results.tsv',H+'tables/HETEROGENEITY_MASTER.tsv',V1+'reviews/statistical_validity_v1.md'],
         ['Uniform pair-level sensitivity table lacks original_estimate, sensitivity_estimate, difference, sensitivitySE/CI, QC_result and interpretation_change','Native prevalence/intercept/SNP-set and calibrated covariance analyses not complete'],
         ['Feasible native sensitivity plans only after scientific/resource gates'], 'Arithmetic correction/rho grids are not native robustness experiments.')
    item('Null findings','PARTIAL','All qualified negative or imprecise findings retained',
         'Original396/1200tables contain every row;18historical validation estimates below threshold retained.',
         [V1+'tables/original_core_396.tsv',V1+'tables/original_extension_1200.tsv',H+'tables/HETEROGENEITY_MASTER.tsv'],
         ['Explicit non-support reason, precision/power/equivalence assessment and admissibility fields not unified','No new primary scientific null exists because no primary test admitted'],
         ['Valid uncertainty and power assessment; do not equate non-significance with equivalence'], 'Historical non-significance is not no shared genetics.')
    item('Failed analyses','PARTIAL','All source failures, QC exclusions, no-go and blocked executions distinguishable',
         '17validationQC+159absent,24core sensitivity rows; failedHDLsequential acquisition retained and successfulreplacementseparate; rejectedprostate container receipt retained.',
         [V1+'REPLICATION_RESULTS.tsv',V2+'logs/hdl_source_acquisition_receipt_v2.json',V3+'logs/prostate_cancer_acquisition_receipt_v3.json',V1+'HYPOTHESIS_FEASIBILITY_MATRIX.tsv'],
         ['One unified ledger needs stage, attempted/not_attempted, timestamp, exit/error, resource guard, retained outputs, scientific exclusion reason and retry/replacement link'],
         ['Join existing failure receipts without overwriting or counting blocked work as attempted estimator runs'], 'Acquisition failure, software stop and scientific QC failure have distinct meanings.')
    item('Prior-literature comparisons','PARTIAL','Pair-specific comparator and claim novelty qualification',
         '1200row novelty master,115selectedprior-estimatecrosswalk rows;1200adequately_supported_novel_result=NO;521literaturecomparisonunresolved.',
         [H+'tables/NOVELTY_MASTER_1200.tsv',H+'tables/novelty_crosswalk.tsv',H+'tables/NOVELTY_REPLICATED_23.tsv'],
         ['No admittedquestion-specific current search, primary comparator harmonization, quantitatively defensible distinctcontribution and human novelty verdict completed'],
         ['Focused original-source literature review tied to an admissible primary question'], 'Crosswalk coverage or source hashes do not establish novelty.')
    item('Claim-to-evidence mapping','PARTIAL','Every proposed scientific claim bound to input/result/QC/uncertainty/review',
         '1200historical extension claim rows plus qualified v1/v3 status reports; frozen core and newcandidate status in separate artifacts.',
         [H+'tables/full_claim_evidence_ledger.tsv',V3+'tables/CURRENT_EVIDENCE_STATUS_V3.json'],
         ['Single versioned ledger for396core+1200extension+41historicalvalidation+anynewclaim lacks current native receipt, joint uncertainty, source status, rejectedclaim, human approval and finalfigure/table keys together'],
         ['Join audited family tables without claiming scientific completion; preserve historical allowed_interpretation'], 'A technical ledger entry is not a biological discovery.')
    item('Data-sharing permissions','BLOCKED_HUMAN_REVIEW','Dataset-specific analytical-use and redistribution clearance; repository release/ethics approvals',
         '8resource-level license evidence rows and public release boundary;15human approval checklist items all unsupplied; MVP accession-specific terms preserved separately.',
         [H+'sources/provenance_license_evidence.tsv',H+'tables/human_approval_checklist.tsv',H+'16_PUBLIC_RELEASE_AND_LICENSE_AUDIT.md',V2+'source_definition_review/provenance_v2_definition_contract.json'],
         ['Per-source access_license, analytical_use_scope, redistribution_scope, restrictions, acknowledgement, checked_date, responsible_author, approval_evidence and release_boundary are not a complete joined machine-readable ledger for all source roles','No investigator ethics determination, source-owner permission or general redistribution grant documented'],
         ['Human owners/source terms verification; exact accession exceptions and third-party cache exclusion'], 'Public HTTP availability and hash success do not supply new rights or ethical approval.')
    write_json(V4/'statistical_validation/independent_supplementary_completeness_v4.json',{'baseline_commit':'2ef91a5e','reviewed_utc':now,'domains':evidence,'current_review_scope':'Existing sealed v1-v3 evidence; later v4 additions need their own receipt. All table field inventories include numeric coverage and explicit missing tokens. No original tables changed.'})
    write_tsv(V4/'statistical_validation/independent_supplementary_completeness_v4.tsv',[{'evidence_domain':e['evidence_domain'],'status':e['status'],'expected_scope':e['expected_scope'],'observed_coverage':e['observed_coverage'],'artifact_paths':[a['path'] for a in e['artifacts']],'missing_columns_or_values':e['missing_columns_or_values'],'dependencies_to_complete':e['dependencies_to_complete'],'interpretation_limit':e['interpretation_limit']} for e in evidence])
    gates=[
      (1,'Full native analysis reproduction','PARTIAL','190planned/0campaign receipts; processed-input pilot1; six newer preprocessing units','3GiB internal floor and complete per-job/source-to-result execution'),
      (2,'Complete relevant source provenance','PARTIAL','45raw projections,42/43originalcorecontainers;100extension+13validationprocessed; source metadata/current statuses in multiple ledgers','Original prostate container, exact relevant full raw chains and unified current provenance'),
      (3,'Statistical QC','PARTIAL','All36lipidpairwiseintercepts>1.2; liabilityh2>1; unknown joint covariance','Matched-method diagnostics, effectiveN/prevalence review and validated covariance'),
      (4,'Valid primary inference','NOT_ADMISSIBLE','No admittedprimaryquestion/result','Scientific/source/power/covariance gates and prospective plan'),
      (5,'Sleep-phenotype interpretation','PARTIAL','12measurementcategories and reviews; no validmeasurement-specific contrast','Expert phenotype/clinical review plus formal power-aware inference'),
      (6,'Novelty','NOT_ADMISSIBLE','No adequatelysupportednovelresult in1200ledger; currentcandidate primaryquestionunadmitted','Question-specific original-source comparator and human novelty review'),
      (7,'External validation','BLOCKED_EXTERNAL_DATA','41historicalestimatesreuse sleep; bothtraitindependent0; MVPpreprocessingonly','Eligiblebothtrait sources/overlap/phenotype/QC/power and native execution'),
      (8,'Sensitivity analyses','PARTIAL','Historical arithmetic/rho/subset audits preserved; no fullnative sensitivity set','Target-specific estimate/uncertainty/QC comparison'),
      (9,'Figure and table consistency','COMPLETE_WITH_QUALIFICATIONS','Historical numericalcounts/BH and all24exporthashes agree; finalstory/format/layout not complete','Addressprint/duplicate-label/current-status metadata; nofinalprimaryfig6'),
      (10,'Human scientific review','BLOCKED_HUMAN_REVIEW','Automated reviewers are not human approval;15author checklistitems unsupplied','Named statisticalgenetics/sleep/investigator review'),
      (11,'Source permissions','BLOCKED_HUMAN_REVIEW','8resourcelevel terms rows with explicit unresolvedscopes','Source-specific analytical/release rights, ownership,ethics/acknowledgements'),
      (12,'SLEEP-specific reporting preparation','PARTIAL','AdaptedSTROBE/STREGA evidence map; pagepointers authorpending; figure factualmetadata nowbuilt','Research gates/human review; finalformat/alttext/reporting decisions by authors')]
    write_json(V4/'reviews/independent_requirement_audit_v4.json',{'baseline_commit':'2ef91a5e','reviewed_utc':now,'overall_verdict':'NOT_READY_FOR_HUMAN_MANUSCRIPT_WRITING','research_complete':False,'human_scientific_review_completed':False,'gate_statuses':[dict(zip(['gate','requirement','status','evidence','completion_dependency'],g)) for g in gates],'reviewer_scope':'Independent automated reproducibility and figure/table readiness audit; no full native estimator, no human specialist approval, no new scientific finding.'})
    write_tsv(V4/'reviews/independent_requirement_audit_v4.tsv',[dict(zip(['gate','requirement','status','evidence','completion_dependency'],g)) for g in gates])
    print(json.dumps({'figure_metadata':len(descriptions),'supplementary_domains':len(evidence),'requirements':len(gates),'missing_artifact_paths':[a['path'] for e in evidence for a in e['artifacts'] if a['sha256'] is None]},indent=2))

if __name__=='__main__':main()
