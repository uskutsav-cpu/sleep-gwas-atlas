#!/usr/bin/env python3
"""Build a result-free clinical phenotype/source review from bounded metadata.

Reads historical h2 QC and locked candidate metadata; estimates no new GWAS
quantity and deliberately omits all pair association results from its outputs.
"""
from pathlib import Path
import csv
import hashlib
import json
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
V4 = ROOT / "sleep_unified_research_v4"


def read_tsv(rel):
    with (ROOT / rel).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rel, rows):
    path = V4 / rel
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


core_path = "sleep_unified_research_v1/reviews/sleep_phenotyping_v1_core.tsv"
window_path = "sleep_unified_research_v1/reviews/separate_sleep_phenotyping_v1_evidence.tsv"
h2_path = "sleep_unified_research_v1/sources/recovered/results/tables/h2_summary.tsv"
family_path = "sleep_unified_research_v1/REPLICATION_SOURCE_LEDGER.tsv"
outcome_h2_path = "discovery_extension/results/replication/replication_source_h2.tsv"
core = read_tsv(core_path)
windows = {r["trait_id"]: r for r in read_tsv(window_path)}
h2 = {r["trait"]: r for r in read_tsv(h2_path)}
family = read_tsv(family_path)
outcome_h2 = {r["study_accession"]: r for r in read_tsv(outcome_h2_path)}
assert len(core) == 12 and len(family) == 217

wording = {
    "insomnia": "Field 1200: trouble initiating sleep at night or waking in the middle of the night; Usually versus Never/rarely or Sometimes.",
    "sleepdur": "Field 1160: reported hours of sleep in every 24 hours, including naps; whole-hour response.",
    "shortsleep": "Same field 1160; <=6 hours versus 7 or 8 hours; >=9 hours omitted.",
    "longsleep": "Same field 1160; >=9 hours versus 7 or 8 hours; <=6 hours omitted.",
    "chronotype": "Field 1180: self-classification as definitely/more morning versus definitely/more evening; uncertain/prefer-not responses omitted.",
    "sleepiness": "Field 1220: likelihood of unintended daytime dozing or sleep; four ordered categories fitted numerically 1-4.",
    "napping": "Field 1190: daytime nap frequency; Never/rarely, Sometimes, Usually fitted numerically 1-3.",
    "snoring": "Field 1210: participant reports a partner/relative/friend complaint about snoring; Yes versus No.",
    "sleep_apnea": "R9 G6_SLEEPAPNO: hospital discharge/cause-of-death G473 or 3472; control exclusions G6_NARCOCATA and G6_SLEEPDISOTH; no minimum-event count.",
    "sleep_efficiency": "Axivity AX3/GGIR: summed >=5-minute low-z-angle-change sleep episodes divided by elapsed main SPT window; selected ACC_SLEEP_EFF_RAW_SIN.",
    "accel_sleep_duration": "Axivity AX3/GGIR: summed detected sleep episodes within main SPT window, averaged across wear; selected ACC_SLEEP_DUR_RAW_SIN.",
    "sleep_timing": "Axivity AX3/GGIR: midpoint of first-to-last detected sleep episode defining main SPT window, averaged across wear; selected ACC_SLEEP_MIDP_RAW_SIN.",
}
reliability = {
    "insomnia": "Single complaint item; no required chronicity, nights/week, daytime impairment, adequate opportunity or diagnostic exclusion. Controls can report occasional symptoms. No accession-specific test-retest coefficient verified.",
    "sleepdur": "Coarse whole-hour recall of 24h sleep; includes naps. Source excludes sleep medication and <3/>18h responses. Recall and reporting error can differ from actigraphy. No exact-source reliability coefficient verified.",
    "shortsleep": "Dichotomized recall; clinical insufficient sleep syndrome not established. Shares exactly 305742 normal-duration controls with longsleep; no independent-mechanism claim.",
    "longsleep": "Dichotomized recall; hypersomnia diagnosis not established. Shares exactly 305742 normal-duration controls with shortsleep; no independent-mechanism claim.",
    "chronotype": "Preference item affected by interpretation/schedules; not biological circadian phase. Four-week help text if variable. Higher morningness is opposite to higher/later sleep midpoint.",
    "sleepiness": "Single ordinal unintended-dozing item; not fatigue, Epworth score, MSLT, narcolepsy or measured sleep minutes. Official field removed All of the time during initial assessment; release coding/version matters.",
    "napping": "Single ordinal behavior item; no nap length, clock time, intention or opportunity. Sleepiness and napping cannot be equated. Exact source is UKB no-BMI member.",
    "snoring": "Participant is respondent; partner is a reported information source. No measured nightly frequency, respiratory event index, severity or diagnosed apnea. Contact/awareness can affect capture.",
    "sleep_apnea": "Registry capture depends on recognition and health care. No AHI, universal PSG or pure OSA subtype established. ICD9(3472A) has official narcolepsy/cataplexy label anomaly; 1574 code-related cases overlap other rows and are not unique misclassified cases.",
    "sleep_efficiency": "Activity inference can score quiet wake as sleep; no EEG staging. Shared numerator with duration and common SPT boundaries. No exact-source test-retest coefficient verified; algorithm-validation evidence is not error-free individual measurement.",
    "accel_sleep_duration": "Main-period inferred sleep excludes sleep outside window; differs from reported 24h duration including naps. Quiet wake and motion can affect estimates; not elapsed SPT length or PSG sleep stages.",
    "sleep_timing": "Behavioral timing of main sleep, affected by schedules and algorithm boundaries; not melatonin phase, L5 or M10 timing. Transformed effect cannot be read as clock minutes.",
}
clinical_question = {
    "insomnia": "Approve symptom-liability wording and whether clinical transport is meaningful without claiming chronic insomnia disorder equivalence.",
    "sleepdur": "Confirm that 24h duration including naps is the intended comparator, rather than nocturnal-only sleep or sleep opportunity.",
    "shortsleep": "Confirm selected duration-tail interpretation; no diagnosis or U-shaped causal inference from tail results alone.",
    "longsleep": "Confirm selected duration-tail interpretation; no hypersomnia diagnosis or mechanism attribution.",
    "chronotype": "Confirm preference-versus-realized-timing comparison and opposite effect directions; do not infer circadian phase.",
    "sleepiness": "Reconcile response-version change with archived four-category coding; retain propensity-to-doze construct.",
    "napping": "Confirm frequency-behavior scope without equating naps with pathological sleepiness or a physiological subtype.",
    "snoring": "Approve participant-reported close-contact complaint label; do not certify OSA from this item.",
    "sleep_apnea": "Finnish coding expert: resolve 3472/3472A historical code meaning, source label and unique contribution; distinguish special outpatient HILMO from primary outpatient AVOHILMO.",
    "sleep_efficiency": "Measurement expert: compare algorithm/device/window and ascertainment before biological contrast; preserve quiet-wake limitation.",
    "accel_sleep_duration": "Measurement expert: decide whether related 24h-report versus main-period-device transport is useful with time lag and shared UKB samples.",
    "sleep_timing": "Measurement expert: verify clock wrapping, realized-timing interpretation and compatibility of independent sources.",
}
crosswalk = []
for r in core:
    t = r["trait_id"]
    hh = h2[t]
    overlap = ("UKB questionnaire sample overlap and actigraphy-subset overlap; exact intersections not quantified here." if r["cohort"] == "UK Biobank" else "FinnGen R9 participants overlap later FinnGen releases; individual cross-enrollment with UKB/MVP unknown.")
    window = windows[t]["time_window"]
    if t in {"sleep_efficiency", "accel_sleep_duration", "sleep_timing"}:
        window = "Up to 7 continuous wear days, 2.8-9.7 years after baseline; valid days averaged; selected transformed UKB European subset."
    crosswalk.append({
        "trait_id": t, "source_id": r["source_id"], "source_release": r["source_release"],
        "construct": r["reviewed_construct"], "wording_device_and_coding": wording[t],
        "assessment_window": window, "source_effect_units": r["effect_convention"],
        "effect_allele": r["effect_allele"], "effect_column": r["effect_column"],
        "SE_column": r["standard_error_column"], "P_column": r["p_column"],
        "positive_direction": r["higher_effect_interpretation"], "N": r["n_total"],
        "cases": r["ncase"], "controls": r["ncontrol"], "cohort": r["cohort"],
        "ancestry": "EUR" if r["cohort"] == "UK Biobank" else "Finnish European",
        "build": windows[t]["build"], "historical_h2_scale": hh["scale"],
        "historical_h2_rounded": hh["h2"], "historical_h2_SE_rounded": hh["se"],
        "historical_h2_Z": hh["z"], "historical_intercept": hh["intercept"],
        "historical_h2_QC": hh["verdict"],
        "h2_evidence_boundary": "Recovered rounded historical table; not v4 estimate, reliability coefficient or cross-scale phenotype-comparability test.",
        "reliability_ascertainment": reliability[t], "sample_overlap": overlap,
        "remaining_human_question": clinical_question[t],
        "primary_source_urls": windows[t]["primary_source_urls"],
        "primary_definition_evidence": "Existing source caches/locked schema audited; R9 endpoint and code anomaly freshly reverified.",
    })
write_tsv("tables/sleep_clinician_core_crosswalk_v4.tsv", crosswalk)

manifest = json.loads((V4 / "source_provenance/sleep_clinician_finngen_R13_manifest_selection_v4.json").read_text())
fg = next(r for r in manifest["selected_rows"] if r["phenocode"] == "F5_INSOMNIA")
assert fg["num_cases"] == "51643" and fg["num_controls"] == "446273"
locked = [r for r in family if r["sleep_trait"] == "insomnia" and r["replication_study_accession"] in {"GCST90479148", "GCST90479330"}]
assert len(locked) == 2
admissibility = []
for r in locked:
    accession = r["replication_study_accession"]
    oh = outcome_h2[accession]
    historical_receipt = json.loads((ROOT / ("discovery_extension/provenance/replication_streaming_receipts/gwas_catalog_" + accession + ".json")).read_text())
    prospective = "finngen_R13_F5_INSOMNIA__" + accession
    admissibility.append({
        "prospective_source_pair_id": prospective, "locked_historical_pair_id": r["pair_id"],
        "historical_candidate_family": "217; preserve full family, report all selected/QC-failed candidates; no outcome-based family narrowing.",
        "sleep_exact_accession": "FINNGEN_R13_F5_INSOMNIA (release-specific endpoint identity, no GCST accession asserted)",
        "sleep_source_url": fg["path_https"], "sleep_cases": fg["num_cases"], "sleep_controls": fg["num_controls"],
        "sleep_primary_metadata_URLs": "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv;https://r13.risteys.finngen.fi/endpoints/F5_INSOMNIA;https://finngen.gitbook.io/documentation/data-description;https://finngen.gitbook.io/documentation/methods/phewas",
        "sleep_N_total": str(int(fg["num_cases"]) + int(fg["num_controls"])),
        "sleep_N_eff_formula": "4/(1/Ncase+1/Ncontrol)",
        "sleep_N_eff_derived": format(4/(1/int(fg["num_cases"])+1/int(fg["num_controls"])), ".12g"),
        "sleep_phenotype": "Clinical registry insomnia: F510|G470 in primary outpatient, hospital discharge or death registry; no minimum-event rule; controls exclude F5_SLEEP.",
        "sleep_controls_detail": "F5_SLEEP includes F51/3074/3064 registry rules plus F5_INSOMNIA,HYPERSOMNIA,SLEEPWAKE,SLEEPTERRORS,NIGHTMARES,SLEEP_NOS; do not import F5_SLEEP's own F5_BEHAVE control rule into F5_INSOMNIA.",
        "sleep_window": "Longitudinal linked national registers; precise endpoint lookback/end horizon and chronicity/spacing not certified by inspected R13 definition.",
        "sleep_effect_schema_primary": "ALT effect; BETA log(OR); SEBETA effect SE; PVAL regenie P; AF_ALT/case/control AF; GRCh38.",
        "sleep_schema_boundary": "Documented standard R13 schema; literal file/header/completeness/finite-row diagnostics not body-verified here.",
        "sleep_N_schema_boundary": "Endpoint case/control counts explicit in official manifest; no per-variant N in documented summary schema; constant effective N is a frozen assumption to audit.",
        "sleep_INFO_gate": "No INFO in documented GWAS schema; current source QC docs describe INFO-distribution checks without numeric release threshold. Per-variant >0.9 not certified; do not fill INFO=1 or borrow another release threshold.",
        "sleep_release_identity": "R13 regenie3.3; official public release timeline Q3 2026; metadata/body object created May2026 is not proof of public release date.",
        "sleep_remote_bytes": "809346932", "sleep_remote_MD5": "f16c21acbf8b9ebfeb0f26a19f0ecfe3", "sleep_remote_generation": "1777989563097164",
        "sleep_body_SHA256": "NOT_DOWNLOADED_NOT_ESTABLISHED", "sleep_h2_intercept_power": "NOT_ESTIMATED; source h2 Z>=4/intercept<=1.2 and prospective pair power remain required.",
        "outcome_accession": accession, "outcome_phenotype": r["replication_phenotype_definition"],
        "outcome_source_url": r["replication_source_url"], "outcome_cases": r["cases"], "outcome_controls": r["controls"],
        "outcome_primary_metadata_URLs": "https://www.ebi.ac.uk/gwas/rest/api/v2/studies/" + accession + ";" + r["replication_source_url"] + "-meta.yaml",
        "outcome_build": r["build"], "outcome_ancestry": r["ancestry"], "outcome_source_MD5": r["replication_checksum"],
        "outcome_remote_bytes": r["replication_content_length_bytes"], "outcome_remote_ETag": r["replication_etag"],
        "outcome_generation_identity": "GWAS_CATALOG_ACCESSION_AND_MD5; no GCS generation applicable",
        "outcome_historical_source_SHA256": historical_receipt["source_verification"]["observed_sha256"],
        "outcome_historical_processed_SHA256": historical_receipt["munged_output_sha256"],
        "outcome_uncertainty_field_mapping": "effect_allele/other_allele; odds_ratio; ci_lower/ci_upper; p_value; effect_allele_frequency; r2. Script47 Z=log(OR)/[(log(ci_upper)-log(ci_lower))/(2*1.959963984540054)], constant N_eff from cases/controls. Supplied standard_error is not required/used by script47; its outcome-file presence/missingness not inspected here.",
        "outcome_reused_wave": "Historical MVP GIA EUR public phenome-wide wave PMID39024449; previously processed outcome reused. No new MVP wave/participant independence from historical validation.",
        "outcome_historical_h2": oh["h2"], "outcome_historical_h2_SE": oh["h2_se"], "outcome_historical_h2_Z": oh["h2_z"],
        "outcome_historical_intercept": oh["LDSC_intercept"], "outcome_historical_h2_QC": oh["primary_status"],
        "outcome_uncertainty_gate": "Historical script47 reconstructs SE from log-CI width assuming95% Wald CI; outcome-specific effect/CI semantics need confirmation for new inferential admission. Historical h2 pass does not establish this.",
        "outcome_phenotype_gate": "Abdominal symptom group versus PanUKB same named PheCode; complete ascertainment/map not identical." if accession == "GCST90479148" else "History of drug allergy PheCode977 versus PanUKB ICD10Z88; cross-taxonomy exactness unverified; historical low-power QC fails.",
        "cohort_independence": "FinnGen, MVP and UKB are distinct named sources; individual cross-enrollment unknown. No certificate of zero intersections; avoid FinnGen sleep x FinnGen outcomes.",
        "data_rights": "Public R13 summary research reuse supported by official FinnGen access pages; official workflow describes online form/instructions, not submitted here. Cite/acknowledge. Explicit raw-source redistribution license not established. MVP EBI terms apply; keep source bodies out of Git.",
        "source_terms_URLs": "https://www.finngen.fi/en/researchers/accessing;https://finngen.gitbook.io/documentation/data-download;https://www.ebi.ac.uk/about/terms-of-use/",
        "source_terms_boundary": "EMBL-EBI imposes no additional restrictions beyond original owners; resource/owner-specific terms may govern. Public hosting is not an explicit original-owner raw-file redistribution certificate.",
        "scientific_scope_judgment": "Exact questionnaire equivalence is unnecessary for a prospectively declared symptom-to-registry transport test; not exact insomnia-disorder replication or a treatment/mechanism test.",
        "human_clinical_approval": "NOT_OBTAINED; clinical interpretation and any chronic-insomnia equivalence claim require qualified human review.",
        "metadata_acquisition_preprocessing_judgment": "CONDITIONAL_CANDIDATE: may freeze lawful bounded acquisition/preprocessing/source-h2 feasibility under explicit source-specific INFO/N limitations without pair outcome inspection; historical stages first. This does not admit rg.",
        "new_pair_inference_judgment": "CONDITIONAL_NOT_ADMITTED: outcome CI semantics, source-specific INFO/N handling, body identity/harmonization, sleep h2/intercept, power and217-family protocol gates outstanding." if accession == "GCST90479148" else "QC_INELIGIBLE: retain locked row; historical outcome h2Z<4; do not rescue by better sleep source or lower threshold.",
        "new_results": "NONE; no body download, new h2 or rg, significance or replication claim.",
    })
write_tsv("tables/sleep_clinician_finngen_MVP_admissibility_v4.tsv", admissibility)

new_sleep = [
    {"source_accession":"GCST90475826 / pha010373.1 / Phe_327_4.EUR.GIA", "cohort":"MVP European GIA", "clinical_construct":"Recorded EHR insomnia PheCode327.4", "question_device_rule":">=2 mapped ICD instances cases;0 instances controls; no related-PheCode exclusions; exact codes in v2 contract/CIPHER14926", "window":"EHR through September2020; minimum lookback/observation density and spacing/distinct-date semantics unknown", "effects_units":"Explicit effect_allele OR; historical processing takes log(OR); reference/ALT labels do not by themselves define effect allele", "uncertainty_fields":"ci_lower/ci_upper used to derive SE assuming95% log-Wald interval; supplied standard_error absent in v3 insomnia audit; p_value/r2/EAF present", "sample":"78566 cases;329572 controls;408138 total", "build":"GRCh38; unique pinned point liftover and HM3 coordinate/allele agreement independently verified in v3", "h2_reliability":"No source h2/intercept/power admitted; algorithm validation not reported in CIPHER; no chronic insomnia reliability coefficient verified", "overlap":"MVP named source distinct from UKB/FinnGen; individual intersection unknown; repeated MVP waves can share participants", "human_judgment":"Meaningful related clinical transport can be prespecified without equivalence; timing/chronicity/uncertainty and exact implementation qualifications retained", "status":"PUBLIC_GIA_DEFINITION_RESOLVED_QUALIFIED; existing sleep x FinnGen10 inference remains stopped", "primary_URLs":"https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90475826;https://phenomics.va.ornl.gov/web/api/phenotype/14926;https://pmc-oa-opendata.s3.amazonaws.com/PMC12857194.1/NIHMS2091849-supplement-Supplemental_Text.pdf"},
    {"source_accession":"FINNGEN_R13_F5_INSOMNIA", "cohort":"FinnGen R13 Finnish European", "clinical_construct":"Recorded registry insomnia", "question_device_rule":"F510|G470 primary outpatient/hospital/death; no minimum-event rule; controls exclude F5_SLEEP; not a symptom questionnaire", "window":"Longitudinal national register capture; exact endpoint start/lookback/end horizon and code prefix implementation not certified here", "effects_units":"ALT log-OR BETA; GRCh38 reference/ALT schema explicitly documented", "uncertainty_fields":"SEBETA documented effect SE; PVAL regenie P; AF_ALT/cases/controls documented; no INFO or per-variantN in standard schema; literal body not inspected", "sample":"51643 cases;446273 controls;497916 total; female cases36693/male14950 from endpoint metadata", "build":"GRCh38; no liftover/materialization performed", "h2_reliability":"No source h2/intercept/power estimated; no chronicity, severity or chart-adjudication reliability certified", "overlap":"Distinct named cohort from MVP/UKB, individual cross-enrollment unknown; later FinnGen releases overlap earlier FinnGen", "human_judgment":"Registry-to-symptom transport is legitimate scope for prespecification; human clinical signoff absent and exact chronic-insomnia equivalence unsupported", "status":"CONDITIONAL_METADATA_ACQUISITION_PREPROCESSING_SOURCE_H2_FEASIBILITY; not pair rg admitted", "primary_URLs":"https://r13.risteys.finngen.fi/endpoints/F5_INSOMNIA;https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv;https://finngen.gitbook.io/documentation/data-description"},
]
write_tsv("tables/sleep_clinician_new_sleep_crosswalk_v4.tsv", new_sleep)

questions = [
    ("CLINICAL_01", "Sleep clinician", "Do the usual UKB initiation/maintenance complaints and FinnGen/MVP recorded insomnia support a useful prespecified transport question while retaining different ascertainment and control definitions?", "Human scope/interpretation judgment; equivalence is not needed for an explicitly related-construct exploratory test.", "No approved clinical equivalence claim."),
    ("CLINICAL_02", "Sleep clinician", "Which statements can be made without chronicity, symptom-night frequency, adequate opportunity or daytime impairment in the source rules?", "Public definitions do not establish chronic insomnia disorder, severity, subtype or treatment relevance.", "Symptom liability and registry-coded insomnia labels only."),
    ("MVP_03", "MVP source/data expert", "For Phe_327_4.EUR.GIA/GCST90475826, were the >=2 codes on distinct dates or separated by a required interval; what minimum VA observation history was required?", "GIA occurrence thresholds/code list/September2020 horizon resolved in v2; timing/lookback/executable implementation remain unknown.", "Do not reopen all public GIA documentation as missing."),
    ("MVP_04", "MVP source/statistical expert", "For GCST90475826 and the reused outcome files GCST90479148/GCST90479330, are deposited CIs95% log-scale Wald intervals for the explicit effect allele, and how were EBI effects/SE converted?", "Script47 uses log-CI-width SE; no source-specific independent uncertainty certification. MVP sleep whole-file CI diagnostics do not certify outcome files.", "New confirmatory pair inference remains conditional."),
    ("MVP_05", "MVP source expert", "Which stage/release explains accession MAC<30 versus final supplement association MAC>40 and genotype MAC<20?", "Documented discrepancy; stage equivalence not established.", "Do not silently normalize MAC thresholds."),
    ("FINNGEN_06", "Finnish clinical coding expert", "For the exact R9 G6_SLEEPAPNO extraction, is ICD9 Finland3472/3472A an intended sleep-apnea historical code or an endpoint/code-table error, and how many unique analyzed cases depend only on it?", "R9 code table freshly confirms1574 code-related cases labelled cataplexy/narcolepsy; counts overlap.", "No pure OSA, misclassification fraction or subtype claim."),
    ("FINNGEN_07", "FinnGen source expert", "For exact R13 F5_INSOMNIA, what observation horizon/lookback, code-prefix expansion and association-stage INFO/MAC filters applied; is release-matched per-variant imputation quality retrievable?", "No minimum-event rule and exact code prefixes/control endpoint are public; standard GWAS schema has no INFO/N. Separate29.77GB annotation object metadata verified, contents not inspected.", "No >0.9 assumption or INFO=1 fill."),
    ("MEASUREMENT_08", "Sleep measurement expert", "Does the source's four-category sleepiness coding span the initial-assessment removal of All of the time; would the response change affect the proposed comparison?", "Official field-version note versus publication coding; no magnitude of bias assumed.", "Retain exact ordinal dozing construct."),
    ("MEASUREMENT_09", "Sleep measurement expert", "Are the chosen device algorithm, valid-night rules, main-period definitions and timing lag acceptable for the intended contrast with24h self-report?", "Core actigraphy shared UKB subset and numerator/window, transformed effects, quiet-wake limitation; no equivalent sampling or reliability established.", "No significance-count biological inference."),
    ("COHORT_10", "Cohort epidemiologist", "Can exact participant/intersection records or defensible source-specific overlap evidence establish UKB-FinnGen-MVP intersections and between-MVP-wave reuse?", "Distinct cohort names alone do not establish individual disjointness; historical outcome participants are reused.", "UNKNOWN_INTERSECTION until verified."),
    ("RIGHTS_11", "Investigator/data steward", "Does the investigator have the required FinnGen public-summary access instructions and applicable EBI rights, and what raw-source redistribution permission is explicit?", "Official free/public research access and acknowledgement/citation verified; no form/DUA signed or communication sent here.", "Research reuse support does not certify repository redistribution."),
]
write_tsv("tables/sleep_clinician_human_source_questions_v4.tsv", [dict(zip(["question_id", "qualified_reviewer", "exact_unsent_question", "current_primary_fact", "admission_consequence"], row)) for row in questions])

objective = []
for cohort,n,device,window in [
    ("Whitehall II", "2144", "GENEActiv", "Up to7 days"),
    ("CoLaus", "2257", "GENEActiv", "Up to14 days"),
    ("Rotterdam Study RS-I/II/III subset", "1418", "Actiwatch AW4", ">=96h valid recording; no week of daylight-saving change"),
]:
    objective.append({"candidate":cohort,"exact_primary_study":"Jones2019 PMID30952852 DOI10.1038/s41467-019-09576-1; non-UKB follow-up cohorts","source_url":"https://www.nature.com/articles/s41467-019-09576-1","clinical_relevance":"Related main-period device sleep duration/efficiency/timing; study methods apply common derivation/exclusions where available.","N":n,"device":device,"window":window,"GWAS_scope_verified":"Targeted follow-up of47 UKB-selected associations; public full-genome non-UKB summary source/accession not established.","effect_units":"Original and inverse-normal analyses where available; do not assume availability of exact core RAW_SIN release.","build_SE_P_AF_INFO_rights":"Full non-UKB source schema/build/alleles/SE/P/AF/INFO/rights not established; UKB public availability statement does not extend automatically.","participant_disjointness":"Non-UKB named source; no individual intersection certificate with UKB/MVP/FinnGen.","admissibility":"NOT_ADMITTED_FOR_LDSC; potentially useful investigator-request target only after full genome source/rights/ancestry/QC/power evidence; no request sent."})
objective.append({"candidate":"HCHS/SOL Sueño/MESA actigraphy resources","exact_primary_study":"NSRR study-specific device documentation; no exact eligible GWAS accession identified by bounded search","source_url":"https://sleepdata.org/datasets/hchs/pages/actigraphy-introduction.md;https://sleepdata.org/datasets/mesa/variables?folder=Sleep+and+Circadian+Studies%2FActigraphy%2FSleep-wake+patterns","clinical_relevance":"Main sleep duration, efficiency/midpoint potentially relevant; participant-level sleep data availability is not a genome-wide summary release.","N":"Sueño recruited2252; not GWAS N","device":"Actiwatch Spectrum documented for Sueño; MESA actigraphy definitions available","window":"Sueño acquisition2010-2013; exact release GWAS window unverified","GWAS_scope_verified":"No eligible public cohort-disjoint full summary GWAS certified.","effect_units":"Unknown GWAS transform/units","build_SE_P_AF_INFO_rights":"Unknown exact summary schema/ancestry/build/rights; multi-ancestry data require appropriate LD/reference rather than forced EUR application.","participant_disjointness":"Named non-UKB studies; exact cohort intersections unknown.","admissibility":"UNAVAILABLE_FOR_CURRENT_TEST; do not substitute questionnaire results, selected SNP tables, PRS scoring files or individual-level NSRR data for full GWAS."})
objective.append({"candidate":"2026 device sleep stage/duration GWAS","exact_primary_study":"Nature Communications DOI10.1038/s41467-026-71252-y","source_url":"https://www.nature.com/articles/s41467-026-71252-y","clinical_relevance":"New algorithm-derived REM/NREM/night sleep/efficiency; not EEG-defined PSG stages or exact core phenotypes.","N":"80013","device":"UKB wrist accelerometry with trained classification model","window":"UKB source; exact window not audited for admission","GWAS_scope_verified":"Primary source explicitly UKB; same cohort category as core actigraphy.","effect_units":"Not source-audited for prospective admission","build_SE_P_AF_INFO_rights":"Not pursued because disjoint-cohort prerequisite fails; no summary-body fetch","participant_disjointness":"UKB source reuses discovery cohort category; not independent objective sleep replication.","admissibility":"PARTIAL_SAMPLE_OVERLAP/PHENOTYPE_SENSITIVITY_ONLY; cannot repair fully independent sleep-side validation."})
write_tsv("tables/sleep_clinician_objective_source_screen_v4.tsv", objective)

paths = [core_path, window_path, h2_path, family_path, outcome_h2_path,
    "config/gwas_schemas.tsv", "sleep_unified_research_v2/source_definition_review/provenance_v2_definition_contract.json",
    "sleep_unified_research_v2/source_definition_review/provenance_v2_source_definition_addendum.md",
    "sleep_unified_research_v3/FROZEN_MVP_PREPROCESSING_PROTOCOL_v4.md",
    "sleep_unified_research_v3/tables/MVP_DOWNSTREAM_GATES_V3.tsv",
    "discovery_extension/scripts/47_stream_replication_sources.py",
    "discovery_extension/provenance/replication_streaming_receipts/gwas_catalog_GCST90479148.json",
    "discovery_extension/provenance/replication_streaming_receipts/gwas_catalog_GCST90479330.json",
    "work/separate_sleep_phenotyping_review/actigraphy.html"]
receipt = {"reviewer":"sleep_clinician_v4 automated independent-role review; not a licensed-human signoff", "created_utc":datetime.now(timezone.utc).isoformat(),
    "scope":"Clinical crosswalk and source admissibility only; no manuscript, treatment advice, inferred subtype/mechanism, contacts or commits.",
    "new_GWAS_body_downloads":0,"new_h2_runs":0,"new_rg_runs":0,"new_pair_results_inspected":False,
    "historical_tables_inspected":"Locked candidate metadata and existing h2 QC; initial broad historical ledger read includes historical discovery columns, not new test outcomes. Final candidate choice is restricted to already217-locked insomnia-MVP rows.",
    "old_packages_modified":False,"core_crosswalk_rows":len(crosswalk),"prospective_pair_metadata_rows":len(admissibility),
    "clinical_scope_equivalence_required":False,"human_clinical_approval_obtained":False,
    "new_inference_admitted":False,"MVP_public_GIA_definition_gate":"RESOLVED_AT_DOCUMENTED_DEFINITION_LEVEL; implementation timing/uncertainty qualifications retained",
    "input_artifacts":[{"path":p,"bytes":(ROOT/p).stat().st_size,"sha256":hashlib.sha256((ROOT/p).read_bytes()).hexdigest()} for p in paths]}
(V4/"reviews/sleep_clinician_review_receipt_v4.json").write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"crosswalk_rows":len(crosswalk),"source_candidate_rows":len(admissibility),"questions":len(questions),"objective_screen_rows":len(objective),"new_inference_admitted":False},indent=2))
