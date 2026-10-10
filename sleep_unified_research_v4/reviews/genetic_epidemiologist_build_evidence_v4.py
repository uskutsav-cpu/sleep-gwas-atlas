#!/usr/bin/env python3
"""Bounded metadata-only genetic epidemiology review builder.

No GWAS body, estimator, shell subprocess, original output, or pair rg is read
for source eligibility/power planning. Existing counts are inherited explicitly.
"""
from pathlib import Path
from datetime import datetime, timezone
from statistics import NormalDist
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parents[2]
V4 = ROOT / "sleep_unified_research_v4"
PREFIX = "genetic_epidemiologist_"
INPUTS = [
    "config/analysis_panel.tsv",
    "sleep_unified_research_v1/FINAL_EVIDENCE_HANDOFF.md",
    "sleep_unified_research_v1/FROZEN_NEW_ANALYSIS_PROTOCOL.md",
    "sleep_unified_research_v1/FROZEN_NEW_ANALYSIS_PROTOCOL_v2.md",
    "sleep_unified_research_v1/SHARED_ESTIMATOR_COVARIANCE.md",
    "sleep_unified_research_v1/reviews/cohort_independence_v1.md",
    "sleep_unified_research_v1/reviews/cohort_independence_v1.json",
    "sleep_unified_research_v1/native/core_pilot_v1/rg_insomnia__bmi.full_precision.json",
    "sleep_unified_research_v1/tables/original_core_396.tsv",
    "sleep_unified_research_v1/tables/original_extension_1200.tsv",
    "sleep_unified_research_v1/REPLICATION_RESULTS.tsv",
    "sleep_unified_research_v2/FROZEN_NEW_ANALYSIS_PROTOCOL_v3.md",
    "sleep_unified_research_v2/source_definition_review/provenance_v2_definition_contract.json",
    "sleep_unified_research_v3/tables/MVP_DOWNSTREAM_GATES_V3.tsv",
    "sleep_unified_research_v3/tables/MVP_RETAINED_P_CI_DIAGNOSTICS.tsv",
    "sleep_unified_research_v3/tables/MVP_WHOLE_SOURCE_CI_CENTERING.tsv",
    "sleep_unified_research_v4/reviews/sleep_geneticist_hypothesis_assessment_v4.md",
    "sleep_unified_research_v4/reviews/sleep_geneticist_primary_admission_matrix_v4.tsv",
    "sleep_unified_research_v4/reviews/statistical_geneticist_critical_review_v4.md",
    "sleep_unified_research_v4/reviews/statistical_geneticist_critical_review_v4.json",
    "sleep_unified_research_v4/reviews/genomicsem_specialist_independent_review_v4.md",
    "sleep_unified_research_v4/reviews/independent_reproducibility_audit_v4.md",
    "sleep_unified_research_v4/statistical_validation/processed_N_coding_diagnostics_v4.tsv",
    "sleep_unified_research_v4/statistical_validation/liability_input_sensitivity_v4.tsv",
    "sleep_unified_research_v4/statistical_validation/statistical_claim_admission_v4.tsv",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_finngen_R13_manifest_selection_v4.json",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_R13_insomnia_primary_excerpt_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_R13_sleep_exclusion_primary_excerpt_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_R13_data_schema_primary_excerpt_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_R13_GWAS_methods_primary_metadata_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_R13_imputation_QC_primary_metadata_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_finngen_R13_insomnia_HEAD_v4.json",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_GCST90479148_metadata_primary_metadata_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_GCST90479148_yaml_primary_metadata_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_GCST90479330_metadata_primary_metadata_v4.txt",
    "sleep_unified_research_v4/source_provenance/sleep_clinician_GCST90479330_yaml_primary_metadata_v4.txt",
    "discovery_extension/results/replication/replication_source_h2.tsv",
    "sleep_unified_research_v4/source_provenance/genetic_epidemiologist_live_primary_access_v4.json",
    "sleep_unified_research_v4/reviews/sleep_clinician_review_v4.md",
    "sleep_unified_research_v4/tables/sleep_clinician_finngen_MVP_admissibility_v4.tsv",
    "sleep_unified_research_v4/tables/sleep_clinician_objective_source_screen_v4.tsv",
]

def seal(path):
    data = path.read_bytes()
    return {"path": str(path), "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}

def read_tsv(rel):
    with (ROOT / rel).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(area, stem, rows):
    path = V4 / area / (PREFIX + stem + "_v4.tsv")
    assert path.parent in {V4 / "reviews", V4 / "source_provenance", V4 / "statistical_validation"}
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path

before = [seal(ROOT / rel) for rel in INPUTS]
panel = read_tsv("config/analysis_panel.tsv")
definition = json.loads((ROOT / INPUTS[12]).read_text())
manifest = json.loads((ROOT / INPUTS[25]).read_text())
fg = next(x for x in manifest["selected_rows"] if x["phenocode"] == "F5_INSOMNIA")
outcome_h2 = {r["study_accession"]: r for r in read_tsv(INPUTS[36])}
normal = NormalDist()
alpha = .05 / 217
zcrit = normal.inv_cdf(1 - alpha / 2)

def power(effect, se):
    ncp = abs(effect) / se
    return 1 - normal.cdf(zcrit - ncp) + normal.cdf(-zcrit - ncp)

lo, hi = .0001, 1.0
for _ in range(100):
    mid = (lo + hi) / 2
    if power(.15, mid) >= .80:
        lo = mid
    else:
        hi = mid
se_ceiling = (lo + hi) / 2
power_rows = [{
    "purpose": "ILLUSTRATIVE_PLANNING_NOT_OBSERVED_POWER", "family_size": 217,
    "two_sided_alpha": format(alpha, ".17g"), "target_abs_rg": .15,
    "assumed_calibrated_SE": se, "two_sided_normal_power": format(power(.15, se), ".17g"),
    "MDE_80_normal_approx": format((zcrit + normal.inv_cdf(.8)) * se, ".17g"),
    "assumptions": "Known calibrated SE; normal estimator; fixed outcome-blind alternative; no conditioning or selection correction inferred",
} for se in [.02, .03, .04, .05, .08, .10]]
power_path = write_tsv("statistical_validation", "prospective_precision_grid", power_rows)

n_rows = []
def add_n(sid, cohort, phenotype, cases, controls, total, provenance, semantics, phenotype_key=None):
    binary = cases is not None
    assert total > 0 and (not binary or cases + controls == total)
    neff = 4 * cases * controls / total if binary else None
    n_rows.append({"phenotype_key": phenotype_key or sid, "source_id": sid, "cohort_population": cohort, "phenotype": phenotype,
        "metadata_total_N": total, "metadata_cases": cases if binary else "NA",
        "metadata_controls": controls if binary else "NA",
        "case_fraction": format(cases / total, ".17g") if binary else "NA",
        "derived_balanced_Neff": format(neff, ".17g") if binary else "NA",
        "Neff_fraction_total": format(neff / total, ".17g") if binary else "NA",
        "metadata_evidence": provenance, "estimator_N_semantics": semantics,
        "interpretation": "Metadata-only count arithmetic; not per-SNP N, h2, liability validation, power or new association"})

for row in panel:
    if row["domain"] != "sleep":
        continue
    cases = None if row["ncase"] == "NA" else int(row["ncase"])
    controls = None if cases is None else int(row["ncontrol"])
    add_n(row["source_id"], "Finnish EUR" if row["trait_id"] == "sleep_apnea" else "UKB EUR",
        row["phenotype_definition"], cases, controls, int(row["n_total"]),
        "config/analysis_panel.tsv; immutable historical metadata",
        "Actual processed per-SNP N not reread; derived Neff is not an estimator-input certification", phenotype_key=row["trait_id"])
add_n(definition["gwas_accession"], "MVP EUR GIA", definition["ascertainment"],
    definition["n_cases"], definition["n_controls"], definition["n_total"], INPUTS[12],
    "v3 frozen constant effective-N rule; source effect/CI calibration unresolved; body not reread")
add_n("FINNGEN_R13_F5_INSOMNIA", "Finnish EUR", "Registry F51.0/G47.0; F5_SLEEP-excluding controls",
    int(fg["num_cases"]), int(fg["num_controls"]), int(fg["num_cases"]) + int(fg["num_controls"]), INPUTS[25],
    "No per-variant N in documented schema; constant effective N requires explicit prospective assumption")
for accession, c, k in [("GCST90479148", 89262, 299053), ("GCST90479330", 535, 448351)]:
    rel = f"sleep_unified_research_v4/source_provenance/sleep_clinician_{accession}_metadata_primary_metadata_v4.txt"
    metadata = json.loads((ROOT / rel).read_text())
    assert metadata["accession_id"] == accession
    assert f"{c:,}" in metadata["initial_sample_size"] and f"{k:,}" in metadata["initial_sample_size"]
    add_n(accession, "MVP EUR", metadata["disease_trait"], c, k, c + k, rel,
        "Historical outcome materializer semantics and CI-derived Z need accession-specific review; body not reread")
n_path = write_tsv("source_provenance", "metadata_N_semantics", n_rows)

design_specs = [
    ("HISTORICAL_41", "ONE_SIDE_EXTERNAL_VALIDATION", "Discovery sleep reused in every estimate", "Named validation outcomes distinct; individual intersections unknown", "Historical 217 selection and .05/217 rule", "Current qualified labels; no both-trait replication or calibrated heterogeneity"),
    ("MVP_INSOMNIA_FINNGEN_10", "PHENOTYPE_SENSITIVITY_CONDITIONAL_TRANSPORT", "Both discovery named sources changed; FinnGen outcomes reused from prior validation", "Four cross-stage named-source distinctions; unknown persons/relatives", "Frozen source-specific v1/v2/v3 amendments; keep ten candidates and 217", "CI effect/test semantics; source QC; pair precision; clinical transport and overlap qualification"),
    ("FINNGEN_INSOMNIA_MVP_ABDOMINAL", "PHENOTYPE_SENSITIVITY_FEASIBILITY_ONLY", "Both discovery named sources changed; MVP outcome reused from historical validation", "FinnGen sleep overlaps broader FinnGen validation wave; exact cross-cohort intersections unknown", "NOT_ALREADY_SOURCE_SPECIFIC_FROZEN; alternate-source multiplicity amendment required", "Sleep body/INFO/N; outcome CI semantics; source h2/intercept; power; source identity; amendment"),
    ("FINNGEN_INSOMNIA_MVP_ALLERGY", "QC_INELIGIBLE_HISTORICAL_OUTCOME", "Both discovery named sources changed; MVP outcome reused", "Same cohort/source risks as abdominal; no pair estimated here", "New source-specific amendment needed if reconsidered; do not select another favorable release", "Outcome h2 Z below 4; retain unestimated status, not negative control or null"),
    ("R9_APNEA_R13_OUTCOME", "PARTIAL_SAMPLE_OVERLAP_SOURCE_DESIGN", "FinnGen release/source reused", "R9/R13 membership intersection unmeasured; no strict nesting assumption", "Frozen exclusion preserved", "Explicit disjoint sample evidence required for independent claim"),
    ("UKB_DEVICE_SELF_REPORT", "MEASUREMENT_CONTRAST_NOT_IDENTIFIED_AS_PURE_MODALITY", "Device subset, shared disease outcome and UKB samples", "Construct/time/selection differences; exact sample intersections unknown", "New target/contrast family not frozen", "Full joint uncertainty; target phenotype rationale; independent device source; prospective contrast power"),
    ("JOINT_BMI_SMOKING", "NEW_PRIMARY_NO_GO_CURRENT", "Exact smoking construct/outcome independent source not frozen", "Exposure meta-analysis intersections with target/source unknown", "No primary admitted; freeze before outcomes", "Clinical target; source identities; complete S/V; residual stability; power; exact novelty gap; independent test"),
]
design_rows = []
for d in design_specs:
    row = dict(zip(["design_id", "current_class", "source_reuse", "intersection_boundary", "family_status", "remaining_gate"], d))
    row["new_pair_outcomes_accessed"] = False
    row["admitted_by_this_review"] = False
    design_rows.append(row)
design_path = write_tsv("reviews", "source_design_admission", design_rows)

claim_specs = [
    ("C01", "Original 396 and 1200 families/counts as dated historical findings", "ADMITTED_HISTORICAL_QUALIFIED", "Current fully reproduced discovery total", "Native full-precision comparisons and all scientific QC decisions remain separate"),
    ("C02", "23 qualified outcome-side external positives among 217", "ADMITTED_QUALIFIED", "23 independent both-trait replications", "Discovery sleep input reused; zero fully independent two-trait successes"),
    ("C03", "Distinct named-cohort source design", "ADMITTED_WITH_PERSON_UNKNOWN", "Proven zero/negligible participant or relative intersection", "Membership/linkage evidence or justified bound; statistical intercept alone insufficient"),
    ("C04", "R9/R13 overlapping source risk", "ADMITTED_DESIGN_BOUNDARY", "New release proves independent replication", "Explicit disjoint participants/relatedness and phenotype/source contract"),
    ("C05", "Source-specific common-variant rg sharing", "QUALIFIED_ESTIMAND", "Sleep causes disease or predicts clinical benefit", "Causal/intervention/risk estimands are not supplied by global rg"),
    ("C06", "Clinical-versus-symptom related-phenotype transport", "CONDITIONAL_NOT_ESTIMATED", "Exact chronic-insomnia replication", "Clinical adjudication and compatible construct/window/control rules not established"),
    ("C07", "GIA documented two-instance case rule resolved", "ADMITTED_DOCUMENTED_RULE", "Distinct-day, chronicity and historical executable implementation certified", "Exact extraction/lookback/spacing evidence needed"),
    ("C08", "Metadata-derived balanced case/control N", "ADMITTED_ARITHMETIC_ONLY", "Total N proves source or pair power", "SNP information, h2, uncertainty/overlap and source-specific N calibration required"),
    ("C09", "MVP receipt-verified preprocessing", "ADMITTED_PREPROCESSING_ONLY", "Upstream CI-derived SE is calibrated", "Accession CI level/test/scale and EBI effect conversion evidence required"),
    ("C10", "FinnGen explicit SE documented", "ADMITTED_SCHEMA_ONLY", "F5 insomnia satisfies inherited INFO/N/raw eligibility", "Literal body/finite rows, INFO policy and N semantics unpassed; source QC != per-variant INFO"),
    ("C11", "MVP allergy source has historical h2 failure", "ADMITTED_QC_FAILURE", "Absent or failed pair proves no sharing", "535 cases; unchanged h2 Z .706 fails; no artificial negative control"),
    ("C12", "Lipid intercept/control-flow diagnostics", "DIAGNOSTIC_ONLY", "Resolved lipid disease interpretation", "All 36 flagged original pair fits need source/LD/model adjudication; native effect unestimated"),
    ("C13", "MS/melanoma fixed-fit liability sensitivity", "ADMITTED_ARITHMETIC_ONLY", "Validated liability h2 or restored melanoma eligibility", "Prevalence/ascertainment/N convention unresolved; melanoma Z unchanged below 4"),
    ("C14", "Device/report construct-specific descriptions", "ADMITTED_QUALIFIED", "Significance counts prove modality-specific biology", "Joint contrast variance, selection/measurement comparability and prospective power required"),
    ("C15", "Selected validation magnitudes/directions are descriptive", "ADMITTED_DESCRIPTIVE", "Independent binomial replication or explained winner's-curse attenuation", "Correlated selection, shared sleep/denominators and ascertainment confound such inference"),
    ("C16", "Covariance-aware difference is required", "ADMITTED_REQUIREMENT", "Calibrated discovery/validation heterogeneity", "Shared genomic boundaries, nonlinear denominator uncertainty and selection remain unpassed"),
    ("C17", "Fixed .15/.05/217/.8 precision target", "ADMITTED_PLANNING_ASSUMPTIONS", "Achieved power or clinical-effect threshold", "Calibrated source/pair planning SE and investigator meaningful-effect decision required"),
    ("C18", "Eligible imprecise estimate is inconclusive", "ADMITTED_INTERPRETATION_RULE", "Nonsignificance establishes absence or equivalence", "Prespecified equivalence/bounded-effect interval and calibrated CI required"),
    ("C19", "Reverse clinical-source route is a feasibility lead", "NO_TEST_ADMISSION", "Old 217 denominator automatically covers alternate-source/model searches", "Freeze source-specific amendment, all attempts, primary choice/alpha allocation before outcomes"),
    ("C20", "Joint BMI/personal-smoking projection is a future lead", "NO_NEW_PRIMARY", "Residualization establishes mediation or novel mechanism", "Exact constructs, target, full uncertainty, power, independent test and novelty gap required"),
    ("C21", "European/Finnish source-specific research scope", "ADMITTED_BOUNDARY", "Universal transport across ancestries, ages or sexes", "No matching independently powered target-population/sex-specific validation established"),
    ("C22", "Automated review package is complete for this perspective", "ADMITTED_REVIEW_ONLY", "Human manuscript/journal clearance", "Qualified human decisions and complete overall native/scientific readiness gates required"),
    ("C23", "Named separate frozen discovery/test families", "ADMITTED_FAMILY_BOUNDARY", "603 independent biological confirmations or union-wide FDR guarantee", "Correlated/nested phenotypes, source selection and precision matter; preserve original corrections and define new inferential families prospectively"),
]
claim_path = write_tsv("reviews", "claim_admission_matrix", [dict(zip(
    ["claim_id", "bounded_claim", "admission", "stronger_claim_blocked", "exact_remaining_requirement"], row)) for row in claim_specs])

questions = [
    ("GE01", "Cohort/genetic epidemiology reviewer", "For exact Ds/Do/Vs/Vo releases, what membership or defensible quantitative bound supports participant/relative intersection, and which remain unknown?", "Source-role matrix including within-stage and prior-wave reuse; retain qualified unknown where no evidence exists", "No proven participant-independent label"),
    ("GE02", "MVP source/statistical reviewer", "For GCST90475826 and GCST90479148/GCST90479330 individually, what are CI level/construction, tested allele, effect scale, P test and EBI conversion contract?", "Accession-specific upstream/author contract; no imported neighboring-accession semantics or result-selected P/CI filters", "No new h2/rg admission from preprocessing alone"),
    ("GE03", "FinnGen/source and statistical reviewer", "Can R13 F5_INSOMNIA satisfy inherited INFO>0.9 and per-variant-N semantics using documented evidence, or is a new prospective source policy scientifically necessary?", "Release-specific evidence or explicit result-free amendment; no INFO=1 or borrowed release threshold", "No body/INFO/N scientific eligibility asserted"),
    ("GE04", "MVP phenotype reviewer", "Are two ICD instances required on distinct days; what spacing, start lookback, observation-duration and executable historical extraction rules apply?", "Exact GIA/CIPHER historical implementation/window evidence; chronicity not inferred", "Documented clinical rule only, not adjudicated chronic-insomnia phenotype"),
    ("GE05", "Sleep clinician/epidemiologist", "Is symptom-to-MVP/FinnGen clinical transport a useful disease-specific question, given control exclusions and health-care ascertainment?", "Signed exact phenotype/window/control comparison and qualified clinical target", "No exact questionnaire/chronic-insomnia or clinical-actionability claim"),
    ("GE06", "Finnish epidemiologist/statistical geneticist", "What age/time/population target and matching prevalence support MS/melanoma liability scale and melanoma all-cancer-excluding controls?", "Applicable epidemiologic source and coherent N/P convention approved prospectively", "Fixed-fit rescaling remains diagnostic; melanoma h2-Z exclusion preserved"),
    ("GE07", "Investigators/statistical geneticist", "If reverse FinnGen sleep sources or multiple releases/models are attempted, what is the primary-source selection rule and complete test/within-hypothesis alpha allocation?", "Source-specific protocol before outcomes; all attempts/stopping/failures retained; historical families unchanged", "No automatic fresh source test at historical alpha"),
    ("GE08", "Clinical investigator/statistical geneticist", "What smallest meaningful rg/residual/contrast and precision/equivalence margin answer the intended research question?", "Outcome-blind justification, estimand-specific planning uncertainty and correction; .15 is not a patient benefit threshold", "No clinically meaningful null or power certificate"),
    ("GE09", "Covariance/statistical geneticist", "Does the proposed contrast include complete h2-denominator covariance, verified common genomic deletions and an explicit selection-aware interpretation?", "Independent calibration and fixed contrast/multiplicity; no index alignment or Cov=0 assumed", "No calibrated heterogeneity or modality contrast"),
    ("GE10", "Sleep/source specialist and investigator", "For joint BMI/smoking sharing, which exact own-smoking construct, source release and clinically justified outcome are fixed before results, with a distinct prior-art gap?", "Exact four-trait estimand, stable residuals/full covariance, outcome-blind power and independent design", "No new confirmatory primary or mediation claim"),
]
question_path = write_tsv("reviews", "human_questions", [dict(zip(
    ["question_id", "decision_owner", "exact_question", "required_evidence", "until_resolved"], row)) for row in questions])

after = [seal(ROOT / rel) for rel in INPUTS]
checks = [
    {"name": "all_small_review_inputs_unchanged", "pass": before == after},
    {"name": "positive_Neff_not_above_total", "pass": all(float(r["derived_balanced_Neff"]) <= r["metadata_total_N"] for r in n_rows if r["derived_balanced_Neff"] != "NA")},
    {"name": "N_table_has_12_core_sleep_and_4_clinical_sources", "pass": len(n_rows) == 16},
    {"name": "exact_normal_power_at_SE_ceiling", "pass": math.isclose(power(.15, se_ceiling), .8, abs_tol=1e-12)},
    {"name": "power_decreases_with_SE", "pass": all(float(a["two_sided_normal_power"]) > float(b["two_sided_normal_power"]) for a,b in zip(power_rows,power_rows[1:]))},
    {"name": "approximate_and_exact_SE_ceiling_agree", "pass": math.isclose(se_ceiling, .15 / (zcrit + normal.inv_cdf(.8)), rel_tol=1e-8)},
    {"name": "MVP_allergy_historical_h2_gate_fails", "pass": float(outcome_h2["GCST90479330"]["h2_z"]) < 4 and outcome_h2["GCST90479330"]["primary_status"] == "FAIL"},
    {"name": "MVP_abdominal_historical_QC_pass_not_new_source_admission", "pass": outcome_h2["GCST90479148"]["primary_status"] == "PASS"},
    {"name": "all_source_designs_unadmitted", "pass": not any(r["admitted_by_this_review"] for r in design_rows)},
]
assert all(c["pass"] for c in checks), checks
outputs = [power_path, n_path, design_path, claim_path, question_path,
    Path(__file__), V4 / "reviews/genetic_epidemiologist_critical_review_v4.md",
    V4 / "source_provenance/genetic_epidemiologist_live_primary_access_v4.json"]
receipt = {
    "schema": "genetic_epidemiologist_review_v4.1", "completed_utc": datetime.now(timezone.utc).isoformat(),
    "verdict": "QUALIFIED_RESOURCE_AND_OUTCOME_SIDE_EVIDENCE_ONLY_NO_NEW_PRIMARY_OR_TWO_TRAIT_REPLICATION",
    "inherited_historical_counts_not_reaudited": {"core_family":396, "core_primary_positive":153, "core_sensitivity_positive":8,
        "extension_family":1200, "extension_FDR_positive":603, "validation_family":217,
        "validation_categories":[23,18,17,159], "validation_pairs_estimated":41, "fully_independent_two_trait_replication":0},
    "planning_only": {"alpha":alpha, "abs_rg":.15, "power":.8, "zcrit":zcrit,
        "exact_normal_SE_ceiling":se_ceiling, "MDE_multiplier_80_normal_approx":zcrit+normal.inv_cdf(.8),
        "achieved_source_or_pair_power":None, "clinical_meaningful_effect_certified":False},
    "inputs_before": before, "inputs_after": after, "outputs": [seal(p) for p in outputs], "checks": checks,
    "no_GWAS_bodies_read_or_downloaded":True, "no_native_estimators_launched":True,
    "no_source_or_Git_mutations":True, "no_new_pair_outcomes_read":True,
    "new_review_arithmetic_only":True, "automated_review_not_human_approval":True,
    "no_test_family_amendment_implemented":True,
    "campaign_status": "Live separately; this static perspective does not certify later native-job completion",
}
receipt_path = V4 / "reviews/genetic_epidemiologist_evidence_receipt_v4.json"
receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
hash_path = V4 / "reviews/genetic_epidemiologist_v4_SHA256SUMS"
hash_path.write_text("".join(f"{seal(p)['sha256']}  {p.relative_to(ROOT)}\n" for p in outputs + [receipt_path]))
print(json.dumps({"checks_passed":len(checks), "inputs_unchanged":len(before), "SE_ceiling":se_ceiling,
    "outputs":len(outputs)+2, "receipt":str(receipt_path)}, indent=2))
