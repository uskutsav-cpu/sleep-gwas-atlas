# Draft research text for an NCDS genetic-data application

**Status:** Working draft, not an application, not signed and not submitted. The [current CLS application form and instructions](https://cls.ucl.ac.uk/data-access-training/data-access/accessing-data-directly-from-cls/) require an editable Word application and a completed CLS variable dictionary. The latest form should be downloaded before submission. If the applicant is a student, the form requires the project supervisor to be the main applicant. This draft deliberately leaves institutional facts and commitments for the research team to supply.

## Proposed title

Genetic susceptibility to insomnia and frailty and physical function in the 1958 British birth cohort

## Plain-language abstract

Sleep problems and frailty often occur together, but it is unclear whether their genetic influences also overlap. We will examine research-only genetic scores previously developed from published studies, using data from the 1958 British birth cohort. We plan to study whether these scores are associated with participants' reported sleep problems, a measure of accumulated health problems at age 50, and measured hand strength and walking speed in later adulthood. We will compare findings using different health measures and carefully account for questions about sleep that appear within some frailty measures. We will assess whether the scores add information beyond basic demographic and genetic-background factors and report the uncertainty of our findings. The scores are experimental and will not be used to make predictions or health decisions for participants.

## Keywords

National Child Development Study; insomnia; sleep; frailty index; physical function; genetics; polygenic score.

## Research purpose and methods

The work evaluates fixed, previously selected 15-marker and 22-marker research scores (with separate alternative weights and an equal-weight count) against measured outcomes. The primary score/outcome pair and multiplicity strategy need agreement by the research team **before** examining target associations. Define the 34-deficit age-50 FI from documented items if the original construction can be reproduced, and specify sleep- and fatigue-excluded exploratory versions with their own denominators and missingness thresholds. Analyze sleep-initiation and maintenance symptoms separately. For the age-62 visit, assess grip and walking speed as separate physical-function outcomes; they are not the full Fried frailty phenotype. If a reproducible physical-frailty phenotype and adequate events become available, specify it separately. Repeated sleep questions permit a further prespecified longitudinal analysis.

Request the combined QC-filtered and imputed NCDS genomic data and the minimum survey, demographic, ancestry and technical fields needed for the agreed outcomes and covariates. Confirm with CLS that age-62 measures can be included in the project-specific linked extract. The [variable inventory](ncds_genetic_sharing_dictionary_audit.md), [FI field map](ncds_age50_frailty_candidate_questionnaire_map.md) and [age-62 field review](ncds_age62_grip_gait_genetic_validation_lead.md) document candidate requests; they are not a complete validated phenotype or the required marked CLS workbook. [Marker request](ncds_score_variant_request.md) gives provisional genome-build mappings and separate counted-allele conventions.

**Working field inventory:** `ncds_request/candidate_variable_inventory.xlsx` collects 54 candidate names into 55 rows matched against the current v11 CLS dictionary, including the six unresolved composite domains and a sleep-maintenance question. `N8HLTHGN` appears in both `ncds8_imagine.sav` and `ncds_2008_followup.sav`; both are shown so that the appropriate dataset can be confirmed. This compact inventory is **not** CLS's completed application dictionary: it lacks confirmed role codes, final FI logic, demographic/technical variables and potentially relevant later-wave fields. The team must review and mark the current official dictionary, rather than submit the inventory as if it were the required form.

Before analyses, assess sample and relative overlap with all discovery and candidate-selection sources; exclude overlap where feasible or report that it cannot be excluded. Check genome build, exact alleles, multiallelic representation, imputation quality, missingness, ancestry and relatedness. Freeze the scoring rules, document variants unavailable after QC, and do not call a truncated score the complete prototype. Fit the same predeclared models and covariates for each score in the same eligible individuals; report sample flow, coefficients, confidence intervals and model assumptions. Where event counts and outcome type permit, report incremental predictive performance and calibration with uncertainty, without claiming that a retrospective score association is clinical prediction. Consider study nonresponse, survivorship and timing, and do not interpret associations as causal.

## AI and data handling — to be finalized by the team

The application asks whether AI tools will be used and requests tool, technical and methodological details. The researchers must specify their actual intended use, the applicable CLS policy and institutionally approved environment. This chat has used public documentation and aggregate evidence only; that does not authorize transmitting CLS participant records or restricted outputs to an external AI service. If machine learning or AlphaGenome is proposed, describe the intended inputs, purpose, validation and approved execution environment; do not make commitments about data handling or permission until verified with CLS and the institution. The core score evaluation can be conducted using conventional models.

The [CLS Research Data Access Framework, version 7, section 7.6](https://cls.ucl.ac.uk/wp-content/uploads/2017/02/CLS_Data_Access_Framework.pdf), requires disclosure of proposed AI use and prohibits processing potentially disclosive individual-level CLS data using cloud AI without explicit written CLS agreement. It permits AI-assisted coding using data-free code, public metadata, and error output without participant values; local and trusted-environment AI use is subject to the specific controls in that section. The team must check the current policy, the project's data agreement and its institutional arrangements before choosing any AI workflow.

## Required team-specific information before an application

- Main applicant/supervisor, affiliations, institutional emails, collaborators and access needs.
- Approved organizational computing, storage and security arrangements; required data-sharing agreement and training.
- Ethics status and institutional determination, incidental-findings plan and any clinical expert, conflicts of interest, funding and proposed dates. Do not infer or fabricate any of these.
- Final precise variable request in the current CLS dictionary, including how to request age-62 measures absent from version 11, and whether the genotyped-sample option and linked biometrics are permitted.
- The exact primary outcome, score, analysis plan, and intended AI use; check whether genotype and phenotype data can be used for that project and whether any data-management costs apply.

The CLS form states its plain-language abstract is published online. Its source and instructions are [available here](https://cls.ucl.ac.uk/wp-content/uploads/2017/02/CLS_data_access_guidelines_and_application.docx); the public copy was inspected locally. No form, dictionary or application has been submitted and no participant-level NCDS data were accessed.
