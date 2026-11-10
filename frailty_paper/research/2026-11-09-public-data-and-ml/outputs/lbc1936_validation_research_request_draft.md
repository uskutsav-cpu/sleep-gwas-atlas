# Draft research proposal: sleep and frailty genetic-score evaluation in LBC1936

**Working draft; not submitted.** Complete investigator, institution, ethics and access details only with accurate information. Check current EGA and LBC policies; no authorization or no-fee status is assumed.

## Provisional title

Evaluation of exploratory insomnia–frailty genetic scores against sleep symptoms and physical frailty in the Lothian Birth Cohort 1936

## Scientific rationale

Insomnia and frailty share some genetic associations, but existing frailty indices can include sleep or fatigue items. Physical frailty and a deficit-accumulation index measure overlapping, distinct constructs. Our candidate genetic sums were selected and weighted using prior GWAS results and have not been independently validated. LBC1936 has previously described Wave 3 sleep questions and a Fried-type physical frailty phenotype, along with genome-wide genetics. This study would evaluate the scores in a cohort with observed outcomes, rather than infer patient prediction from marginal GWAS associations.

## Proposed data and eligibility

Request permission to link eligible LBC1936 genetic data—including the joint WGS call release [EGAD00001005478](https://www.ega-archive.org/datasets/EGAD00001005478), or a suitable current genotyped/imputed alternative—to authorized Wave 3 phenotype and covariate data. Confirm whether a project-specific request can combine them and whether any access, transfer, computing or other fees apply. The [candidate marker crosswalk](ncds_score_variant_request.csv) is a reusable project-wide reference of the 33 distinct variants across three score files, **not** an LBC marker inventory; exact assembly, alleles and source contrast must be verified against the approved release.

Request subject-level or approved analysis-environment access, as applicable, to:

- Wave 3 sleep quality, usual nocturnal awakening, weeknight/weekend sleep latency and duration, recorded response codes and dates; do not equate these single symptoms with clinical insomnia.
- Published Wave 3 Fried-type physical frailty status and all available component values and coding, including the walking, grip, activity, exhaustion and weight-change inputs, plus prior-wave measures needed for the published algorithm.
- A deficit FI and item-level definitions, including identifiable sleep and fatigue items, where available. Derive any sleep-excluded or sleep-and-fatigue-excluded sensitivity indices only where item content, scoring and data rights are confirmed; these modified indices are not automatically validated as equivalent to the published FI.
- Genotype data, sample/variant QC and variant-specific imputation or sequencing quality, genome build, ploidy, multiallelic representation and exact allele contrasts; sample and variant missingness; ancestry PCs and relatedness information where permitted.
- Age at assessment, sex and relevant clinical covariates with documentation; dates and missingness necessary for a cohort flow diagram and permitted small-cell-safe summaries.

The public [variable crosswalk](lbc1936_dictionary_wave_crosswalk.csv) lists candidate Wave 3 names; field coding and availability require authoritative confirmation. The 2017 paper's 85 frail participants reflect a different complete-case analysis and must not be assumed to be our event count. Verify sample intersection, attrition, ancestry, score coverage and relatedness before assessing feasibility.

## Proposed analysis

Before inspecting outcomes, lock each evaluated score definition and its intended primary outcome; keep the 15-marker three-weight prototype, 22-marker equal-weight count, and 22-marker three-weight prototype separate. Do not treat the different weight sets as fitted joint-outcome scores. Confirm the FFS source-effect convention and every biallelic/multiallelic contrast before scoring. Document all prior exploratory selection, including use of MVP insomnia, FinnGen hospital frailty and the CHARGE gait results; evaluate whether this cohort contributed participants or relatives to any discovery or selection source. If independence cannot be established, describe the evaluation as exploratory and qualify any replication claim.

Primary analysis: a prespecified Wave 3 Fried-type physical phenotype (binary or ordinal only after evaluating outcome counts and the published scoring definition). Secondary analyses: sleep symptom measures, deficit FI and available sleep/fatigue-excluded sensitivities. Genotype remains fixed, so genotype-to-Wave 3 outcome analysis is possible conceptually, but concurrent sleep–frailty associations do not establish temporal causation.

Compare each score against the **same** baseline covariates on the **same** eligible people. Estimate association and uncertainty, and, where sample size/events support it, compare discrimination, calibration and incremental performance, using an appropriate internal bootstrap for uncertainty and no refitting of the candidate weights to the test outcomes. Fix a plan for missing variants and data before model fitting; a partial score is a different score. Avoid treating p-value differences across outcomes as a test of different effects. Restrict flexible models and subgroup claims given limited expected events. All results are research-only.

## Questions requiring a custodian response

1. Can the relevant genetic release and Wave 3 questionnaire, FI and Fried-component data be linked for this specific project, and what written approvals are required under the current agreement?
2. Is the requested data package free for our institution and use, including processing/transfer costs? What eligibility, IRB and secure-computing conditions apply?
3. How many eligible individuals and disclosure-safe frailty events have jointly usable sleep, physical-frailty and genotypes? Are PCs and relatedness data available?
4. Which of the requested 33 variants have acceptable call rate, exact alleles and sequencing or imputation quality in the release? Are full genome files necessary or is approved region/variant extraction supported?
5. Are the published FI and physical-frailty components and item-level codes available at Wave 3, including separate sleep and fatigue items? Can the participant or relative overlap with discovery/selection cohorts be checked within approved governance?

This outline is for review and access feasibility only; it has not been sent, signed, approved, or used to analyze restricted data.

The [published FI item audit](lbc1936_published_fi_item_audit.md) finds no explicit sleep/fatigue term among 30 top-level labels, notes shared physical-function items, and identifies source-table coding ambiguities to resolve before reconstruction.
