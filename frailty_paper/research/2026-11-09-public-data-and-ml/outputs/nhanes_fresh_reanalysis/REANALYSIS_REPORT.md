# NHANES sleep–frailty reanalysis

Rerun on 1 November 2026 from public CDC NHANES 2011–2014 inputs already preserved in the study package. This is a retrospective observational analysis, not a genetic-score validation.

## Reproduction and quality

The complete reconstruction and analysis pipeline was rerun in `study/`. Its 249 independent source-item and score checks reported zero mismatches; the pipeline verification passed. The newly generated participant table (3,632 rows), score table (3,632), survey-results table (109), and supplement model table (40) exactly matched the previously audited versions after CSV parsing. The pipeline made no changes to the original study outputs. The CDC source directory in this working copy is a symbolic link, so the working directory alone is not a portable reconstruction package. Use the separately supplied original full ZIP, which includes the source files, to reproduce on another computer.

This successful pipeline rerun is a computational consistency check, not independent proof that every clinical coding choice is correct. The earlier independent audit separately verified source-file integrity against archived manifests, checked source exposure and design values, reconstructed score denominators, and reran the survey inference; see `../nhanes_manuscript_audit/audit_and_readiness.md`. The original package's source-hash verifier regenerates its manifest before verifying it, so do not interpret its own hash check as an independent authentication of the originals. An earlier verification checked the archived manifest against the original ZIP.

## Primary observational finding

Of 3,632 interviewed adults aged 60 or older, 3,472 had positive MEC weights, and 3,336 had sleep-disclosure data and met the ≥80% observed-item rule for each of four FI versions. Sleep exposure was whether a participant **ever told a health professional about trouble sleeping**, not diagnosed insomnia. Two NHANES cycles were combined using MEC weights divided by two; the full survey design was retained before domain subsetting. A linear survey model adjusted for age, sex and cycle gave:

| Adapted score, same participants | Sleep-disclosure coefficient | 95% CI |
|---|---:|---:|
| 49 items, including sleep and fatigue | 0.0736 | 0.0566 to 0.0907 |
| 47 items, excluding direct sleep and fatigue | 0.0698 | 0.0528 to 0.0869 |

For the 47-item analysis p=1.54×10⁻⁹ using 32 design degrees of freedom; the default residual-df sensitivity in the independent audit gives p=4.42×10⁻⁹ with the same coefficient and standard error. The paired coefficient difference was 0.00382 (95% CI 0.00290–0.00474); item sets and denominators both change, so this does not isolate a causal effect of overlap. The FI and score comparison were constructed retrospectively; do not characterize this analysis as prospectively preregistered.

In the prespecified-as-requested but retrospectively executed age-65+ sensitivity, n=2,326, the coefficient was 0.05954 (95% CI 0.04287–0.07621; unadjusted p=2.87×10⁻⁸, exploratory-family Holm p=3.73×10⁻⁷). Both survey-cycle estimates were positive; the cycle interaction was inconclusive. The survey-weighted correlation between the 49- and 47-item scores on the common 3,336 participants was 0.9976; this high correlation reflects largely shared items and does not establish equivalent constructs or prediction. The full set of models, correlations, completeness rules, figures and exploratory analyses is in `study/`.

The index is a NHANES adaptation of a published 49-item index, not the Atkins UK Biobank index. The adapted physical score uses proxy measures and does not reconstruct the original Fried phenotype. Remaining FI components can still be indirectly associated with sleep and health-care use. Cross-sectional associations, the disclosure-based exposure, nonresponse and confounding preclude conclusions about direction or causality.

## Genetics and access

No linked genotypes were analyzed; these results do not estimate genetic correlation, provide a sleep-excluded frailty GWAS, or validate either research genetic score. Standardized 15- and 22-variant source files are provided; FI and insomnia conventions have documented upstream support, while the FFS A1 convention is supported but its complete transformation is not author-confirmed. Multiallelic contrasts require source and target verification; see the independent audit and `../stage6_score_allele_input_gate.md`. An independent genotype–phenotype cohort and appropriate approvals remain necessary for genetic prediction testing.

The additional cohorts have not been analyzed at participant level: ACT, BLSA and HRS require requests or registration and project-specific permissions, and no adequate unrestricted linked participant data have been established in this workspace. See `../independent_outcome_validation_cohort_feasibility.md` and the ACT/BLSA-specific feasibility files. The user-requested public observational NHANES analysis is complete within its stated scope; its findings do not resolve the separate genetic questions.

## Reproduction environment

For the rerun, the project virtual environment at `.../investigate-whether-we-can-directly-reconstruct-2/work/venv` and R `survey` package in that project's `work/finalization/Rlib` were used. See `study/README.md` for dependencies and commands, `reproduction.log` for execution output, and `study/source_checksums.csv` for checksums. The final script returned exit code 0.
