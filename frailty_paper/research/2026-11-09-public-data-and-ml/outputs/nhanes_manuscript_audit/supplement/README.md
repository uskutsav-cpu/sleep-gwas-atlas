# Supplementary materials: NHANES audit

These files accompany `../manuscript_with_nhanes.md` and `../audit_and_readiness.md`. They are a review copy of selected outputs from the user-provided reconstruction ZIP and the independent source comparison, not a replacement for the full package.

* `supplement_49_item_audit.csv` gives item definitions, missingness and the direct sleep/fatigue exclusion flag; `supplement_item_audit_checks.csv` contains the original package's independent-recode comparisons.
* `supplement_survey_models.csv`, `supplement_flow.csv` and `supplement_missingness_by_exposure_cycle.csv` preserve the reproduced survey estimates, eligibility counts and item missingness. The models use the explicit design-df convention; the independently computed residual-df sensitivity is reported in the manuscript and audit.
* `figure1_fi_associations.pdf`, `figure2_distribution_mortality.pdf`, their source-data CSVs and `figure_legends.md` were reproduced from the user package. The mortality panel is a same-cohort construct check with varying follow-up, not external validation.
* `score_15_variants_standardized.csv` and `score_22_variants_standardized.csv` are the supplied experimental score tables. `raw_score_source_audit.csv` and `raw_score_source_receipt.json` report our separate numerical comparison with upstream GWAS releases. The FFS effect convention and multiallelic contrasts remain subject to the source limitations documented in the audit; the equally weighted 22-marker count is a distinct score.
* `independent_check.json` reports the independent denominator and raw exposure/design-value checks. It does not itself recode all 49 original fields.

The full package, source CDC files and scripts are in the separately located `nhanes_manuscript_sensitivity_final.zip`; the isolated rerun and our scripts are in the current workspace at `work/nhanes_independent_audit/`. Participant data are not duplicated in this review supplement.
