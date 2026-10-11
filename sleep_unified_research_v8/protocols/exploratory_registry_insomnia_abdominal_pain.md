# Frozen exploratory source/construct transport experiment

Date: 2026-10-11. Freeze precedes opening the new pair estimate. Selection used historical results and source feasibility; this is not an outcome-blind confirmatory study. One pair, one sleep source, one outcome source, one stock LDSC specification; no source, phenotype, filter or estimator retries after inspecting the result.

## Question and estimand

Estimate marginal common-variant genetic correlation between FinnGen R13 registry-coded insomnia F5_INSOMNIA and MVP GIA European abdominal pain GCST90479148. This asks whether the previously reported complaint/symptom sharing has corresponding support under a different sleep source and ascertainment. Cohort, ascertainment, recruitment and measurement change together; no pure measurement effect is identified. New clinical-source evidence could support, fail to support or remain too imprecise for a source-transport interpretation. Preserve every outcome.

MVP abdominal pain is the same outcome sample used historically. Its previous estimate was inspected and contributed to selection. FinnGen sleep is newly paired with it; the result is not fully independent two-trait replication, an unbiased replication of a selected effect, a newly discovered relationship, or a calibrated difference from the historical estimate. A related 2025 clinical-versus-habitual-sleep preprint limits broad conceptual novelty.

## Verified existing inputs and qualified assumptions

Use the exact two existing compressed HM3 derivatives in source_review/qualified_native_input_contract.json. The independent source reviewer verified their current SHA-256 against previous receipts. Do not download or reprocess a verified source.

FinnGen: 51,643 cases/446,273 controls, assumed constant effective N=185146.703773, ALT effect, GRCh38 source; the existing SNP identity/harmonization derivative is reused. MVP: assumed constant effective N=274973.347782, existing source allele/OR/log-CI uncertainty derivative reused. MVP supplied SE is unavailable in the inspected sample; uncertainty is reconstructed under the explicit assumption of 95% approximately log-Wald CIs. GWAS-SSF fields alone do not certify confidence level or construction. Supplied SAIGE P can differ from a Wald P; the bounded discrepancy diagnostic is not a validity failure or proof of a source error. No per-variant N or INFO availability is fabricated. Source-specific HM3 fallback and constant-N assumptions remain limitations.

Both named biobank sources differ from UKB. Individual cross-enrollment or relatives between FinnGen/MVP/UKB is not certified. Use the stock estimated cross-trait intercept; do not constrain it to zero. European reference matching is approximate for the Finnish source, not independently proven optimal. Existing standalone h² diagnostics passed necessary historical conditions under their assumptions; pair-specific diagnostics are reported anew. They do not establish prospective power for this rg. No arbitrary minimum detectable effect is asserted.

## Estimator and endpoints

Use the existing pinned Python/LDSC environment, unmodified stock estimator and previously reviewed return-value capture. Standard reference/weight prefix from the input contract, 200 stock block deletions, default filters and unconstrained h²/cross-trait intercepts. Primary output is rg, its stock jackknife SE, marginal 95% interval, z and P, both pair-specific h²/SE/intercepts, gencov/intercept, retained SNP count and warnings. No new estimator or covariance model is developed.

There is one exploratory result. Report its nominal uncertainty as conditional on source/estimator assumptions. Show P<0.05/217 only as a conservative descriptive historical-threshold sensitivity; this is not a new complete 217-pair confirmatory family or a correction for winner's curse, source selection, repeated research or shared outcomes. No calibrated cross-source difference test is computed. Core396/extension1200/historical217 families and all their outputs remain unchanged.

If the estimator fails, yields unstable/nonfinite output, exceeds source QC gates, or triggers a resource guard, retain the failed/qualified result without changing source or method to obtain success. A favorable point sign alone is not validation. Interpretation must state interval width, source qualifications and the historical outcome reuse.

## Execution and preservation

Run one heavy worker under the existing shared native lock and unchanged reviewed process-monitor/worker functions. Keep at least3GiB internal and5GiB SSD free, worker RSS≤2GiB, new native outputs≤4GiB, deadline30minutes. Write new results only under this V8 SSD namespace. Preserve interrupted/failed outputs and inspect live handles before any continuation; no completed native fit is rerun. No new controller version, source download, historical overwrite, main merge or restricted request.

Independent numerical/scientific review of the actual new output is required before retaining any scientific interpretation. This protocol admits this qualified exploratory fit only; it does not clear public release, clinical diagnosis equivalence, causality, mechanism, independent replication or journal acceptance.
