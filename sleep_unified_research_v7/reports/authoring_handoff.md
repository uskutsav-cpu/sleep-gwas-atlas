# Authoring handoff: one integrated sleep GWAS resource paper

## Decision

**Human manuscript authoring:** ready to begin one carefully qualified, integrated resource manuscript using the frozen V4 evidence and the final V7 benchmark. **Nature Communications submission readiness:** not supported by the current evidence. The package establishes no new method, independently replicated two-trait result, calibrated contrast, or major biological advance. No manuscript prose, abstract, introduction, results narrative, discussion, cover letter, or submission material was written.

Keep the 396-row core, 1,200-row extension, and historical 217-candidate validation family distinct. The 41 estimated external rows reuse discovery sleep sources and are outcome-side validation, not independent two-trait replication. Treat V7 as a bounded resource-coverage analysis; it supports neither global novelty nor a human usability benefit.

## Evidence and manuscript-ready assets

- V4 preserves 396 core correlations, 1,200 extension correlations, 158 heritability estimates, 62 sensitivity estimates, and 41 qualified external outcome-side estimates. Statistical results were independently reproduced at the recorded precision. These results predate V7; no V7 GWAS or LDSC analyses ran.
- The interrupted validation source-12 download resumed from its preserved 57,671,680-byte checkpoint. All 13 validation source bodies and all 13 independently reviewed preprocessing/QC outputs are complete across 14,377,180 retained rows. Verified sources were reused rather than downloaded again.
- The machine-readable candidate is built from frozen inputs and contains the 396/1,200/41 estimate tables, 158 h² rows, 62 sensitivity rows, source/phenotype metadata, and schema. Its independent automated consumer check and isolated reproducible build passed. It remains a **candidate for rights and metadata review**, not a cleared public release.
- V4 already has a 23-figure evidence set and a five-panel main-figure selection. Reuse these existing assets; do not rerender unchanged plots. V7 adds `../figures/cross_resource_coverage.svg`, backed by the annotated TSVs and manifest hashes.
- Final V7 benchmark: 6,548 unique annotations (1,637 result records × four accessible snapshots); A/B/C/D/E union = 0/1,183/0/366/88; all 1,247 apparent candidates are B/C/E under the strict source/estimand rule. The 361 qualified match-absent rows (85 core, 276 extension) are D in each named snapshot. Full candidate references are retained.
- Independent V7 statistical and scientific review findings are saved in `../reviews/independent_statistical_review.md` and `../reviews/independent_scientific_review.md`. These automated reviews challenge actual result and source rows; they are not human peer review, empirical replication, or a human user study.

## Seven requested priorities

| Priority | Status and evidence | Boundary before public release or submission |
|---|---|---|
| Benchmark existing resources | Complete for the V4 literature review and V7 comparison of HGA release 3, Morrison S10–S15, Goodman SD24, and Fan SD2/SD4/SD11. | SleepChart remains prior-audit context; CTG-VL has no captured direct export. The 361 count is bounded to the four searched snapshots. |
| Source-aware harmonization | Complete as a research layer: 12-construct dictionary, outcome-source definitions, 83-field crosswalk, source IDs/releases, coding notes, explicit unknowns, and complete V7 candidate-row references. | Investigators must approve unresolved source definitions, participant overlap, and source-specific transformations. |
| Uncertainty and reliability | Complete with reported marginal intervals, fixed QC flags, and an explicit reliability/covariance limitation framework. | No calibrated discovery/validation covariance, phenotype reliability, or between-measurement contrast is available. |
| Biological use cases | Complete as four bounded retrieval examples using existing V4 estimates and limitations. | These examples are not V7 discoveries or causal/mechanistic results. |
| Machine-readable release | Complete as a reproducible candidate package with schema and SHA-256 inventory. | Human investigators must clear source-specific rights, attribution, metadata, repository, and license before deposit. Fan-derived row references require explicit rights review. |
| Independent quality/usability validation | Complete for independent automated candidate-consumer review plus V7 statistical and scientific challenges of actual tables/source pairs. | This is not a human usability study, peer review, or independent cohort replication. |
| Reproducibility and documentation | Complete for the frozen result tables, machine-readable candidate, and V7 annotation build; isolated candidate rebuild passed, and V7 output table hashes verify. | Full raw-to-processed provenance remains partial for optional core/extension replay tails. Validation replay is complete. |

## Existing results that must remain qualified

The 41 outcome-side estimates reuse sleep GWAS sources from discovery. The original 23 positives therefore do not establish independent two-trait replication. All intervals are marginal; the discovery/validation covariance is not calibrated. The 62 sensitivity rows do not test between-fit differences. Single-item questionnaire phenotypes do not establish clinical diagnoses or measured reliability. The 30 HGA bibliographic candidates were downgraded to E because publication metadata did not establish the exact source/estimand; 16/30 showed opposite-signed rg values in the independent source check. This comparison is a source-classification warning, not a new biological finding.

## Exact remaining critical path

1. **Investigator rights decision:** determine whether the candidate fields and derived cross-resource annotations can be shared, with particular review of Fan's CC BY-NC-ND materials, and settle attribution/repository/license terms.
2. **Scientific metadata approval:** confirm phenotype definitions, source versions, participant overlap, ancestry and effect coding; retain explicit unknowns where unavailable.
3. **Institution and author approval:** complete the applicable ethics/data-use determination and approve authorship, contributions, funding, conflicts, AI-use disclosure, and final interpretation.
4. **Submission planning:** decide journal fit after authors review the evidence. Current Nature Communications guidance expects important specialist advances and transparent data availability; it requires code availability when custom code is central. Data/code permissions and the current evidence do not satisfy a submission-readiness conclusion.
5. **Human final review:** investigators should verify all key table/figure interpretations and decide whether the bounded resource contribution merits one integrated paper. No scientific calculation remains on the V7 path.

## Optional work, not a blocker to drafting

- The preregistered 60-task paired AI-answer benchmark was deferred before task sampling. Its hypothesis is untested; running it would not establish human researcher usability or create biological evidence.
- A complete SleepChart row-level rejoin would require its previously inspected workbook, which is absent from the active checkout and was not re-downloaded. CTG-VL has no direct export in this evidence set.
- Source-to-processed comparison tails remain partial: extension 79/100; separate core replay attempts 5/27 and 7/35. These do not change the completed estimates and are optional provenance improvements for the current manuscript scope. No new controller or raw-source replay is admitted without a concrete defect and a defined scientific consequence.
- Additional covariance grids, local-architecture searches, broad association scans, and repeated native calculations are not justified by this evidence and are not on the critical path.

## Running work and new evidence

The final V7 build, figure, and post-fix reviewer checks are complete. No research job is currently running. The validation checkpoint recovery and preprocessing were already complete before this V7 build. V7 generated new descriptive source-coverage annotations and exposed that metadata-only HGA matches cannot be called exact; it generated **no new association, calibrated contrast, independent replication, causal claim, mechanism, or biological discovery**.
