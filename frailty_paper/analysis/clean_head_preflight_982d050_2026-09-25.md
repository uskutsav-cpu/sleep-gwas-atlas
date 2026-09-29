# Clean-archive frailty preflight — source HEAD `982d050`

## Result

`make -C frailty_paper preflight` completed on a clean Git archive of source commit `982d050ca900c56a32a6ee1d2a8afac42c1ccf7e`, using Python 3.11.11. All ten package audits and the 130-test suite passed. The full captured output is `preflight_clean_head_982d050_2026-09-25.log` (SHA-256 `ff8b0e16a0c48915d33975688a4ce05cea8f6548f227b470e5ca97dccd257d4f`).

The archive validation used the configured external data root. It verified the 311 PubMed XML source files by recorded byte count and SHA-256 before materializing a temporary copy inside the clean archive, because the missing-abstract source audit requires XML beneath the review tree. The temporary clean archive and XML copies were removed after the run; no external inputs or the mixed checkout were changed by the preflight.

## Gate results

- Acquired-resource file integrity: 354/354 rows verified by path, bytes and SHA-256; metadata: 354 rows and 22 required columns.
- Analysis plan: locked v1 accepted; 12 sleep traits and 396-pair multiplicity family.
- Cohort-overlap ledger: 34 rows valid; exact participant intersections remain unknown.
- Review structure: 56,117 title/abstract records, zero decisions, zero full-text records.
- Bibliography: PASS; 29 citation uses, 27 unique cited keys, no missing, duplicate or uncited keys.
- Missing-abstract source audit: 935/935 PMIDs found in source XML; none had an abstract in the checked XML; 53 duplicate source records.
- Manuscript quantitative claims: 28/28 passed the current claim audit.
- Latent-factor Q_SNP scan: 32,035,589 source rows; seven factors; zero variants at the locked threshold and zero excluded rows. This validates the specified exclusion only, not equivalence to the source paper's full post-GWAS pipeline.
- Reporting locators: 92 checklist rows, 93 line references, zero range errors.
- Frailty unit suite: 130/130 passed.

## Scope and current-head relation

At the time of this post-validation review, branch HEAD was `2db8efcf928a105f12b6123f966b8ebb8b202f0d`. This historical integrated preflight predates the later narrow Q_SNP audit-script correction (`<` to the author pipeline's inclusive `<=`) and its boundary regression test. After that correction, the full seven-factor Q_SNP audit and the pinned Python 3.11.11 suite passed separately (32,035,589 rows scanned; zero hits/exclusions; 131 tests passed with the untracked Brain6-only test module excluded). The complete integrated `preflight` has not been rerun after the code/test correction, so this report establishes the 982d050 package snapshot only; it is not a full current-HEAD preflight. Subsequent LAVA and source-fidelity reports record the separate targeted validations and current pause/resource state.

This package preflight is not a full scientific end-to-end replay. It does not resolve the 56,117 unscreened review records, missing licensed database exports, unresolved participant overlap or replication, failed LAVA completeness/QC gates, unavailable alternate frailty phenotypes, or downstream locus/molecular analyses. Overall manuscript readiness remains **NO-GO**.
