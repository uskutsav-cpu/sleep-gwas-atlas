# Package validation snapshot — 2026-09-23

This is a bounded validation record for the current repository state. It does not certify manuscript completeness, data access, independent replication, or downstream mechanistic analyses.

## Commands and results

| Check | Result | Scope / limit |
|---|---|---|
| `make -C frailty_paper test` | PASS — 62 tests | Existing unittest suite only; does not establish external scientific validity. |
| `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' make -C frailty_paper verify-manifest` | PASS — 347 rows, 347 files verified | Verifies registered file sizes and SHA-256 at the configured external root. |
| `make -C frailty_paper validate-plan` | PASS — v1, 12 sleep traits, multiplicity 396 | Confirms the frozen analysis-plan lock and recorded input hashes. |
| `make -C frailty_paper validate-overlap` | PASS — 34 rows | Source-level ledger agrees with the frozen 12-trait panel; exact participant intersections remain unknown. |
| `make -C frailty_paper validate-review` | PASS — 56,092 records, 0 decided, 0 full-text decisions | Structural validation only. This confirms the review remains unscreened; it is not a completed systematic review. |
| `git diff --check` | PASS | Whitespace validation for the current worktree at audit time. |

## Interpretation

The package validators and unit/integration tests pass in the current local environment. The scientific and reporting gates remain as listed in `../STATUS.md` and `../BLOCKERS.md`; in particular, cohort overlap is not known exactly, manual database exports and screening are absent, and the manuscript package remains provisional.
