# Integrated frailty preflight with manuscript-claim audit

- Completed: 2026-09-24 06:26 UTC.
- Command: `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' PYTHON='../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python' make -C frailty_paper preflight`
- Python: 3.11.11 isolated environment.
- Starting HEAD: `3cd8436ef07d5a70aba4bcca2bda66bdf9da5c2e`.
- Worktree: mixed and dirty; unrelated Brain6/root paths and pre-existing frailty resource/data-path changes were outside the frailty package checks.
- Full output: `preflight_post_abstract_audit_py311_2026-09-24.log`.

## Results

| Check | Result |
|---|---|
| Acquired-resource verifier | PASS — 351/351 files by size and SHA-256 |
| Acquisition-manifest metadata | PASS — 351 rows; 22 required columns |
| Frozen analysis plan | PASS — v1; 12 sleep traits; multiplicity 396 |
| Cohort-overlap ledger | PASS — 34 rows; exact intersections remain unknown |
| Review-screening structure | PASS — 56,092 records; zero decisions; zero full-text records |
| Bibliography | PASS — 29 citation uses, 27 entries, no errors |
| Missing-abstract source audit | PASS — 935/935 PMIDs in verified retained PubMed XML; 53 duplicates; no abstract text in any occurrence |
| Manuscript quantitative claims | PASS — all 28 selected text/source comparisons |
| Latent-factor Q_SNP audit | PASS — 32,035,589 rows; zero significant Q variants or excluded rows |
| Reporting locators | PASS — 92 rows; 93 line references; zero errors |
| Frailty package suite | PASS — 97/97 tests |

This is an engineering, source-integrity and selected-claim consistency snapshot. It does not establish eligibility, complete manuscript accuracy, participant independence, completed screening, or the downstream analyses gated on eligible and independent inputs.
