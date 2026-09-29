# Post-HFRS integrated preflight record

- Run date: 2026-09-24 UTC (completed after commit `1822c8a5996b252d5ed40a2f7831ea83071663fe`).
- Command: `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' make -C frailty_paper PYTHON=.venv/bin/python preflight`
- Worktree scope: mixed worktree; the preflight read the committed 351-row acquisition manifest and external resource files. No Brain6 results were inputs to this check.
- Record type: condensed summary of the terminal output from the completed command, not a byte-for-byte stdout capture.

## Results

| Check | Result |
|---|---|
| Acquired-resource verifier | PASS — 351 rows; 351 files verified by size and SHA-256 |
| Manifest metadata | PASS — 351 rows; 22 required columns |
| Frozen analysis plan | PASS — v1; 12 sleep traits; multiplicity 396 |
| Cohort overlap ledger | PASS — 34 rows; exact intersections remain unknown |
| Review screening structure | PASS — 56,092 title/abstract records; zero decisions; zero full-text records |
| Bibliography | PASS — 29 citation uses; 27 entries; no errors |
| Latent-factor source Q_SNP audit | PASS — 32,035,589 rows; threshold `5e-8/7`; zero significant variants; zero excluded rows |
| Reporting locators | PASS — 92 checklist rows; 93 line references; zero errors |
| Frailty package tests | PASS — 94/94 |

This confirms package, manifest and recorded-input checks only. It does not resolve sample overlap, independent replication, manual screening, source eligibility or any downstream locus/mechanism gates.
