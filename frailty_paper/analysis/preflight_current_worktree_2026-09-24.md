# Current-worktree integrated preflight

- Completed: 2026-09-24 UTC; command completed by 04:13 UTC.
- Command: `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' make -C frailty_paper PYTHON=.venv/bin/python preflight`
- Branch: `frailty-paper-v1`. The last observed source HEAD when starting was
  `137fc9f`; while the command was running, HEAD advanced to `c6f5527` through
  a Brain6-only commit. `git diff 137fc9f..c6f5527 -- frailty_paper` is empty.
- Worktree: mixed and dirty. The current acquisition manifest and unrelated
  root workflow/environment files were part of the live worktree; no claim is
  made that the unrelated modified workflow/environment paths were validated
  by this package preflight.
- Record type: condensed summary of the successful command output, not a
  byte-for-byte stdout capture.

## Results

| Check | Result |
|---|---|
| Acquired-resource verifier | PASS — 351 rows; 351 files verified by size and SHA-256 against the external root |
| Manifest metadata | PASS — 351 rows; all 22 required columns |
| Frozen analysis plan | PASS — v1; 12 sleep traits; multiplicity 396 |
| Cohort-overlap ledger | PASS — 34 rows; exact participant intersections remain unknown |
| Review-screening structure | PASS — 56,092 title/abstract records; zero decisions and zero full-text records |
| Bibliography | PASS — 29 citation uses; 27 entries; no errors |
| Latent-factor source Q_SNP audit | PASS — 32,035,589 rows; threshold `5e-8/7`; zero significant variants and zero exclusions |
| Reporting locators | PASS — 92 checklist rows; 93 line references; zero errors |
| Frailty package tests | PASS — 94/94 |

This establishes acquisition integrity and workflow-validator/test status for
the live package inputs at this run. It does not resolve sample overlap,
independent replication, manual screening, source eligibility, or any
downstream locus/mechanism gate.
