# Integrated frailty preflight after source-abstract audit

- Completed: 2026-09-24 05:57 UTC.
- Command: `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' make -C frailty_paper PYTHON=.venv/bin/python preflight`
- Starting HEAD: `ef8be115e59ebe4e8840fc887c6c289050cd5ba6` (`Audit missing abstracts against source XML`).
- Worktree: mixed and dirty, including unrelated Brain6/root changes and pre-existing frailty resource-manifest/data-path changes. The preflight does not validate those unrelated files.
- Record type: summary of completed command output; not a byte-for-byte stdout capture.

## Results

| Check | Result |
|---|---|
| Acquired-resource verifier | PASS — 351 manifest rows; all 351 files verified by size and SHA-256 against the mounted external root |
| Acquisition-manifest metadata | PASS — 351 rows; all 22 required columns |
| Frozen analysis plan | PASS — v1; 12 sleep traits; multiplicity 396 |
| Cohort-overlap ledger | PASS — 34 rows; exact participant intersections remain unknown |
| Review-screening structure | PASS — 56,092 title/abstract records; zero decisions and zero full-text records |
| Bibliography audit | PASS — 29 citation uses, 27 entries, no errors |
| Missing-abstract source audit | PASS — all 935 PMIDs found in the checksum-verified 311 retained PubMed XML files; 53 duplicate PMIDs; no duplicate occurrence contains abstract text |
| Latent-factor Q_SNP audit | PASS — 32,035,589 source rows; zero significant Q variants or excluded rows |
| Reporting checklist locators | PASS — 92 rows, 93 line references, no errors |
| Frailty package suite | PASS — 97/97 tests |

This verifies package integrity and software gates at the stated inputs. It does not establish study eligibility or resolve manual licensed-database searches, dual screening/adjudication, exact cohort overlap, source-specific GWAS blockers, independent replication, or downstream locus/mechanism analyses.
