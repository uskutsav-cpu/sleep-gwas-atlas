# Live-HEAD frailty preflight — 2026-09-24

- Branch: `frailty-paper-v1`
- Source HEAD at launch: `16fa568a69e3364dd22d75ad35e3ca9175958fc1`
- Command: `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' make -C frailty_paper PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python preflight`
- Captured stdout/stderr: `preflight_live_head_py311_2026-09-24_1403.log` (SHA-256 `a8b2e53318e8472f1c798dc60a3a4504c65e35dd7ea90d038ab39b953889ab29`)

## Results

- Acquisition manifest: 351/351 files passed path/size/SHA-256 verification; metadata audit passed 351 rows and 22 columns.
- Locked plan: v1, 12 sleep traits, multiplicity denominator 396.
- Cohort overlap: 34 ledger rows valid; exact participant intersections remain unknown.
- Review structure: 56,092 title/abstract records, zero decisions, zero full-text records.
- Bibliography: PASS, 29 citation uses and 27 entries.
- Missing-abstract source audit: PASS, all 935 missing-abstract PMIDs found in source XML; 53 duplicate PMIDs noted; no abstracts recovered.
- Manuscript quantitative claims: 28/28 PASS.
- Latent-factor Q_SNP scan: 32,035,589 source rows, zero significant variants and zero excluded rows under the locked threshold.
- Reporting locators: 92 checklist rows and 93 references, zero errors.
- Pinned Python 3.11.11 test suite: 112/112 PASS.

The run completed successfully on the live mixed worktree. It is not a clean-checkout attestation and does not resolve scientific eligibility, exact sample overlap, manual database exports/screening, independent replication, or downstream local/molecular gates. Four generated audit JSONs were refreshed against the current inputs during this preflight.
