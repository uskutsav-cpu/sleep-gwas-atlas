# Current-head frailty preflight — 2026-09-24 15:25–15:27 UTC

- Source HEAD before the preflight: `b12a53aaa3335ae2b75ea163a57845f2fd06ee3d` on `frailty-paper-v1`.
- Command: `make -C frailty_paper PYTHON=/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03/frailty_paper/.venv/bin/python FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' preflight`
- Interpreter: project environment Python 3.13.13.
- Exit status: 0.
- Complete captured output: `analysis/preflight_current_head_py313_2026-09-24_1525.log` (25,356 bytes; SHA-256 `98b0f7ae70e2e012070ddb49a5b8193f8c026ce119b6dc5b90f03c545887610b`).

## Results

- Acquired-resource manifest: 351/351 files verified by path, size, and SHA-256; all 351 rows passed the 22-column metadata audit.
- Locked analysis plan: version 1; 12 sleep traits; 396-pair multiplicity.
- Cohort ledger: 34 rows valid; exact intersections remain unknown.
- Review queue: 56,092 title/abstract records, zero screening decisions; validator passed.
- Bibliography: 29 citation uses, 27 entries, no missing or uncited keys.
- Missing-abstract provenance: 935/935 PubMed IDs confirmed in source XML with no abstract; 53 duplicate source PMIDs.
- Manuscript quantitative claims: 28/28 passed.
- Latent-factor Q_SNP audit: 32,035,589 source rows, zero significant Q_SNP variants and zero rows excluded under the predeclared threshold.
- Reporting locators: 92 checklist rows, 93 line references, zero range errors.
- Package suite: 112 tests passed.

This validates the frailty package checks in the live mixed worktree. It does not validate unrelated Brain6/root workflow changes, complete manual review screening, phenotype eligibility blockers, or a clean checkout. Preflight-generated audit provenance now records Python 3.13.13 for the Q_SNP run.
