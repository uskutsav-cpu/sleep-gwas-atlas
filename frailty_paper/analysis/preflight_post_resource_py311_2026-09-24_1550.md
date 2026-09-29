# Integrated frailty preflight after resource registration

- Run window: 2026-09-24 15:50–15:56 UTC (captured log mtime 15:56:32 UTC).
- Source HEAD: `f111f9fc334474a73a86ad9cc8fe7db5828c16cc`.
- Worktree: mixed; unrelated Brain6/root changes were outside the audited package scope.
- Interpreter: pinned isolated Python 3.11.11 at `work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python`.
- Command: `make -C frailty_paper preflight PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'`.
- Exit status: 0.
- Results: 351/351 acquired-resource files verified; all 351 manifest rows passed 22-column metadata validation; analysis-plan lock version 1, 12 sleep traits and 396-pair family passed; the 34-row cohort ledger passed with exact intersections still unknown; 56,092 review records remain undecided; bibliography audit passed; 935 missing-abstract PMIDs were found in source XML and none had abstract text; 28/28 manuscript claims passed; Q_SNP audit scanned 32,035,589 source rows with zero significant variants and zero exclusions under the source-defined threshold; 92 checklist rows/93 line references passed; all 112 package tests passed.
- Captured log: `analysis/preflight_post_resource_py311_2026-09-24_1550.log` (SHA-256 `94e3de74a62ec9d55527acb8f5f01f9041c4c45f47375cb8fea33ef826c8ab1a`).
- Limits: validates the frailty package against mounted external inputs, not a clean checkout or unrelated Brain6/root workflows; does not resolve screening, access, overlap or scientific eligibility blockers.
