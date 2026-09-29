# Live-head pinned preflight — 2026-09-24 14:51 UTC

- Branch: `frailty-paper-v1`
- Source HEAD at launch: `ec682d2b69fffbea56bf2b7e88cb35135cc846cf`
- Command: `FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' make -C frailty_paper PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python preflight`
- Captured stdout/stderr: [`preflight_live_head_py311_2026-09-24_1451.log`](preflight_live_head_py311_2026-09-24_1451.log), SHA-256 `e6ba9982b263090dcabb97a015139086f87fcb8636384f0946af63e57af739a3`.

## Results

- Acquired-resource integrity: 351/351 files verified; metadata audit: 351 rows and 22 required columns.
- Frozen plan: pass (v1, 12 sleep traits, multiplicity denominator 396).
- Cohort-overlap ledger: pass (34 rows); exact participant intersections remain unknown.
- Review structure: pass (56,092 records; zero decisions; no full-text records).
- Bibliography audit: pass (29 citation uses, 27 entries).
- Missing-abstract source audit: pass (935 missing-abstract PMIDs confirmed absent from all retained XML occurrences; 53 duplicate PMIDs).
- Manuscript quantitative-claim audit: 28/28 pass.
- Latent-factor Q_SNP audit: 32,035,589 rows, zero significant variants, zero excluded rows; seven-factor audit passes.
- Reporting locators: 92 rows, 93 line references, zero errors.
- Pinned Python 3.11.11 package suite: 112/112 pass.

The checks ran on the live mixed worktree at the recorded HEAD, with mounted external data. The complete output log is retained. The preflight establishes these code/input-integrity and package checks only; it does not resolve screening, licensed database searches, independent replication, exact participant overlap, alternate frailty-source eligibility, or the failed LAVA completeness gate.
