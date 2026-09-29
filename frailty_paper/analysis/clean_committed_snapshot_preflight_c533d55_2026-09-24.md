# Clean committed frailty snapshot preflight

**Source commit:** `c533d55ca49972efa34bf7a5e6083202fdd80a3f` on `frailty-paper-v1`.
**Run date:** 2026-09-24 UTC.
**Snapshot:** Git archive at the source commit, extracted to `/private/tmp/frailty-clean-current-20260924`; no uncommitted checkout files were included.
**Python:** 3.11.11 from `work/conda-envs/frailty-py311/bin/python`.
**Command:** `FRAILTY_EXTERNAL_STORAGE_ROOT="/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1" make -C /private/tmp/frailty-clean-current-20260924/frailty_paper PYTHON=/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03/work/conda-envs/frailty-py311/bin/python preflight`.

The four runtime data roots (`data/raw`, `data/harmonized`, `data/munged`, and `frailty_paper/data`) were linked read-only to the configured external volume. All 311 PubMed XML source files were hard-linked from the prior verified snapshot only after each size and SHA-256 matched this snapshot's `pubmed_source_files.tsv`. The clean tree did not contain its generated review output queues, so the documented `scripts/08_build_review_records.py` builder was run in the isolated snapshot. It reproduced 60,989 source records, 56,092 deduplicated records and 4,897 duplicate occurrences; output-manifest SHA-256 `7301215e3d9fa40cfd479b0afbc4cc73b3cd3e8c1c4c227bda22e6055480379c` exactly matches the committed manifest value.

The complete preflight passed: 351/351 acquired-resource files; 22 required metadata columns; locked 12-trait / 396-pair plan; 34-row cohort ledger; review validation (56,092 records, zero screening decisions); bibliography; 935/935 missing-abstract mappings; 28/28 manuscript quantitative claims; 32,035,589 Q_SNP rows with zero significant variants and zero excluded rows; 92 reporting rows / 93 source locators with zero errors; and 112/112 tests.

Captured command output: `analysis/preflight_clean_committed_snapshot_c533d55_py311_2026-09-24.log`, SHA-256 `01b7f9c105b451157bbc9ad3cd527b3c826991cbf5e9ae05104cea44a75c6568`. The mixed-worktree preflight recorded separately in `preflight_external_root_project_py313_2026-09-24_2027.md` used `frailty_paper/.venv/bin/python` (Python 3.13.13), not Python 3.11.11; its prior interpreter label has been corrected.

This verifies the committed software snapshot against the mounted data and recorded contracts. It is not a data-independent end-to-end scientific replay and does not resolve access, screening, sample-overlap, replication, or failed LAVA-gate blockers.
