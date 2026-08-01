# PROJECT STATE

**Branch** `claude/neuro-immune-cardio` · **HEAD** `6b16671` · pushed to origin

## Machine-derived counts (regenerate: `scripts/23_phase1_reconcile.py`)

| | |
|---|---|
| traits in ledger | 149 |
| Phase 1 complete (munged + parsable h2) | **19** |
| h2 PASS / DROP | 19 / 0 |
| **core sleep anchors complete** | **6 / 6** |
| core_sleep | 6 / 6 |
| disease | 13 / 121 |
| secondary_sleep | 0 / 14 |
| aging | 0 / 8 |
| **frozen primary rg family** | **78 pairs** (6 core sleep x 13 disease) |
| surviving BH-FDR < 0.05 | **32** |

Counts reconcile across filesystem, ledger and summary table — enforced by
`scripts/23_phase1_reconcile.py`, which exits non-zero on drift.

## Done this run

- **insomnia recovered.** The prior "unavailable" verdict was wrong; it tested
  stale CTG URLs. Live download found via the official CNCR page. **All six
  core sleep anchors are now complete.**
- **Count discrepancy resolved.** It was two traits missing from a stale
  h2_summary (`ckd`, `heart_failure`) plus a truncated `atrial_fibrillation`
  h2 log being miscounted as a DROP.
- **Primary FDR family FROZEN** at 78 pairs, table sha256 recorded in
  `results/phase1/fdr_families.tsv`. Earlier 15- and 40-pair matrices were
  interim and are superseded.
- **Phase 2/3 toolchain installed and verified**: R 4.6.1, LAVA 0.1.5, PLACO+
  vendored at a pinned commit. `scripts/21_verify_phase23_tools.sh` →
  TOOLCHAIN_VERIFIED.
- Scientific language corrected throughout; CLAUDE.md's false long-sleep
  expectation replaced with the measured result.

## Exact next actions

1. `bash scripts/20_install_phase23_tools.sh` is done; next download the
   **LAVA LD-block reference** (not yet present) — Phase 2 cannot run without it.
2. `brew install octave`, then run the pleioFDR chr-21 demonstration and record
   SUPPORTED_VALIDATED / RUNS_WITH_DIFFERENCES / UNSUPPORTED.
3. Continue Phase 1 breadth:
   `.ldsc-env/bin/python scripts/14_phase1_batch.py --role disease --limit 3`
   repeated — 108 disease traits remain.
4. Re-freeze the primary family only when the disease set is final.

## Hard constraints

8 GB RAM · disk fluctuates in single-digit GB · background jobs are killed when
the agent window closes, so bulk work must run as repeated foreground batches
or be launched by a human in a terminal.
