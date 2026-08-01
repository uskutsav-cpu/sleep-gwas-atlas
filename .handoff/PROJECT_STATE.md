# PROJECT STATE

**Branch** `claude/neuro-immune-cardio` · **Remote** github.com/uskutsav-cpu/sleep-gwas-atlas

## Stage

Phase 1 **partially complete**. Phase 2 and Phase 3 **not started** — blocked on
toolchain (see BLOCKERS.md).

## What is done

- **Stage A complete.** Environment audited, `results/phase_ledger.tsv` built
  (149 traits × 41 fields, rebuilt from primary evidence each run).
- **5 of 6 core sleep anchors through Phase 1**, all passing the QC gate:
  `sleepdur`, `shortsleep`, `longsleep`, `sleepiness`, `chronotype`.
  `insomnia` is blocked (see BLOCKERS.md).
- **13 traits have h²**; all 13 PASS the Z ≥ 4 and intercept ≤ 1.20 gate.
- **Primary rg matrix: 5 core sleep × 8 disease = 40 pairs, 15 survive BH-FDR.**

## Exact next actions

1. `python3 scripts/14_phase1_batch.py --role disease --limit 3` — repeat until
   the disease count stops rising. Each trait takes 4–16 minutes; the driver is
   resumable and skips anything already munged.
2. Re-run the primary rg matrix after each batch (command in RUNBOOK.md).
3. Install R before attempting Phase 2 (`brew install r`), then LAVA.
4. Resolve `insomnia` — it is the one missing core anchor and it matters more
   than any additional disease trait.

## Counts (verify with `scripts/14_phase1_batch.py --status`)

| | |
|---|---|
| traits in ledger | 149 |
| Phase 1 complete (munged + h²) | 13 |
| core sleep complete | 5 / 6 |
| disease complete | 8 / 121 |
| secondary sleep / aging complete | 0 / 22 (Codex munged 6 in its own tree) |
| documented exclusions + blockers | 24 |

## Hard constraints that shape everything

- **8 GB RAM, single-digit GB free disk.** Harmonization must use
  `--chunksize`; raw and harmonized files are evicted after munging.
- Raw inputs are reproducible from `config/public_gwas_sources.tsv`
  (url + sha256 + bytes). Verified: re-downloads return byte-identical hashes.
- LDSC lives at `ldsc/`, `.ldsc-env/`, `ref/` — **symlinks** into
  `~/Documents/Codex/2026-07-28/here/sleep-gwas-atlas/`. Do not delete that tree.
