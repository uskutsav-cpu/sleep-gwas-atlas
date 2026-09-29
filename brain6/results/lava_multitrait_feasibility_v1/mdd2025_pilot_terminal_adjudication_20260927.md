# MDD2025 LAVA pilot terminal adjudication — 2026-09-27

The reconnected Extreme SSD was mounted at the pinned path. The original
four-worker launch had empty logs, no result shards, and no exit receipt; its
recorded PIDs were now unrelated system processes. All 111 materialized inputs,
worker partitions, frozen configuration/rule hashes, and the pinned Rscript
hash validated before any resume. A first isolated resume attempt exited 255
for every worker without producing a shard: the copied R runtime's temporary
Conda prefix had disappeared. Restoring that prefix as a symlink to the
unchanged SSD environment made R 4.3.3 and LAVA 0.1.5 load. The second
isolated attempt preserved both earlier sets of logs and receipts, used the
same four worker configs, and exited 0 for all four workers.

| Worker | Completed / expected loci | TESTED | NOT_RUN | FAILED |
|---|---:|---:|---:|---:|
| 1 | 28 / 28 | 28 | 0 | 0 |
| 2 | 28 / 28 | 28 | 0 | 0 |
| 3 | 28 / 28 | 27 | 1 | 0 |
| 4 | 27 / 27 | 25 | 2 | 0 |

The [pilot decision](</Volumes/Extreme SSD/brain6-work/lava-multitrait-feasibility-v1/mdd2025_trait_pilot_v1/results/pilot_decision.json>)
is a **numerical trait-only diagnostic pass** under the frozen rule: 108/111
pilot loci were testable versus 83/111 in the same canonical MDD loci, a
25-locus or 22.52 percentage-point gain. Low-local-h² `NOT_RUN` fell from 27
to 3. The audit independently re-read the terminal shards, exact locus
partitions, summaries, `SCREEN_COMPLETE` logs, zero exit codes, launch chain,
and SHA-256 hashes; a second audit verified the existing aggregate and
decision without rewriting them. The machine-readable
`mdd2025_pilot_terminal_adjudication_20260927.json` binds the key receipts.

**Scientific interpretation remains on hold.** The
[pre-outcome sample-size audit](mdd2025_pilot_sample_size_audit_20260927.md)
identified that the pilot passed source `NEFF` as LAVA's binary `N` while
providing one aggregate case fraction. The PGC source header defines `NEFF`
as the sum of cohort-level effective sizes, whereas its `NCAS`/`NCON` are
variant-level counts. A read-only scan of all 7,363,302 source rows found
`NEFF / (4*NCAS*NCON/(NCAS+NCON))` ranging from 0.8593 to 0.9529 (mean
0.8830). This descriptive scan did not use effect, SE, or P-value fields and
does not select a replacement N model after seeing the pilot. The source file
still matches SHA-256
`5d6fc5aee638e73457da703b75e0b4d4dacfb1b87578fd464b771209c0dfb22a`.

**Family gate: FAIL_QC.** The pilot's numerical `advance_to_full_trait_only_screen`
flag cannot override the pre-outcome interpretation hold. Independently, even
a perfect MDD replacement would leave 3,131 `NOT_RUN` cells among the other
six traits, above the frozen whole-family ceiling of 873. No full MDD screen,
full LAVA rescue family, candidate-region promotion, or gated downstream
analysis is authorized by this pilot. The separately versioned
[LAVA rescue-v1 QC](../lava_rescue_v1/lava_rescue_v1_qc.json) remains
`FAIL_QC`. This pilot workflow did not write to canonical v3, v2, or
roundoff outputs.
