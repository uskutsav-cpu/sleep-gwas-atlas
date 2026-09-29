# Brain6 terminal checkpoint — 2026-09-27

## Verified execution

The Extreme SSD is mounted at `/Volumes/Extreme SSD`. The original MDD2025
trait-only pilot launch had no terminal shard or exit receipt. Its four recorded
PIDs no longer belonged to the pilot. The resume preserved the original launch,
logs, configurations, materialized 111-locus inputs, and output directory. The
first isolated attempt exited 255 before work because the pinned R runtime's
temporary Conda prefix was absent; the prefix was restored as a symlink to the
same SSD environment, without changing Rscript or input hashes. The second
attempt completed under the same four pinned worker partitions:

| Worker | Completed / assigned | TESTED | NOT_RUN | FAILED | Exit |
|---|---:|---:|---:|---:|---:|
| 1 | 28 / 28 | 28 | 0 | 0 | 0 |
| 2 | 28 / 28 | 28 | 0 | 0 | 0 |
| 3 | 28 / 28 | 27 | 1 | 0 | 0 |
| 4 | 27 / 27 | 25 | 2 | 0 | 0 |

All 111/111 pilot loci are complete; there are no unfinished shards. The
[terminal adjudication](../lava_multitrait_feasibility_v1/mdd2025_pilot_terminal_adjudication_20260927.md)
binds the source, exit, rule, and decision hashes. The receipt-bound pilot
decision is `PASS_TRAIT_ONLY_DIAGNOSTIC`: 108/111 MDD loci tested versus 83/111
at those same canonical loci, an improvement of 25 loci (22.52 percentage
points). This is a numerical processability diagnostic, not an accepted
replacement analysis.

## Scientific and family gates

The pre-outcome [sample-size audit](../lava_multitrait_feasibility_v1/mdd2025_pilot_sample_size_audit_20260927.md)
identified unresolved PGC `NEFF` versus LAVA binary `N` and aggregate case
fraction semantics. The pilot used source per-variant `NEFF` as `N`, although
the source defines `NEFF` as the sum of cohort-level effective sizes. The
terminal interpretation is therefore
`DIAGNOSTIC_NUMERICAL_PASS_SCIENTIFIC_INTERPRETATION_HOLD`; no full MDD screen
is authorized from this pilot.

The immutable canonical v3 family remains `FAILED_QC_NOT_PROMOTED`:
3,720/17,465 cells are `NOT_RUN`, above the frozen maximum of 873. The
separately versioned LAVA rescue v1 is `FAIL_QC`. Even a perfect MDD
replacement would leave 3,131 `NOT_RUN` cells among the other six traits.
Even a perfect replacement of both sleep traits would leave 1,758; the
observed continuous-duration substitute alone has 876. These are arithmetic
lower bounds, not predicted replacement results. The
[bound calculation](../lava_multitrait_feasibility_v1/lower_bounds.json)
rules out a single-trait rescue under the current family rule.

The approved narrow Track B legacy admission passed the locked preflight,
full-row adapter and independent BH audit, scientific compatibility checks,
and protected-slot hash audit. The five-track PLACO pair family covers
20,370,452 tested rows and passes its correction audit. It yields 25
shared-association candidate loci (19 carried forward, six from Track B),
all blocked from final evidence tiering by the LAVA gate. These are neither
validated shared local signals nor causal loci. The
[five-track gate](../loci/five_track_candidate_region_gate_v1/decisions.tsv)
and [combined report](../brain6_unblock_parallel_v2/combined_unblock_report.md)
contain the receipt-bound details.

## Current decision

**NO-GO for final region tiering and dependent fine-mapping, colocalization,
gene, cell-type, and pathway claims.** No candidate promotion or downstream
analysis was run after the failed LAVA gate. The canonical v3, baseline v2,
roundoff, and original pilot artifacts retain their distinct output roots and
decisions.

A new rescue requires source-verified, scientifically equivalent and better
powered GWAS for multiple traits, with defensible per-variant sample-size
semantics, reference compatibility, and a prospectively frozen pilot. The
current local source inventory has no qualified set that can meet the frozen
873-cell ceiling. Any different family or QC rule would require a separate
prospective protocol; it cannot reclassify canonical v3.
