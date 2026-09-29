# Live FI×sleep LAVA checkpoint audit — 2026-09-26 03:06 UTC

## Live run state

The coordinator session handle `27134` remained pollable at this checkpoint. Its external state file, `runner_state.json`, was updated at 2026-09-26T03:06:21.439916+00:00 and reported 9,580 verified receipts of 29,940 expected, two effective workers, zero launch failures in this invocation, and zero stale claims recovered. The requested count remains six; the adaptive guard is at two workers after its recorded 2,772/3,072 MiB swap threshold. The current two claims are napping locus 549 (PID 78548) and sleep_efficiency locus 974 (PID 78549); they are distinct `(trait, locus)` jobs.

At the same snapshot, 9,580 per-locus receipt JSON files matched the coordinator receipt count. The append-only `parallel_events.jsonl` parsed without malformed JSON and contained 4,706 receipt events with 4,706 unique `(trait, locus_index)` keys. The event log is scoped to this coordinator history and is not expected to contain the carried-forward receipts from before its current invocation. The coordinator reports no stale-claim recovery or launch failure.

## Throughput and resources

From 9,556 receipts at 02:42:57 UTC to 9,580 at 03:06:21 UTC, 24 receipts completed in 23 minutes 24 seconds, approximately 62 receipts/hour. The remaining 20,360 receipts would take about 13.8 days at that short-window rate; this is a rough, highly provisional estimate because per-locus runtimes vary.

The host-wide memory sample showed 37% free and 3,979 free 16-KiB pages (about 62 MiB); the external SSD had 1.6 TiB free. The run remains at two workers under the already-triggered adaptive swap guard. The locked analysis SHA-256 is `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; the prepared-input manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. No scientific parameters, loci, QC thresholds, or inputs changed.

This is a live incremental checkpoint, not a complete receipt-manifest rehash or final LAVA result. The locked first-pair and chronotype QC failures still prohibit local-sharing inference; do not promote local results before the full scheduled family and locked QC complete.
