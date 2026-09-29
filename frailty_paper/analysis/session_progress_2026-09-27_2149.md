# Frailty project throughput update — 2026-09-27 21:49 UTC

The coordinator reports 21,680/29,940 receipts, two active workers and zero launch failures/stale claims. Latest full receipt audit remains the 21:44:36–21:45:10 UTC scan, which validated 21,670; inventory reached 21,672, with no receipt, claim or duplicate issues. All 12 locked gates fail.

The coordinator's first post-reattachment fixed window was 439.58 receipts/hour. Recounting immutable `receipt` events in `parallel_events.jsonl` for 2026-09-27T21:19:17.173033+00:00–2026-09-27T21:49:17.173033+00:00 gives 128 distinct trait-locus completions (256/hour), with zero duplicated events; the latest 15 minutes had 48 completions (192/hour). At 8,260 remaining receipts, the 30-minute rate implies ~32.3 hours, provisional. Evidence: `lava_rolling_throughput_2026-09-27_2149.json`.

The two-worker safeguard remains active following earlier severe swapping. Current host memory free is 40%; CPU and swap-total sampling remain unavailable from this restricted shell. Frozen input/lock hashes remain unchanged. Two independent reviewers are ready, but their delivery route is still pending.
