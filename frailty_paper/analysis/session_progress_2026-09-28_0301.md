# Frailty LAVA checkpoint — 2026-09-28 03:01 UTC

Coordinator state at 2026-09-28T03:00:56.617958+00:00: 22,010/29,940, two active claims, zero launch/stale-claim errors. Full audit 2026-09-28T03:00:02+00:00–2026-09-28T03:00:40+00:00 validated 22,006; inventory reached 22,009; receipt/claim/duplicate issues zero; both workers live. All 12 locked gates fail (3,820 negative variance; 252 no-reference-SNP; 2 other failures).

After the earlier interruption, receipt events have been continuous since 2026-09-28T02:16:02.285836+00:00; the largest interreceipt gap since resumption is 50.3 seconds. The latest rolling rates are 230/hour over 30 minutes (115 events) and 188/hour over 15 minutes; the coordinator rate is 235.44/hour. With 7,930 remaining slots, the 30-minute rate implies ~34.5 hours of continuous runtime, provisional; wall-clock timing depends on the host remaining active. Evidence: `lava_rolling_throughput_2026-09-28_0301.json`.

Figure 4 was refreshed. Host memory free is 54%; external SSD has 1.5 TiB free. CPU and swap-total samples remain unavailable in the restricted shell. Frozen lock/input hashes unchanged. Reviewer delivery remains unresolved.
