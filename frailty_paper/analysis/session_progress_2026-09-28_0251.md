# Frailty LAVA resumed-run checkpoint — 2026-09-28 02:51 UTC

Coordinator state at 2026-09-28T02:51:49.477068+00:00: 21,983/29,940, two active claims, zero launch failures and stale recoveries. Full audit 2026-09-28T02:50:06+00:00–2026-09-28T02:50:56+00:00 validated 21,979; inventory reached 21,980; receipt/claim/duplicate issues zero and both PIDs live. All 12 gates fail (3,812 negative variance; 252 no-reference-SNP; 2 other failures).

The uninterrupted resumed interval is now ~35 minutes. Latest rolling event rates are 264/hour over 30 minutes and 240/hour over 15 minutes; the coordinator's completed window is 280.96/hour. The two measures are close. 7,957 receipts remain, implying ~30.1 hours of continuous runtime at the 30-minute rate, provisional. Current event calculation: `lava_rolling_throughput_2026-09-28_0251.json`.

Figure 4 was refreshed. Host memory free was 32%, external SSD free space is 1.5 TiB; current CPU and swap totals cannot be sampled in the restricted shell. The run's historical inactivity gaps mean wall-clock completion depends on the host staying active. Frozen lock/input hashes remain unchanged; reviewer delivery information is still pending.
