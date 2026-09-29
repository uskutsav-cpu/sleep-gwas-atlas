# Frailty project checkpoint — 2026-09-27 22:05 UTC

Coordinator state at 2026-09-27T22:05:45.422378+00:00: 21,756/29,940 receipts, two active claims (PID 15810 longsleep locus 2097; PID 15811 napping locus 2097), zero launch failures, zero stale recoveries. Full audit 2026-09-27T22:04:58+00:00–2026-09-27T22:05:16+00:00 validated 21,750; inventory reached 21,752; receipt/claim/duplicate issues zero, both workers live. All 12 locked gates fail (3,790 negative variance, 252 no-reference-SNP, 2 other errors).

The latest 30-minute receipt-event window contains 115 distinct trait-locus completions (230/hour) with 0 duplicate events; the trailing 15-minute rate is 280/hour. The coordinator's original post-reattachment rate is 439.58/hour. 8,184 receipt slots remain, implying ~35.6 hours at the 30-minute rate, provisional. Event window/hash are preserved in `lava_rolling_throughput_2026-09-27_2205.json`.

Figure 4 was refreshed from this audit. Host memory free was 55%; external SSD free space is 1.5 TiB. CPU and swap-total sampling remain restricted. Frozen lock/input hashes unchanged. The two reviewer archives are ready; delivery destination remains pending.
