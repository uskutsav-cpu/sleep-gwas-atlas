# Frailty project checkpoint — 2026-09-27 21:55 UTC

Coordinator state at 2026-09-27T21:55:32.889600+00:00: 21,704/29,940 receipts, two active claims, zero launch failures, and zero stale recoveries. Full read-only audit 2026-09-27T21:54:43+00:00–2026-09-27T21:55:04+00:00 validated 21,700; inventory reached 21,702; no receipt, claim or duplicate issues; both workers live. All 12 locked gates fail (3,774 negative variance, 252 no-reference-SNP, 2 other errors).

The latest 30-minute event window had 118 distinct receipt events (236/hour) and no duplicate events; the latest 15-minute segment had 50 events (200/hour). The coordinator's original reattachment window remains 439.58/hour. 8,236 slots remain; the rolling 30-minute rate implies ~34.9 hours, provisional. Machine evidence: `lava_rolling_throughput_2026-09-27_2155.json`.

Figure 4 was regenerated. Host memory free was 63%, external SSD free space 1.5 TiB, and CPU/swap-total sampling remains blocked in the restricted shell. Frozen lock/input hashes unchanged. The two independent reviewers are ready; handoff destination is pending.
