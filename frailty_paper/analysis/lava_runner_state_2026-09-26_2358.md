# FI×sleep LAVA live runner state — 2026-09-26 23:58 UTC

Read-only scheduler and host sample. This is not a full receipt-integrity audit; no fresh OS process-table liveness check was available.

- Runner state updated 2026-09-26 23:58:33 UTC: 18,385 / 29,940 receipts (11,555 remaining).
- Six workers requested; two effective under the existing memory-pressure fallback.
- Active scheduler claims: sleepiness locus 1763 (worker PID 72535) and snoring locus 1763 (worker PID 72534). The claim identities are distinct; this state file alone does not prove process liveness.
- Launch failures: 0. Stale claims recovered: 0.
- Runner-reported swap trigger remains 3,883 / 5,120 MiB; concurrency remains capped at two.
- Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Host sample: 8,622 free pages at 16 KiB/page (~134 MiB); compressor occupied 122,782 pages (~1.87 GiB). A direct `sysctl vm.swapusage` read was denied. External SSD had 1.5 TiB free. Internal-volume latest check remains 2.9 GiB at 23:18 UTC.

At the 23:58 sample, receipts had increased by 116 from the 23:18 state snapshot (40 minutes, ~174/hour; a noisy short-window rate). The 23:07 full integrity audit remains the latest audited checkpoint. No extra LDSC job was started under this memory pressure; no LAVA configuration, inputs, claims, or scientific settings were changed.
