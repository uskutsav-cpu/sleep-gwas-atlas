# FI×sleep LAVA coordinator state — 2026-09-27 00:32 UTC

Read-only coordinator-state snapshot; this is not a full receipt-integrity audit.

- Updated at 2026-09-27T00:32:09 UTC: 18,500/29,940 verified receipts; 11,440 remain by coordinator count.
- Six workers requested; two effective. The scheduler reports the existing swap safeguard at 3,883/5,120 MiB (downshift threshold 3,584 MiB used).
- Active claims at snapshot: sleep_apnea locus 1775 (worker PID 72534) and sleep_efficiency locus 1775 (worker PID 72535). No fresh process-liveness check was made in this state-only read.
- Launch failures this invocation: 0. Stale claims recovered: 0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- From the 00:26:12 snapshot at 18,465 to 00:32:09 at 18,500, 35 receipts arrived in 5m57s (~353/hour). This short-interval rate is not a reliable completion estimate.
- Host capacity sample at 00:32 UTC: 2.8 GiB internal free space and 1.5 TiB external free space. Earlier memory-pressure sample at 00:23 UTC reported 38% system-wide free; `sysctl vm.swapusage` was denied. The runner's 3,883/5,120 MiB swap reading continues to govern its concurrency safeguard.
- The latest full receipt-integrity audit remains 00:15 UTC (18,420 validated, inventory/runner reached 18,422, no receipt/claim/duplicate issues; both worker PIDs passed liveness then). All 12 frozen 1% trait gates fail; no local-sharing inference is admissible.

Sources: external `runner_state.json` and append-only `parallel_events.jsonl` under `.../lava_sensitivity_v1/`, plus read-only host samples described above. No runner, claim, receipt, input, or analysis setting was changed.
