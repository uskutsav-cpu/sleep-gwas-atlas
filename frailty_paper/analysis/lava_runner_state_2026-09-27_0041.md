# FI×sleep LAVA coordinator state — 2026-09-27 00:42:09 UTC

Read-only state snapshot. This is not a receipt-integrity audit or a fresh process-liveness check.

- Coordinator reports 18,557/29,940 verified receipts (11,383 remain by runner count).
- Requested concurrency: 6; effective concurrency: 2; active claim files: 2.
- Runner-reported active jobs: `[{"claimed_utc": "2026-09-27T00:41:50.343100+00:00", "locus_index": 1781, "pair_order": 4, "process_group": 72535, "trait": "longsleep", "worker_pid": 72535}, {"claimed_utc": "2026-09-27T00:42:04.926944+00:00", "locus_index": 1781, "pair_order": 5, "process_group": 72534, "trait": "napping", "worker_pid": 72534}]`
- Claim files: `[{"claimed_utc": "2026-09-27T00:41:50.343100+00:00", "locus_index": 1781, "pair_order": 4, "process_group": 72535, "trait": "longsleep", "worker_pid": 72535}, {"claimed_utc": "2026-09-27T00:42:04.926944+00:00", "locus_index": 1781, "pair_order": 5, "process_group": 72534, "trait": "napping", "worker_pid": 72534}]`
- Launch failures this invocation: 0; stale claims recovered: 0.
- Swap safeguard: swap use reached 3883/5120 MiB (≥3584 MiB used); reducing three workers to two; control file retains desired_workers=6.
- The controller is already at two workers after recorded downshifts. No concurrency control, runner process, claim, receipt, input or scientific setting was changed.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- From the 00:32:09 UTC state (18,500 receipts) to 00:42:09 UTC, 57 receipts arrived (~341/hour). This short interval is not a reliable completion estimate.
- Host storage sample: 2.8 GiB free on the internal volume and 1.5 TiB free on the external SSD. `memory_pressure` showed 5,559 free pages of 524,288 (16 KiB pages); process-list access and `sysctl vm.swapusage` are restricted in this environment. The runner-reported swap safeguard remains authoritative.
- The latest full receipt-integrity audit remains 00:15 UTC (18,420 validated; runner/inventory reached 18,422 while scanning; zero receipt/claim/duplicate issues). All 12 frozen 1% gates fail; no local-sharing inference is admissible.

No runner, claim, receipt, input, or analysis setting was modified.
