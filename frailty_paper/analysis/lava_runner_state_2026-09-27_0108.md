# FI×sleep LAVA coordinator state — 2026-09-27T01:08:38.149114+00:00

Read-only state snapshot; not a full receipt-integrity audit.

- Runner count: 18,699/29,940; 11,241 remain.
- Six workers requested; 2 effective; live claims: `[{"claimed_utc": "2026-09-27T01:08:21.540201+00:00", "locus_index": 1795, "pair_order": 6, "process_group": 72535, "process_group_alive": true, "trait": "shortsleep", "worker_pid": 72535}, {"claimed_utc": "2026-09-27T01:08:32.808559+00:00", "locus_index": 1795, "pair_order": 7, "process_group": 72534, "process_group_alive": true, "trait": "sleep_apnea", "worker_pid": 72534}]`.
- Launch failures this invocation: 0; stale claims recovered: 0.
- Persisted fallback reason: swap use reached 3883/5120 MiB (≥3584 MiB used); reducing three workers to two; control file desired_workers=6.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Latest full receipt audit is 00:45:59–00:46:24 UTC: 18,575 validated, runner reached 18,577 during scan, zero receipt/claim/duplicate issues; all 12 frozen 1% gates failed.
- No runner, worker-control, claim, receipt, input, or scientific setting was modified.
