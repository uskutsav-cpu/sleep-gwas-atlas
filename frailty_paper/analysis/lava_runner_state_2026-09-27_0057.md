# FI×sleep LAVA coordinator state — 2026-09-27T00:57:46.940767+00:00

Read-only coordinator snapshot. This is not a full receipt-integrity audit.

- Runner count: 18,639/29,940; 11,301 remain by runner state.
- Six workers requested; 2 effective; `concurrency.json` retains desired_workers=6.
- Active process groups/locus claims: [{"claimed_utc": "2026-09-27T00:57:29.951111+00:00", "locus_index": 1789, "pair_order": 6, "process_group": 72535, "process_group_alive": true, "trait": "shortsleep", "worker_pid": 72535}, {"claimed_utc": "2026-09-27T00:57:36.610139+00:00", "locus_index": 1789, "pair_order": 7, "process_group": 72534, "process_group_alive": true, "trait": "sleep_apnea", "worker_pid": 72534}]
- Launch failures this invocation: 0; stale claims recovered: 0.
- Existing resource fallback: swap use reached 3883/5120 MiB (≥3584 MiB used); reducing three workers to two.
- Frozen analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` (embedded lock SHA `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`).
- No runner, control, claim, receipt, input, or scientific setting was modified by this check.
- Latest full receipt-integrity audit: 2026-09-27 00:45:59–00:46:24 UTC, 18,575 receipts validated, two arrivals during scan; zero receipt/claim/duplicate issues. All 12 frozen 1% gates failed.
