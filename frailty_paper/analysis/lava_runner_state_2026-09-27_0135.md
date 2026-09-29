# FI×sleep LAVA coordinator state — 2026-09-27T01:35:41.549668+00:00

Read-only state snapshot; not a full receipt-integrity audit.

- Runner state: 18,899/29,940; 11,041 remain.
- Six workers requested; 2 effective. Active worker PIDs [72534, 72535] passed liveness checks.
- Active claims: `[{"claimed_utc": "2026-09-27T01:35:20.534184+00:00", "locus_index": 1815, "trait": "shortsleep", "worker_pid": 72534}, {"claimed_utc": "2026-09-27T01:35:21.879448+00:00", "locus_index": 1815, "trait": "sleep_apnea", "worker_pid": 72535}]`.
- Launch failures this invocation: 0; stale claims recovered: 0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Latest full receipt audit: 2026-09-27T01:35:25+00:00–2026-09-27T01:35:43+00:00 UTC; 18,899 validated, no receipt/claim/duplicate issues. All 12 frozen 1% gates fail.
- Direct 01:35 UTC swap sample: 4,546/5,120 MiB used. The existing adaptive fallback remains at two workers.
- No runner, worker-control, claim, receipt, input, or scientific setting was modified.
