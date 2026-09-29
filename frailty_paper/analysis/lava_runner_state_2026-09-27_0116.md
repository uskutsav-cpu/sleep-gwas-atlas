# FI×sleep LAVA coordinator state — 2026-09-27T01:16:24.965857+00:00

Read-only state snapshot; not a full receipt-integrity audit.

- Runner count: 18,736/29,940; 11,204 remain.
- Six workers requested; 2 effective. Active worker PIDs [72534, 72535] passed liveness checks.
- Active claims: `[{"claimed_utc": "2026-09-27T01:16:23.530618+00:00", "locus_index": 1799, "trait": "longsleep", "worker_pid": 72535}, {"claimed_utc": "2026-09-27T01:16:25.834092+00:00", "locus_index": 1799, "trait": "napping", "worker_pid": 72534}]`.
- Launch failures this invocation: 0; stale claims recovered: 0.
- Persisted fallback: swap use reached 3883/5120 MiB (≥3584 MiB used); reducing three workers to two. Direct host sample at 01:14 UTC: 5,021.38/6,144 MiB swap used (1,122.62 MiB free). The existing resource fallback remains necessary; six workers are requested, but only two are effective.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Latest full receipt audit: 01:16:19–01:16:29 UTC; 18,735 validated, two arrived during scan, zero receipt/claim/duplicate issues. All 12 frozen 1% gates fail.
- No runner, worker-control, claim, receipt, input, or scientific setting was modified.
