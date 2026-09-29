# FI×sleep LAVA coordinator state — 2026-09-27T01:50:04+00:00

Read-only state snapshot at the end of the full receipt audit; not a separate process-command inspection.

- Runner state: 19,011/29,940; 10,929 remain.
- Six workers requested; 2 effective. Coordinator state refreshed at 2026-09-27T01:50:04.272997+00:00.
- Active claims: sleep_efficiency locus 1826 (PID 72534, alive=True), sleep_timing locus 1826 (PID 72535, alive=True).
- Launch failures this invocation: 0; stale claims recovered: 0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Adaptive resource fallback: swap use reached 3883/5120 MiB (≥3584 MiB used); reducing three workers to two. CPU and direct process-command telemetry were unavailable in this sandbox; the full auditor's PID liveness check returned true for both active worker PIDs.
- Full read-only receipt audit: 19,009 valid receipts, zero receipt/claim/duplicate issues; all 12 frozen 1% gates failed.
- No runner, worker-control, claim, receipt, input, or scientific setting was modified.
