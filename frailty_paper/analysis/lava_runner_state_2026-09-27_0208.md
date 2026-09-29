# FI×sleep LAVA coordinator state — 2026-09-27T02:08:18.971017+00:00

Read-only state snapshot at the end of the full receipt audit; not a separate command-line audit.

- Runner state: 19,117/29,940; 10,823 remain.
- Six workers requested; 2 effective.
- Active claims: napping locus 1837 (PID 72534, alive=True), shortsleep locus 1837 (PID 72535, alive=True).
- Launch failures this invocation: 0; stale claims recovered: 0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Adaptive resource fallback: swap use reached 3883/5120 MiB (≥3584 MiB used); reducing three workers to two. Auditor PID-liveness checks returned true for both active workers; `ps`/`top` command-line telemetry was denied by the sandbox.
- Full read-only receipt audit: 19,116 valid receipts, zero receipt/claim/duplicate issues; all 12 frozen 1% gates failed.
- No runner, worker-control, claim, receipt, input, or scientific setting was modified.
