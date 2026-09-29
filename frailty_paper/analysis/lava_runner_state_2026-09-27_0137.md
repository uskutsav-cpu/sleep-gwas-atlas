# FI×sleep LAVA coordinator state — 2026-09-27T01:37:12Z

Read-only state snapshot; not a full receipt-integrity audit.

- Runner state: 18,913/29,940; 11,027 remain.
- Six workers requested; two effective. Coordinator PID 72428 and workers 72534/72535 are live.
- Active claims at the sample: worker 72534 on `sleepdur` locus 1816 and worker 72535 on `sleepiness` locus 1816. Their R child command lines matched these claims.
- Launch failures this invocation: 0; stale claims recovered: 0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc433db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Direct swap sample: 5,048.50/6,144 MiB used. The controller's existing fallback remains necessary; no worker settings were changed.
- Latest full receipt audit: 01:35:25–01:35:43 UTC, validated 18,899 files, stable inventory/state and zero integrity/claim/duplicate issues. All 12 frozen 1% gates failed.
