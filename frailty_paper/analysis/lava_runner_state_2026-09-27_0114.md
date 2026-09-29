# FI×sleep LAVA coordinator state — 2026-09-27T01:14:23Z

Read-only state snapshot; not a full receipt-integrity audit.

- Runner count: 18,723/29,940; 11,217 remain.
- Six workers requested; two effective. Coordinator PID 72428 and worker PIDs 72534/72535 were present in the process table.
- Live claims at 01:13:52–01:13:56 UTC: `sleepdur` locus 1797 (PID 72535) and `sleepiness` locus 1797 (PID 72534), distinct claims.
- Launch failures this invocation: 0; stale claims recovered: 0.
- Persisted fallback: swap use reached 3,883/5,120 MiB (at least 3,584 MiB), reducing three workers to two. Direct host sample at 01:14 UTC: 5,021.38/6,144 MiB swap used (1,122.62 MiB free). The existing resource fallback remains necessary; six workers are requested, but only two are effective.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Latest full receipt audit remains 00:45:59–00:46:24 UTC: 18,575 validated, runner reached 18,577 during scan, zero receipt/claim/duplicate issues; all 12 frozen 1% gates failed.
- No runner, worker-control, claim, receipt, input, or scientific setting was modified.
