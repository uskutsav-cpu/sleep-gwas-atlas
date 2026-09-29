# FI×sleep LAVA live runner state — 2026-09-26 23:18 UTC

Read-only snapshot of the externally stored coordinator state at 2026-09-26 23:18:39 UTC. This records scheduler-reported state only; no OS process-table or fresh CPU/RAM/swap telemetry was available in this snapshot.

- Verified receipts: 18,269 / 29,940; 11,671 not yet verified.
- Requested/effective workers: 6 / 2. The runner reports swap at 3,883 / 5,120 MiB and the existing resource guard is holding concurrency at two.
- Active claims: `shortsleep` locus 1653 (worker PID 72534) and `sleep_apnea` locus 1653 (worker PID 72535). These are distinct trait/locus identities.
- Launch failures this invocation: 0. Stale claims recovered: 0.
- Frozen analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Coordinator-reported PID map: worker slots 1–4 have PIDs 72534–72537; only two workers are effective. Do not treat the PID map as proof that all four slots are currently active.

No scheduler, claim, locus, input, parameter, threshold, or output change was made. In particular, the six-worker setting was not forced past the established swap safeguard.
