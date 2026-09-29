# FI×sleep LAVA parallel event-ledger audit — 07:07 UTC

Two snapshots were captured at 2026-09-26T07:08:04+00:00 and 2026-09-26T07:08:04+00:00; stable across reads: True. Coordinator state was 10,203 at both snapshots.

- Event lines: 5,384; receipt events: 5,329; unique `(trait, locus_index)` receipt identities: 5,329.
- Duplicate receipt identities: 0; malformed JSON lines: 0; event types: receipt=5329, worker_finished=54, scheduler_paused=1.
- Active claims: sleep_efficiency locus 1006 (PID 78549); sleep_timing locus 579 (PID 78548); PID liveness: {'78549': True, '78548': True}. Requested/effective workers: 6/2; launch failures: 0; stale recoveries: 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; event-stream SHA-256: `54e6a8dfcc34054372396cf72ff47ad55eea8517464c7cc576b4afbaacaf3d21`.

This audits only the current append-only coordinator event history; it is not a full receipt-file audit. The latest full receipt-file audit validated 10,194 files at 2026-09-26T07:04:50+00:00. Machine-readable details: `lava_event_ledger_audit_2026-09-26_0707.json`.
