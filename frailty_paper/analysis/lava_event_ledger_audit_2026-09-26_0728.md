# FI×sleep LAVA parallel event-ledger audit — 07:28 UTC

Two snapshots at 2026-09-26T07:29:23+00:00 and 2026-09-26T07:29:23+00:00 were identical for the event stream and coordinator receipt count: 10,273.

- Event lines: 5,454; receipt events: 5,399; unique `(trait, locus_index)` identities: 5,399.
- Duplicate identities: 0; malformed lines: 0; event types: receipt=5399, worker_finished=54, scheduler_paused=1.
- Active claims: napping locus 583 (PID 78548); snoring locus 1009 (PID 78549); both PIDs passed liveness checks. Requested/effective workers: 6/2; launch failures: 0; stale recoveries: 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; event-stream SHA-256: `9087fae731209c43b3e4fab256f982f8ff8e19966c183d9f49b7dfe4e7ac2fe4`.

This event-ledger scan covers the current coordinator history only. The companion full receipt audit validated 10,259 files, with two later receipts added while it ran. Machine-readable details: `lava_event_ledger_audit_2026-09-26_0728.json`.
