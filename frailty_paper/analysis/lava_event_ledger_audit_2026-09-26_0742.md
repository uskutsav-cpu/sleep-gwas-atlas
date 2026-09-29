# FI×sleep LAVA event-ledger audit — 2026-09-26 07:42 UTC

Two reads of the append-only coordinator event ledger were identical at 2026-09-26T07:43:12.371473+00:00.

- Event lines: 5,646; event types: receipt=5591, scheduler_paused=1, worker_finished=54.
- Receipt events: 5,591; unique `(trait,locus_index)` identities: 5,591; duplicates: 0; malformed lines: 0.
- Coordinator state: 10,461/29,940; two active jobs (longsleep locus 592, sleepiness locus 1019); requested/effective workers 6/2; launch failures 0; stale recoveries 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Event-stream SHA-256: `f72462f85ba6c3631d69510ec641b4349aab71be7b3fb7c5de4b3c98e6e49b54`.

This ledger covers receipt events recorded by the current coordinator and is not a complete validation of all on-disk receipt files. The latest full receipt audit is the 07:38 snapshot. Machine-readable details: `lava_event_ledger_audit_2026-09-26_0742.json`.
