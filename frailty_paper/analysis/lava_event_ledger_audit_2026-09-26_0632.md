# FI×sleep LAVA parallel event-ledger audit — 06:32 UTC

Audit time: 2026-09-26T06:33:17+00:00. Coordinator state was 10,103 at both snapshots; the event file and runner count were stable.

- Event lines: 5,285; receipt events: 5,230, mapping to 5,230 unique `(trait, locus_index)` pairs.
- Duplicate receipt identities: 0; malformed JSON lines: 0; event types: receipt=5230, worker_finished=54, scheduler_paused=1.
- Active claims: sleep_apnea locus 1001 (worker 78549); sleepdur locus 574 (worker 78548). Requested/effective concurrency: 6/2; launch failures: 0; stale recoveries: 0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; event-stream SHA-256: `d8302a9530e527b594d1fead3239646f8d6c49895dd71764b98a3cb35dd998df`.

This append-only event-ledger check covers the current coordinator history and does not validate every carried-forward receipt file. The companion full audit validated 10,101 files through 2026-09-26T06:32:25+00:00. Machine-readable details: `lava_event_ledger_audit_2026-09-26_0632.json`.
