# FI×sleep LAVA parallel event-ledger audit — 06:09 UTC

Audit time: 2026-09-26T06:10:15+00:00. The coordinator state was 10,036 for both snapshots; event file stable: True; runner count stable: True.

- Append-only event lines: 5,217; receipt events: 5,162, all unique `(trait, locus_index)` identities: 5,162.
- Duplicate receipt identities: 0; malformed JSON lines: 0. Event types: receipt=5162, worker_finished=54, scheduler_paused=1.
- Active claims: sleepdur locus 571 (worker PID 78548); sleepdur locus 997 (worker PID 78549). Requested/effective workers: 6/2; launch failures: 0; stale recoveries: 0.
- Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Event-stream SHA-256: `67c08a644d9c96fd69998830f0b1ab2caf4e6169d6bbf7822533cfec09dd66ee`.

This event-ledger scan covers only the current parallel-resume event history; it is not a full validation of the carried-forward receipt-file family. The companion full receipt audit validated 10,031 files at 2026-09-26T06:09:05+00:00. Machine-readable details: `lava_event_ledger_audit_2026-09-26_0609.json`.
