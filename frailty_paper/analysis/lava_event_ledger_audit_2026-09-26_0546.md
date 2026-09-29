# FI×sleep LAVA parallel event-ledger audit — 05:46 UTC

Audit time: 2026-09-26T05:46:47.370745+00:00. The coordinator state was 9,968 at both the start and end of the event scan. The event ledger contained 5,149 lines at scan and post-scan checks.

- Receipt events: 5,094, mapping to 5,094 unique `(trait, locus_index)` pairs.
- Duplicate receipt pairs: 0.
- Malformed JSON lines: 0.
- Other events: 54 `worker_finished`, 1 `scheduler_paused`.
- Two current coordinator claims: shortsleep loci 568 and 994 (workers 78548 and 78549); launch failures and stale recoveries remain zero.
- Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Event-stream SHA-256: `3ed7a630646dbd946ac3fca5d343365cb5049d0d0a10a44caec82c788dccbb8b`.

This append-only event-ledger check is not a full receipt-file validation. The latest full receipt-file audit available after this event scan is the 05:43 audit, which validated 9,960 files; eight later coordinator receipts were outside that audit. Machine-readable details: `lava_event_ledger_audit_2026-09-26_0546.json`.
