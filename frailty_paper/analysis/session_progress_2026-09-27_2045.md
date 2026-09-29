# Frailty project checkpoint — 2026-09-27 20:45 UTC

The existing FI × sleep LAVA run remains active on the resource safeguard's two-worker setting (six requested). At 2026-09-27T20:46:45.577129+00:00, runner state reports 21,341/29,940 receipts, zero launch failures and zero stale claims. The most recent rolling rate field remains 711.20/hour over 30 minutes, but the window overlaps the three-to-two downshift at 20:38:53, so it is not a steady-state two-worker rate. 8,599 receipt slots remain; no ETA is asserted.

Full audit at 20:45:20–20:45:47 validated 21,336 receipts; inventory reached 21,337, receipt/claim/duplicate checks were clean, both workers live, and all 12 gates fail (3,720 negative-variance, 252 no-reference-SNP, two other failures). Audit/state: `lava_full_receipt_audit_2026-09-27_2045.json`, `lava_runner_state_snapshot_2026-09-27_2045.json`. The current Figure 4 was regenerated from this audit. The two unclassified process failures remain documented; no retry or mutation was made.

The 56,117-record review queue remains undecided; two distinct reviewer archives are ready. Licensed search exports/full texts, human screening/adjudication, exact HFRS access and exact cohort intersections remain outstanding.
