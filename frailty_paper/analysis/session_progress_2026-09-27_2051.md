# Frailty project checkpoint — 2026-09-27 20:51 UTC

Runner state at 2026-09-27T20:54:29.315191+00:00 reports 21,366/29,940 receipts, two effective workers (six requested), zero launch failures and zero stale claims. The last rolling 30-minute rate (512.70/hour) spans the 20:38:53 downshift; it is not a steady two-worker rate. 8,574 receipt slots remain; no ETA is asserted.

The full read-only audit ran 20:51:16–20:52:27 (71 seconds), validated 21,354 receipts and observed inventory advance to 21,358; receipt, claim and duplicate checks were clean and both workers were live. All 12 gates fail (3,723 negative-variance, 252 no-reference-SNP, two other failures). Evidence: `lava_full_receipt_audit_2026-09-27_2051.json`, `lava_runner_state_snapshot_2026-09-27_2051.json`. Figure 4 was refreshed from this scan. To reduce added I/O during the safeguard downshift, the next full scan will be in about 10 minutes; runner-state checks continue between audits.

The frozen 56,117-record screening queue has zero decisions; two separate reviewer archives are ready and verified. Human screening/adjudication, licensed-source exports/full texts, exact HFRS access and exact cohort intersections remain open.
