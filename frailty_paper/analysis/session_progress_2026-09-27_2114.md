# Frailty project checkpoint — 2026-09-27 21:14 UTC

The reattached LAVA coordinator reports 21,536/29,940 receipts at 2026-09-27T21:16:11.830904+00:00, with two active distinct worker claims, zero launch failures and zero stale-claim recoveries. Its throughput window reset during reattachment and remains unset; 8,404 slots remain, with no ETA yet.

The independent full audit at 21:14:16–21:15:02 validated 21,519 receipts; inventory advanced to 21,526 during scanning. Receipt, claim and duplicate checks were clean; both PIDs 15810/15811 were live. All 12 gates fail (3,733 negative-variance, 252 no-reference-SNP, two other failures). Evidence: `lava_full_receipt_audit_2026-09-27_2114.json`, `lava_runner_state_snapshot_2026-09-27_2114.json`. Figure 4 was updated from this audit. The 46-second scan is followed by a 10-minute full-audit interval to limit added I/O; runner-state monitoring continues.

The two unclassified process failures remain documented with hashes; neither was retried. The screening queue remains 56,117 records with zero decisions and two reviewers ready. Licensed search exports/full texts, human decisions/adjudication, exact HFRS access and exact cohort intersections remain outstanding.
