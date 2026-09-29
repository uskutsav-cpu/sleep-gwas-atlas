# Frailty project checkpoint — 2026-09-27 20:38 UTC

The existing FI × sleep LAVA run remains active. At 2026-09-27T20:41:11.351423+00:00, the runner reports 21,314/29,940 receipts and two effective workers (six requested), with zero launch failures or stale claims. The last rolling 30-minute rate is 711.20/hour, measured before the current downshift; it is not a valid new two-worker rate. The rough 12.1-hour remainder estimate is therefore provisional. A resource guard automatically downshifted three to two workers when swap use reached 6,310/7,168 MiB (5,800 MiB threshold). Active distinct claims remain on PIDs 15810 and 15811.

The full audit at 20:38:39–20:38:52 validated 21,301 receipts; inventory reached 21,302, with no receipt/claim/duplicate issues. All 12 gates fail (3,719 negative-variance, 252 no-reference-SNP, two other failures). Evidence: `lava_full_receipt_audit_2026-09-27_2038.json`, `lava_runner_state_snapshot_2026-09-27_2038.json`. The two unclassified failures are recorded separately; no receipts were edited or retried. Figure 4 was regenerated from this audit.

The title/abstract queue remains 56,117 records with zero decisions; two reviewer archives are prepared/verified. Human review, licensed-source exports/full texts, exact HFRS access and exact cohort intersections remain outstanding.
