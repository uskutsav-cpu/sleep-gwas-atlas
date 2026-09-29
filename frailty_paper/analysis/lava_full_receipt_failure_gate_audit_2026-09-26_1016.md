# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 10:16 UTC

Read-only full-family audit of the existing FI×sleep LAVA run. No settings, claims, receipts, thresholds, or inputs were changed.

- Audit window: 10:16:18–10:16:24 UTC.
- Validated 12,446 of 29,940 receipt files (41.60%). The inventory grew from 12,446 to 12,447 while scanning; runner state advanced from 12,445 to 12,447.
- Receipt integrity, logs, claims, duplicate claim identities, and duplicate trait/locus identities: zero issues. Both active worker PIDs (78548, 78549) were live on distinct claims. Launch failures and stale-claim recoveries: zero.
- Scheduler requests six workers; two are effective because swap usage reached 2,772/3,072 MiB (90%), triggering the existing resource safeguard.
- Frozen failure gate: all 12 traits fail the 1% maximum process-failure criterion. Failures total 2,384: 2,142 all-phenotype negative-variance and 242 no-specified-SNPs-in-reference. Local-sharing inference and PLACO remain unjustified.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Ordered receipt hash-stream SHA-256: `23bd94407bf3d85dc36c23d2d6eeab120f2c9e884e67bd09cb1d2665b1ec4c81`.

Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1016.json`. The audit is an advancing snapshot, not a completion check; continue the existing run only under the persisted safeguards and retain the frozen gate.
