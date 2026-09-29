# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 07:17 UTC

Read-only scan interval: 2026-09-26T07:17:46+00:00 to 2026-09-26T07:17:51+00:00.

- Validated receipts: 10,227/29,940 (34.16%); receipt files and runner state matched at both scan boundaries.
- Structural, identity, claim, duplicate, and matching-log issues: zero.
- Active jobs: 2 live PIDs (78548, 78549); effective workers 2; requested workers 6 under the existing 90%-swap fallback.
- Launch failures and stale recoveries: zero.
- Process failures: 1,955 (1,713 all-phenotype negative variance; 242 no reference SNPs). All 12 traits fail the locked 1% family gate.
- Locked analysis SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- This audit did not mutate run outputs, receipts, thresholds, or scheduler state.

Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_0717.json`.
