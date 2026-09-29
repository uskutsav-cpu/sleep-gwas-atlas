# FI×sleep LAVA full-family audit — 2026-09-26 09:30 UTC

The full read-only auditor validated 11,842/29,940 receipt files from 09:30:20 to 09:30:25 UTC. Inventory was stable at 11,842 through the scan. Both active worker PIDs (78548, 78549) passed liveness checks and were assigned distinct claims. Receipt, claim, duplicate-identity, and duplicate-claim checks were clean; launch failures and stale recoveries were zero.

The campaign remains at 6 requested / 2 effective workers. The persisted adaptive fallback says: “swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two”. Lock SHA-256 is `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 is `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

From the 09:25 audit (11,783 validated at 09:25:22) to this audit (11,842 validated at 09:30:20), throughput was approximately 725 receipts/hour over a short interval. This is provisional; approximately 18,098 receipts remain. Process-failure categories are 2,017 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% failure gates still fail; no local-sharing inference is admissible.

No inputs, analysis settings, thresholds, or existing receipts were modified by this audit. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_0930.json`.
