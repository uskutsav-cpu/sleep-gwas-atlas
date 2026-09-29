# FI×sleep LAVA full-family audit — 2026-09-26 10:12 UTC

The full read-only auditor validated 12,401/29,940 receipts from 10:12:54 to 10:12:59 UTC. Inventory advanced 12,401→12,402 during the scan; runner state advanced 12,399→12,402.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims: longsleep locus 1119 and snoring locus 686. Receipt/claim/identity/duplicate issues, launch failures, and stale recoveries were zero. Six workers remain requested and two effective under the existing swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.

From the 10:07:47 audit's 12,332 runner-verified receipts to this audit's 12,401 validated receipts, observed throughput was approximately 805 receipts/hour. This is provisional; 17,538 receipts remained after scan-end inventory. Failure categories total 2,136 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% gates remain failed, so no local-sharing inference is admissible.

This audit is read-only and does not modify settings, inputs, thresholds, or receipts. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_1013.json`.
