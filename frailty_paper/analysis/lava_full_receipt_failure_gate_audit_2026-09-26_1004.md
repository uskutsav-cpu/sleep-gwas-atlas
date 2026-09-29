# FI×sleep LAVA full-family audit — 2026-09-26 10:04 UTC

The full read-only auditor validated 12,284/29,940 receipts from 10:04:18 to 10:04:22 UTC. Inventory advanced 12,284→12,285 during the scan; runner state remained at 12,283.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims: insomnia locus 1113 and shortsleep locus 681. Receipt/claim/identity/duplicate issues, launch failures, and stale recoveries were zero. Six workers remain requested and two effective under the existing swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.

From the 09:58:46–09:58:51 audit's 12,201 validated receipts to this audit's 12,284, observed throughput was approximately 912 receipts/hour. This is provisional; 17,656 receipts remained. Failure categories total 2,113 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% gates remain failed, so no local-sharing inference is admissible.

This audit is read-only and does not modify settings, inputs, thresholds, or receipts. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_1004.json`.
