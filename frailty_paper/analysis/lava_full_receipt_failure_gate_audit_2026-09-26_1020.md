# FI×sleep LAVA full-family audit — 2026-09-26 10:20 UTC

The full read-only auditor validated 12,504/29,940 receipts from 10:20:21 to 10:20:26 UTC. Inventory advanced 12,504→12,506; runner state advanced 12,502→12,504.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims: longsleep locus 692 and shortsleep locus 1124. Receipt/claim/identity/duplicate issues, launch failures, and stale recoveries were zero. Six workers remain requested and two effective under the existing swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.

Compared with the 10:18:57 scan at 12,484 validated receipts, the 10:20:21 scan validated 12,504, an observed short-window rate of approximately 857 receipts/hour. This estimate is provisional; 17,434 receipts remained after inventory scan. Failure categories total 2,151 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% gates remain failed, so no local-sharing inference is admissible.

This audit is read-only and does not modify settings, inputs, thresholds, or receipts. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_1020.json`.
