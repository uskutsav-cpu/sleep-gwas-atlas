# FI×sleep LAVA full-family audit — 2026-09-26 09:54 UTC

The full read-only auditor validated 12,160/29,940 receipts from 09:54:41 to 09:54:46 UTC. Inventory advanced 12,160→12,162 during the scan; runner state remained 12,160 at both reads.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims: sleep_timing locus 675 and sleep_apnea locus 1106. Receipt/claim/identity/duplicate issues, launch failures, and stale recoveries were zero. Six workers remain requested and two effective under the existing swap safeguard. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.

From the 09:51 audit's 12,119 validated receipts to this audit's 12,160, observed throughput was approximately 780 receipts/hour. This short-window estimate is provisional; 17,778 receipts remained after the scan. Failure categories total 2,074 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% gates remain failed, so no local-sharing inference is admissible.

This audit is read-only and does not modify settings, inputs, thresholds, or receipts. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_0954_41_recheck.json`.
