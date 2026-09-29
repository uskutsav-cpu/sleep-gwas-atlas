# FI×sleep LAVA full-family audit — 2026-09-26 09:45 UTC

The full read-only auditor validated 12,035/29,940 receipts from 09:45:29 to 09:45:34 UTC. Inventory advanced 12,035→12,036; runner state advanced 12,033→12,035.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims. Receipt/claim/identity/duplicate issues, launch failures, and stale recoveries were zero. Six workers remain requested and two effective under the existing swap safeguard: “swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two”. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.

Relative to the 09:41:58 scan (11,988 validated), this scan began at 12,035, an observed short-window rate of approximately 802 receipts/hour. This is provisional; 17,904 receipts remained after scanning. Failure categories total 2,051 all-phenotype negative-variance and 242 no-reference-SNP. All 12 frozen 1% gates remain failed; no local-sharing inference is admissible.

This audit is read-only and does not modify settings, inputs, thresholds, or receipts. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_0945.json`.
