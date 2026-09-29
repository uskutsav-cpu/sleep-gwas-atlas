# FI×sleep LAVA full-family audit — 2026-09-26 10:01 UTC

The full read-only auditor validated 12,240/29,940 receipts from 10:01:34 to 10:01:40 UTC. Inventory advanced 12,240→12,242; runner state remained 12,239 at both scan boundaries.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims. Receipt/claim/identity/duplicate issues, launch failures, and stale recoveries were zero. Six workers remain requested and two effective under the existing swap safeguard: “swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two”. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.

Relative to the 09:58:46 scan start (12,201 validated), this scan began at 12,240, an observed short-window rate of approximately 836 receipts/hour. This is provisional; 17,698 receipts remained after scanning. Failure categories total 2,097 all-phenotype negative-variance and 242 no-reference-SNP. All 12 frozen 1% gates remain failed; no local-sharing inference is admissible.

This audit is read-only and does not modify settings, inputs, thresholds, or receipts. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_1001.json`.
