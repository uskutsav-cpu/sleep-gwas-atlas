# FI×sleep LAVA full-family audit — 2026-09-26 09:39 UTC

The full read-only auditor validated 11,958/29,940 receipts from 09:39:32 to 09:39:37 UTC. Inventory advanced 11,958→11,960; runner state advanced 11,957→11,959.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims. Receipt/claim/identity/duplicate issues, launch failures, and stale recoveries were zero. Six workers remain requested and two effective under the existing swap safeguard: “swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two”. Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.

Relative to the 09:36 audit (11,924 validated at 09:36:58), the 09:39 scan began at 11,958, an observed short-window rate of approximately 795 receipts/hour. This is provisional; 17,980 receipts remain after scanning. Failure categories total 2,034 all-phenotype negative-variance and 242 no-reference-SNP. All 12 frozen 1% gates remain failed; no local-sharing inference is admissible.

This audit is read-only and does not modify settings, inputs, thresholds, or receipts. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_0939.json`.
