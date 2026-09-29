# FI×sleep LAVA full-family audit — 2026-09-26 09:36 UTC

The full read-only auditor validated 11,924/29,940 receipt files from 09:36:58 to 09:37:03 UTC. Inventory advanced 11,924→11,925; runner state was 11,923 at both read boundaries.

Both worker PIDs (78548, 78549) passed liveness checks on distinct claims. Receipt integrity, duplicate identities, claim issues and duplicate claims were clean; launch failures and stale recoveries were zero. The existing scheduler remains at 6 requested/2 effective. Its persisted fallback is “swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two”. Analysis-lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged.

Between the 09:33:16 audit (11,875 validated) and this audit start (11,924 validated at 09:36:58), the short-window rate was approximately 795 receipts/hour; this estimate is noisy and provisional. 18,015 receipts remain after the scan. The failure categories total 2,033 all-phenotype negative-variance and 242 no-reference-SNP outcomes. All 12 frozen 1% gates still fail, so no local-sharing inference is admissible.

This audit is read-only; it changes no input, setting, threshold, or receipt. Machine-readable details: `lava_full_receipt_integrity_2026-09-26_0936.json`.
