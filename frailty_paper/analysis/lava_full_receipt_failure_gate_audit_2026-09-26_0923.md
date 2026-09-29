# FI×sleep LAVA full-family audit — 2026-09-26 09:23 UTC

The read-only full-family auditor validated 11,767/29,940 immutable receipt files from 09:23:57 to 09:24:03 UTC. The inventory grew from 11,767 to 11,768 during the scan. Remaining after the scan: 18,172 receipts.

Receipt integrity, duplicate identities, claim issues, and duplicate claims were clean. Active worker PIDs 78548, 78549 passed liveness checks. Launch failures and stale recoveries were zero. Six workers remain requested and 2 effective under the existing resource safeguard. Lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; input-manifest SHA-256 remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

From the 09:22 audit (11,745 validated at 09:22:16 UTC) to this audit (11,767 validated at 09:23:57 UTC), the short-interval rate was approximately 784 receipts/hour. This is noisy and provisional. There were 2,010 all-phenotype negative-variance and 242 no-reference-SNP process outcomes. All 12 frozen 1% failure gates remain failed, so no local-sharing result is admissible.

No run settings, analysis inputs, thresholds, or receipts were modified by the audit. Machine-readable audit: `lava_full_receipt_integrity_2026-09-26_0923.json`.
