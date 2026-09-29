# FI×sleep LAVA full-family audit — 2026-09-26 10:48 UTC

The read-only full-family auditor scanned 2026-09-26T10:47:18+00:00 to 2026-09-26T10:47:24+00:00. It validated 12,867/29,940 receipt files. The inventory grew from 12,867 to 12,869 during scanning; runner state advanced 12,866→12,869. Receipt, claim, duplicate trait/locus, and duplicate claim issues were zero. Both active worker PIDs (78548, 78549) passed liveness checks and held distinct claims: insomnia locus 1142 and sleepdur locus 710. Launch failures and stale-claim recoveries were zero.

FI×sleep remains at 6 requested / 2 effective workers under its recorded adaptive swap safeguard. 17,071 receipts remain (57.02%). Process failure categories total 2,485: 2,243 all-phenotype negative-variance and 242 no-reference-SNP. All 12 frozen trait-level 1% gates fail, so local-sharing inference remains inadmissible.

Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged. The audit is read-only and did not alter scientific inputs, scheduler settings, receipts, or checkpoints. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1048.json`.
