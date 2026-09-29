# FI×sleep LAVA full-family audit — 2026-09-26 10:50 UTC

The read-only full-family auditor scanned 2026-09-26T10:49:46+00:00–2026-09-26T10:49:52+00:00 UTC. It validated 12,898/29,940 receipt files. The filesystem contained 12,898 at scan start and 12,900 at scan end, while runner state held at 12,898→12,898. Receipt, claim, duplicate trait/locus, and duplicate claim issues were zero. Both active worker PIDs (78548, 78549) passed liveness checks on distinct claims: sleep_apnea locus 712 and sleep_apnea locus 1143. Launch failures and stale-claim recoveries were zero.

FI×sleep remains at 6 requested / 2 effective workers under its existing adaptive swap safeguard. 17,040/29,940 receipt slots remain. Process failure categories total 2,489: 2,247 all-phenotype negative-variance and 242 no-reference-SNP. All 12 frozen trait-level 1% gates fail, so local-sharing inference remains inadmissible.

Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and input-manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` are unchanged. The audit is read-only; scientific inputs, scheduler settings, receipts, and checkpoints were not changed. Machine-readable evidence: `lava_full_receipt_integrity_2026-09-26_1050.json`.
