# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:46 UTC

Read-only full-family audit completed 2026-09-26T08:46:13+00:00–08:46:20+00:00 UTC.

- Expected receipts: 29,940; validated: 11,305; the inventory remained stable at 11,305 throughout the scan.
- Runner state advanced 11,304→11,306 during scanning. A post-audit coordinator read at 08:46:25 reached 11,307/29,940.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Active PIDs 78548 and 78549 passed auditor liveness checks. Requested/effective workers: 6/2 under the unchanged ≥90%-swap fallback. Launch failures: 0; stale claims recovered: 0.
- Process failures: 2,156 (1,914 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both unchanged; audit mutation flag `False`.

From 07:17:51 (10,227) to 08:46:25 (11,307), the provisional measured rate is ~732 receipts/hour over 88m34s; 18,633 remain (~25.5 hours at this rate). No run settings or thresholds changed. Failed gates continue to prohibit downstream local-sharing, PLACO, fine-mapping and colocalization claims.
