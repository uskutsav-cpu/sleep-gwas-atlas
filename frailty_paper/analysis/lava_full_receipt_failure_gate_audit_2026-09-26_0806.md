# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:06 UTC

Read-only full-family audit completed 2026-09-26T08:06:50+00:00–2026-09-26T08:06:54+00:00 UTC.

- Expected receipts: 29,940; validated: 10,803; inventory grew 10,803→10,804 while scanning.
- Runner count advanced 10,801→10,804; latest runner snapshot: 10,804/29,940. Both active worker PIDs 78548, 78549 passed auditor liveness checks.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Requested/effective workers: 6/2. Existing fallback reason remains: swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two. Launch failures: 0; stale claims recovered: 0.
- Process failures: 2,037 (1,795 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both remain unchanged. Audit mutation flag: `False`.

The file inventory grew by one while the audit ran, and the coordinator count advanced by three; the post-scan inventory and coordinator both report 10,804. This audit did not alter the run, claims, receipts, input, or thresholds. No downstream local-sharing, PLACO, fine-mapping, or colocalization inference is permitted by the failed gates.
