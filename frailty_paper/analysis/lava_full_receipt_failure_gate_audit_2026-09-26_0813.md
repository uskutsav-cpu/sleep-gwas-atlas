# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:13 UTC

Read-only full-family audit completed 2026-09-26T08:13:54+00:00–2026-09-26T08:13:59+00:00 UTC.

- Expected receipts: 29,940; validated: 10,895; file inventory grew 10,895→10,897 during the scan.
- Runner state was 10,894 at both audit boundaries. An immediate post-audit coordinator read at 08:13:59 reported 10,897/29,940; that snapshot matched the post-scan file inventory. Two active worker PIDs 78548, 78549 passed auditor liveness checks.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Requested/effective workers: 6/2. Existing fallback remains: swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two. Launch failures: 0; stale claims recovered: 0.
- Process failures: 2,057 (1,815 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both remain unchanged. Audit mutation flag: `False`.

From the 07:17:51 checkpoint (10,227) to the 08:13:59 runner snapshot (10,897), the measured rate is ~716/hour over 56m08s; 19,043 remain (~26.6 hours at this provisional rate). The audit did not alter the run, claims, receipts, inputs, or thresholds. No downstream local-sharing, PLACO, fine-mapping, or colocalization inference is permitted by the failed gates.
