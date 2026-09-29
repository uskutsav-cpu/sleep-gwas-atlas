# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:23 UTC

Read-only full-family audit completed 2026-09-26T08:23:40+00:00–08:23:45+00:00 UTC.

- Expected receipts: 29,940; validated: 11,025; the file inventory grew 11,025→11,026 during scanning.
- Runner state was 11,023 at both scan boundaries. The immediate 08:23:46 coordinator snapshot reconciled with the post-scan inventory at 11,026/29,940.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Active worker PIDs 78548 and 78549 passed auditor liveness checks. Requested/effective workers: 6/2. Existing swap fallback remains at 2,772/3,072 MiB (≥90%); launch failures: 0; stale claims recovered: 0.
- Process failures: 2,067 (1,825 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both remain unchanged. Audit mutation flag: `False`.

From 07:17:51 (10,227) to 08:23:46 (11,026), the provisional measured rate is ~727 receipts/hour over 65m55s; 18,914 receipts remain (~26.0 hours at this rate). The audit did not alter the run, claims, receipts, inputs, or thresholds. No downstream local-sharing, PLACO, fine-mapping, or colocalization inference is permitted by the failed gates.
