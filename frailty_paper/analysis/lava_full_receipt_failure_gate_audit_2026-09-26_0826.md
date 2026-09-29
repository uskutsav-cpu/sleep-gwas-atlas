# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:26 UTC

Read-only full-family audit completed 2026-09-26T08:26:30+00:00–08:26:36+00:00 UTC.

- Expected receipts: 29,940; validated: 11,065; inventory grew 11,065→11,067 during scanning.
- Runner state was 11,065 at both audit boundaries. An immediate 08:26:38 coordinator snapshot matched the post-scan inventory at 11,067/29,940.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Active worker PIDs 78548 and 78549 passed auditor liveness checks. Requested/effective workers: 6/2. Existing swap fallback remains at 2,772/3,072 MiB (≥90%); launch failures: 0; stale claims recovered: 0.
- Process failures: 2,074 (1,832 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both remain unchanged. Audit mutation flag: `False`.

A verified 60-second observation recorded coordinator progress from 11,041 to 11,063 receipts, with active assignments changing. From 07:17:51 (10,227) to 08:26:38 (11,067), the provisional rate is ~733 receipts/hour over 68m47s; 18,873 remain (~25.8 hours at this rate). The audit did not alter run state or thresholds. No downstream local-sharing, PLACO, fine-mapping, or colocalization inference is permitted under the failed gates.
