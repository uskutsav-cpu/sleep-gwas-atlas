# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:36 UTC

Read-only full-family audit completed 2026-09-26T08:36:11+00:00–08:36:16+00:00 UTC.

- Expected receipts: 29,940; validated: 11,196; inventory grew 11,196→11,197 during scanning.
- Runner state advanced 11,194→11,197 across the audit. The 08:36:15 coordinator snapshot reconciled at 11,197/29,940.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Active worker PIDs 78548 and 78549 passed auditor liveness checks. Requested/effective workers: 6/2; swap safeguard remains at 2,772/3,072 MiB (≥90%). Launch failures: 0; stale claims recovered: 0.
- Process failures: 2,112 (1,870 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both remain unchanged. Audit mutation flag: `False`.

From 07:17:51 (10,227) to 08:36:15 (11,197), the provisional measured rate is ~742 receipts/hour over 78m24s; 18,743 remain (~25.3 hours at this rate). No run settings or thresholds changed. The failed gates do not permit downstream local-sharing, PLACO, fine-mapping, or colocalization inference.
