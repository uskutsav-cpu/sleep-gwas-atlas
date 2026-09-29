# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:34 UTC

Read-only full-family audit completed 2026-09-26T08:34:09+00:00–2026-09-26T08:34:13+00:00 UTC.

- Expected receipts: 29,940; validated: 11,167; inventory grew from 11,167 at scan start to 11,168 at scan end.
- Runner count advanced 11,166→11,166 across the scan; an immediate later coordinator snapshot at 08:34:24 reported 11,171/29,940. Both active worker PIDs 78548, 78549 passed auditor liveness checks.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Requested/effective workers: 6/2. Existing swap fallback remains active. Launch failures: 0; stale claims recovered: 0.
- Process failures: 2,109 (1,867 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both unchanged. Audit mutation flag: `False`.

From the 07:17:51 checkpoint (10,227) to the 08:34:24 snapshot (11,171), measured throughput is ~740/hour over 76m33s; 18,769 receipts remain (~25.4 hours at this provisional rate). The audit did not change receipts, claims, inputs, or thresholds. No downstream local-sharing, PLACO, fine-mapping, or colocalization inference is permitted by the failed gates.
