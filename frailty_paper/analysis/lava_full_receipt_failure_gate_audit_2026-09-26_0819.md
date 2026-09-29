# FI×sleep LAVA full receipt and frozen failure-gate audit — 2026-09-26 08:19 UTC

Read-only full-family audit completed 2026-09-26T08:19:37+00:00–2026-09-26T08:19:42+00:00 UTC.

- Expected receipts: 29,940; validated: 10,971; inventory remained stable at 10,971 files during the scan.
- Runner state was 10,969 at both scan boundaries. The subsequent 08:19:43 coordinator snapshot advanced to 10,972/29,940.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claim identities: 0.
- Active worker PIDs 78548 and 78549 passed liveness checks. Requested/effective workers: 6/2. The existing fallback remains because swap reached 2,772/3,072 MiB (≥90%); launch failures: 0; stale claims recovered: 0.
- Process failures: 2,060 (1,818 all-phenotype negative-variance failures; 242 no-reference-SNP failures). All 12 trait-level frozen 1% gates fail.
- Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Both remain unchanged. Audit mutation flag: `False`.

From 07:17:51 (10,227) to the 08:19:43 coordinator snapshot (10,972), the provisional measured rate is approximately 723 receipts/hour over 61m52s; 18,968 receipts remain (~26.3 hours at this short-window rate). The audit did not alter the run, claims, receipts, inputs, or thresholds. No downstream local-sharing, PLACO, fine-mapping, or colocalization inference is permitted by the failed gates.
