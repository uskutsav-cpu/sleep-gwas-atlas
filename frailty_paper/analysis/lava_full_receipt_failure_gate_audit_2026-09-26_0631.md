# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:31 UTC

Read-only scan: 2026-09-26T06:31:15+00:00 to 2026-09-26T06:31:19+00:00.

- Validated 10,098/29,940 receipts. Runner state was 10,098→10,098; receipt-file counts matched at both scan boundaries.
- Receipt/claim integrity errors: 0/0; duplicate receipt identities/claims: 0/0.
- Concurrency: 6 requested, 2 effective. PIDs 78548, 78549 passed liveness checks. Launch failures: 0; stale recoveries: 0.
- Process failures: 1,927 total (1,685 all-phenotype negative variance; 242 no-reference-SNP).
- Traits passing the frozen 1% gate: 0/12.

Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

Read-only integrity check only: receipts, thresholds, input hashes, and output directories are preserved. The FI×sleep local-correlation family remains incomplete and fails its frozen gate for all traits. No local-sharing inference or downstream PLACO/fine-mapping/colocalization is authorized by this family; preserve failures and continue the existing scheduler.
