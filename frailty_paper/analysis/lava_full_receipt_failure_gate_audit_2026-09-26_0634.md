# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:34 UTC

Read-only scan: 2026-09-26T06:33:51+00:00 to 2026-09-26T06:33:55+00:00.

- Validated 10,105/29,940 receipt files; filesystem count remained 10,105 throughout.
- Runner state reported 10,104 at scan start and 10,104 at scan end, one receipt behind the stable filesystem count. No identity, claim, duplicate, receipt, or log-integrity issues were recorded.
- Concurrency: 6 requested, 2 effective under the adaptive swap guard. Two active PIDs 78548, 78549 passed liveness checks; launch failures 0; stale recoveries 0.
- Process failures: 1,927 (1,685 all-phenotype negative variance; 242 no specified SNPs in reference). All 12 traits fail the locked 1% gate.

Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

Read-only audit only; do not replace or mutate any receipts, thresholds, hashes, or output directories. The FI×sleep LAVA family remains incomplete and fails its frozen gate. No local-sharing result or downstream PLACO/fine-mapping/colocalization is admissible from this family.
