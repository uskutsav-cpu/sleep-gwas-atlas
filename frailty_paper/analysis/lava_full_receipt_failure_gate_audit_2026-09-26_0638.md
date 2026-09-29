# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:38 UTC

Read-only scan: 2026-09-26T06:37:52+00:00 to 2026-09-26T06:38:04+00:00.

- Validated 10,115/29,940 receipt files; runner state matched at 10,115 at scan start/end and receipt-file count was stable.
- Receipt, claim, duplicate-identity, duplicate-claim, and matching-log issues: all zero.
- Requested/effective concurrency: 6/2. Active PIDs 78548, 78549 passed liveness checks; launch failures 0; stale recoveries 0.
- Process failures: 1,929 (1,687 all-phenotype negative variance; 242 no-reference-SNP); traits passing frozen 1% gate: 0/12.

Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This was a read-only integrity audit. Existing receipts, checkpoints, hashes, settings, and outputs remain untouched. All 12 traits still fail the frozen family QC; do not report local-sharing results or proceed to downstream PLACO/fine-mapping/colocalization from this family.
