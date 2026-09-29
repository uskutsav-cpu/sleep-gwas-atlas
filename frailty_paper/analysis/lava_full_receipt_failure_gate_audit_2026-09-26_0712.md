# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 07:12 UTC

Read-only scan: 2026-09-26T07:12:55+00:00 to 2026-09-26T07:13:00+00:00.

- Validated 10,215/29,940 receipt files (34.12%). Runner state and file counts both equaled 10,215 at scan start/end.
- Receipt, claim, duplicate-identity, duplicate-claim, and matching-log issue counts: all zero.
- Requested/effective workers: 6/2. Active PIDs 78548, 78549 passed liveness checks; launch failures 0; stale recoveries 0. The adaptive guard reports swap use at 2,772/3,072 MiB (90%) and reduced the requested three workers to two.
- Process failures: 1,954 total (1,712 negative-variance; 242 no-reference-SNP). Traits passing locked 1% gate: 0/12.

Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This was a read-only audit. Existing receipts, completed loci, thresholds, hashes, and output directories remain untouched. Every trait remains outside the frozen family QC gate; no local-sharing inference or downstream PLACO/fine-mapping/colocalization is admissible from this family.
