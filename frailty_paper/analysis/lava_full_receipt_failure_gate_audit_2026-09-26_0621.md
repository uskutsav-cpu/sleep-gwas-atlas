# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:21 UTC

Read-only scan: 2026-09-26T06:20:59+00:00 to 2026-09-26T06:21:05+00:00.

- Validated 10,074/29,940 receipt files. Runner state advanced 10,072→10,074 during scanning; the receipt-file count matched at scan start/end.
- Receipt, claim, duplicate-identity, and duplicate-claim issue counts: 0, 0, 0, 0.
- Concurrency remains 6 requested/2 effective under the existing resource guard. Both active worker PIDs 78548, 78549 passed liveness checks; launch failures 0; stale recoveries 0.
- Process failures: 1,924 total (1,682 all-phenotype negative variance; 242 no specified SNPs in reference).
- Traits passing the locked 1% failure gate: 0/12.

Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

The audit is read-only; no receipts, thresholds, inputs, or output directories were changed. Every trait continues to fail the frozen family gate; no local-sharing estimate or downstream PLACO/fine-mapping/colocalization result is admissible from this family. Keep all failure receipts and continue the same resumable scheduler.
