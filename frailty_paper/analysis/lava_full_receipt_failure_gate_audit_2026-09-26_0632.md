# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:32 UTC

Read-only scan: 2026-09-26T06:32:20+00:00 to 2026-09-26T06:32:25+00:00.

- Validated receipts: 10,101/29,940; runner state was 10,100 at both boundaries.
- Receipt, identity, claim, duplicate and matching-log issues: 0, 0, 0, 0, and 0.
- Concurrency: 6 requested, 2 effective; swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two.
- Active claims: 2; shortsleep locus 1001 (worker PID 78549); sleep_timing locus 574 (worker PID 78548); launch failures: 0; stale claims recovered: 0.
- Process failures: 1,927 total (1,685 all-phenotype negative variance; 242 no-reference-SNP).
- Locked maximum locus failure fraction: 1%; traits passing: 0/12.

Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This read-only audit validates receipt schemas, hashes, identities, matching logs, claims, and frozen failure-gate results. It changed no receipts, completed loci, thresholds, hashes, or output directories. All 12 traits fail the frozen gate, so local-sharing inference and downstream PLACO/fine-mapping/colocalization are not justified. Continue the active atomic-claim scheduler and preserve results.
