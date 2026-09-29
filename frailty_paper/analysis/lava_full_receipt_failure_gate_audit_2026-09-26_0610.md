# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:10 UTC

Read-only scan: 2026-09-26T06:09:55+00:00 to 2026-09-26T06:09:59+00:00.

- Validated receipts: 10,035/29,940; runner state reported 10,034 at scan start and 10,034 at scan end. File count matched across the scan.
- Receipt/claim issues: 0 receipt issues, 0 claim issues, 0 duplicate receipt identities, 0 duplicate claims.
- Concurrency: 6 requested, 2 effective under the existing adaptive swap guard.
- Active claims: 2; PIDs 78548, 78549; both process-liveness probes passed. Launch failures: 0; stale claims recovered: 0.
- Process failures: 1,917 total (1,675 all-phenotype negative variance; 242 no specified SNPs in reference).
- Frozen failure threshold: 1%; traits passing: 0/12.

Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This audit is read-only and preserves the existing receipts, completed loci, thresholds, hashes, and output directories. Every trait still fails the locked family gate; no FI local-sharing inference or downstream PLACO/fine-mapping/colocalization is admissible from this family. Continue the current atomic-claim scheduler and retain all failures.
