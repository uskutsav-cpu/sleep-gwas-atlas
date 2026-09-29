# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:09 UTC

Read-only scan: 2026-09-26T06:09:01+00:00 to 2026-09-26T06:09:05+00:00.

- Validated receipts: 10,031/29,940; runner state advanced 10,030→10,031 during the scan.
- Receipt issues: 0; claim issues: 0; duplicate claims: 0; duplicate trait/locus identities: 0.
- Concurrency: 6 requested, 2 effective; swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two.
- Active claims: 2; PIDs 78548, 78549; launch failures this invocation: 0; stale claims recovered: 0.
- Process failures: 1,917 total (all phenotypes negative variance: 1,675, no specified snps in reference: 242).
- Locked maximum locus failure fraction: 1%; traits passing: 0/12.

Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This audit is read-only and validates receipt schemas, hashes, identities, matching logs, claims, and locked failure-gate results. It does not alter receipts, completed loci, thresholds, hashes, or output directories. The 1% gate fails for all traits, so this LAVA family does not support local-sharing inference or downstream PLACO/fine-mapping/colocalization. Continue the active atomic-claim scheduler and preserve all results.
