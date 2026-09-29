# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:04 UTC

Read-only scan: 2026-09-26T06:04:34+00:00 to 2026-09-26T06:04:39+00:00.

- Validated receipts: 10,016/29,940; runner state advanced 10,015→10,016 during the scan.
- Receipt, identity, claim, and duplicate issues: 0, 0, 0, 0.
- Concurrency: 6 requested, 2 effective; the existing swap guard reports swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two.
- Active claims: 2; PIDs 78548, 78549; launch failures this invocation: 0; stale claims recovered: 0.
- Process failures: 1,915 total (1,673 all-phenotype negative variance; 242 no specified SNPs in reference).
- Locked maximum locus failure fraction: 1%; traits passing: 0/12.

Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This audit is read-only. It does not change receipts, completed loci, thresholds, hashes, or output directories. The frozen 1% failure gate fails for every trait; no local-sharing inference or downstream PLACO/fine-mapping/colocalization is admissible from this family. Continue the current atomic-claim scheduler and retain all failures.
