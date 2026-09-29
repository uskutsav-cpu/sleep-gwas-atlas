# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:13 UTC

Read-only scan: 2026-09-26T06:12:40+00:00 to 2026-09-26T06:12:46+00:00.

- Validated receipt files: 10,045/29,940; runner state reported 10,045 at scan start and 10,045 at scan end. File count remained stable through the scan.
- Integrity: 0 receipt issues, 0 claim issues, 0 duplicate receipt identities, and 0 duplicate claims.
- Concurrency: 6 requested, 2 effective. Two live claims are held by PIDs 78548, 78549; both passed process-liveness checks. Launch failures: 0; stale claims recovered: 0.
- Process failures: 1,917 (1,675 all-phenotype negative variance; 242 no specified SNPs in reference).
- Locked failure threshold: 1%; traits passing: 0/12.

Frozen analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

The audit is read-only. It preserves all receipts, completed loci, thresholds, hashes, and output directories. All 12 traits still fail the frozen gate, so no local-sharing result or downstream PLACO/fine-mapping/colocalization is admissible from this family. Continue the current atomic-claim scheduler and preserve all failures.
