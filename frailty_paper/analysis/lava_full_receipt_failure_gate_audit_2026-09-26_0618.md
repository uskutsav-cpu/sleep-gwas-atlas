# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 06:18 UTC

Read-only scan: 2026-09-26T06:17:58+00:00 to 2026-09-26T06:18:02+00:00.

- Validated receipts: 10,064/29,940; runner state matched at 10,064 before and after the scan, and the receipt-file count remained stable.
- Integrity: no receipt/claim issues, duplicate receipt identities, or duplicate claims.
- Concurrency: 6 requested, 2 effective under the existing adaptive swap guard. Active PIDs 78548, 78549 both passed liveness checks; launch failures 0, stale recoveries 0.
- Process failures: 1,919 total (1,677 negative-variance, 242 no-reference-SNP).
- Traits passing locked 1% failure gate: 0/12.

Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This scan is read-only; receipts, completed loci, thresholds, hashes, and output directories remain unchanged. All traits still fail the frozen family gate, so no local-sharing or downstream locus inference is admissible from this family. Preserve failure receipts and continue the existing atomic-claim scheduler.
