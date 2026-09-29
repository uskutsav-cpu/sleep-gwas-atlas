# FI×sleep LAVA full receipt and failure-gate audit — 2026-09-26 07:28 UTC

Read-only scan: 2026-09-26T07:28:22+00:00 to 2026-09-26T07:28:32+00:00.

- Fixed file snapshot: 10,259 receipts validated. The file count grew 10,259→10,261 during the scan, and runner state grew 10,259→10,261; the two later files are outside this audit snapshot.
- Receipt issues: 0; claim issues: 0; duplicate receipt identities: 0; duplicate claims: 0.
- Requested/effective workers: 6/2; active PIDs 78548, 78549; launch failures: 0; stale recoveries: 0.
- Process failures: 1,956 (1,714 all-phenotype negative variance; 242 no-reference-SNP). All 12 traits fail the locked 1% gate.

Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

This read-only audit validates the fixed receipt-file snapshot and does not modify receipts, claims, settings, or outputs. All traits remain outside the frozen QC gate; no local-sharing inference or downstream PLACO/fine-mapping/colocalization is supported.
