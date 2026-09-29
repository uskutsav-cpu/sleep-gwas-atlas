# FI×sleep LAVA full receipt audit — 2026-09-26 07:38 UTC

Read-only audit interval: 2026-09-26T07:38:40+00:00 to 2026-09-26T07:38:47+00:00.

- Fixed snapshot: 10,399/29,940 receipt files validated. File count changed during scan: True (10,399 at start; 10,400 after). Runner state reported 10,398 at scan start and 10,398 at scan end, so the file inventory was ahead by one and two receipts respectively.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipt identities: 0; duplicate claims: 0.
- Active claims: 2; PIDs 78548, 78549, each process-verified live by the auditor. Requested/effective workers: 6/2; launch failures: 0; stale recoveries: 0.
- Process failures: 1,984 (1,742 all-phenotype negative variance; 242 no specified SNPs in reference). All 12 traits fail the locked 1% gate.

Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Receipt hash-stream SHA-256: `866cad54cbd61c88b3ba40e5483209b022d3abb8cd421614048205a48eab63bd`. Auditor SHA-256: `45f12b71ba52d0b4112e414bfa4b3a724ee9b56a39bfb8bd13ded27d654e0307`.

This is a read-only fixed-snapshot audit. Scientific settings, receipts, claims, and outputs were not changed. No local-sharing inference or downstream PLACO/fine-mapping/colocalization is supported. Machine-readable details: `lava_full_receipt_integrity_2026-09-26_0738.json`.
