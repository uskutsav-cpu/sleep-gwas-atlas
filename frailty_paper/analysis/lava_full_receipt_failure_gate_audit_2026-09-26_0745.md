# FI×sleep LAVA full receipt audit — 2026-09-26 07:45 UTC

Read-only scan: 2026-09-26T07:45:02+00:00 to 2026-09-26T07:45:06+00:00.

- Fixed snapshot validated 10,483/29,940 receipt files. File count 10,483 at start and 10,484 afterward; runner state reported 10,482 at both boundaries. The list changed during scanning, so one newer receipt was outside the validated snapshot.
- Receipt issues: 0; claim issues: 0; duplicate trait/locus receipts: 0; duplicate claims: 0.
- Active worker PIDs: 78548, 78549; auditor liveness checks: {'78548': True, '78549': True}. Requested/effective workers: 6/2; launch failures: 0; stale recoveries: 0.
- Process failures: 1,995 (1,753 negative variance; 242 no-reference-SNP). All 12 traits fail the locked 1% gate.

Analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Input-manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. Receipt stream SHA-256: `6e92766f76f1886cff8ce53fdb987afc7c38ec174994e7aa3f17ffee7a09ab03`. Auditor SHA-256: `45f12b71ba52d0b4112e414bfa4b3a724ee9b56a39bfb8bd13ded27d654e0307`.

The scan was read-only and leaves scientific settings, receipts, claims, and outputs unchanged. No local-sharing inference or downstream PLACO/fine-mapping/colocalization is supported. Machine-readable details: `lava_full_receipt_integrity_2026-09-26_0745.json`.
