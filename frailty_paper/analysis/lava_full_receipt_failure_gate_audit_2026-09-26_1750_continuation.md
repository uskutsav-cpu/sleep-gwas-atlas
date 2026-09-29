# FI×sleep LAVA full receipt audit — 2026-09-26 17:52 UTC

- Audit interval: 17:50:43–17:51:59 UTC.
- Validated: 17,040/29,940 receipt files (56.95%). Two additional files appeared during the scan; the end inventory had 17,042 files and runner state reported 17,041.
- Live jobs: PID 78549 on `napping` locus 1530 and PID 78548 on `shortsleep` locus 1530. Both process groups passed liveness checks.
- Integrity: no receipt, claim, identity, or duplicate issues; zero launch failures and zero stale recoveries.
- Concurrency: six requested, two effective under the persisted swap fallback.
- Failures: 3,023 all-phenotypes-negative-variance and 252 no-reference-SNP outcomes; all 12 trait gates fail the frozen 1% rule. No local-sharing inference is supported.
- Frozen hashes: lock `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared inputs `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Since the prior 17:38:26 end-state (17,000), runner state advanced 41 receipts in 13m33s (~181/hour); this short-window rate is not used as a stable completion estimate.

Raw audit: `lava_full_receipt_integrity_2026-09-26_1750_continuation.json`. The receipt set changed during scanning, so the audit distinguishes 17,040 validated files from the end inventory of 17,042.
