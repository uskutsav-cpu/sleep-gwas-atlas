# FI×sleep LAVA full receipt audit — 2026-09-26 17:38 UTC

- Audit interval: 17:37:03–17:38:26 UTC.
- Validated: 16,995/29,940 receipt files (56.76%). Five additional files appeared during the scan; end inventory and runner state both reported 17,000.
- Live jobs: PID 78548 on `longsleep` locus 1526 and PID 78549 on `snoring` locus 1525. Both process groups passed liveness checks.
- Integrity: no receipt, claim, identity, or duplicate issues; zero launch failures and zero stale recoveries.
- Concurrency: six requested, two effective under the persisted swap fallback.
- Failures: 3,010 all-phenotypes-negative-variance and 252 no-reference-SNP outcomes; all 12 trait gates fail the frozen 1% rule. No local-sharing inference is supported.
- Frozen hashes: lock `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared inputs `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Since the prior 17:30:14 state (16,975), the end state advanced 25 receipts in 8m12s (~183/hour); this short-window rate is not used as a stable completion estimate.

Raw audit: `lava_full_receipt_integrity_2026-09-26_1736_continuation.json`. The receipt set changed during scanning, so the audit correctly distinguishes 16,995 validated files from the 17,000 end-state count.
