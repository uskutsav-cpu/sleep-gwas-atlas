# Full-family FI×sleep LAVA receipt and gate audit — 05:40 UTC

Audit window: 2026-09-26T05:40:18+00:00 to 2026-09-26T05:40:23+00:00.

## Receipt integrity

The scan validated **9,945** receipts out of 29,940 expected. The file list was stable; runner state matched at 9,945 at scan start and 9,945 at scan end. All receipts matched the frozen trait/pair/locus identities and coordinates and had nonempty logs. Structural or schema issues: **0**; duplicate receipt identities: **0**; duplicate claim identities: **0**; malformed claims: **0**.

Active worker PIDs 78548, 78549 were process-verified alive. Requested/effective concurrency: 6/2; launch failures: 0; stale claim recoveries: 0. Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

The auditor reads the live runner and receipt/log files without mutating them. The ordered receipt hash stream is recorded in the machine-readable JSON.

## Frozen 1% locus-failure gate

A locus fails if process status is failed or either univariate FI/sleep status is not `TESTED`. The denominator remains the frozen 2,495 loci per trait. Since the campaign is incomplete, the observed failures are lower bounds; all 12 traits already exceed the locked 1% threshold.

| Sleep trait | Process failures | Observed failed loci (lower bound) | Fraction of 2,495 | Passes? |
|---|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 1786 | 71.58% | No |
| chronotype | 304 | 1326 | 53.15% | No |
| insomnia | 80 | 298 | 11.94% | No |
| longsleep | 118 | 362 | 14.51% | No |
| napping | 86 | 299 | 11.98% | No |
| shortsleep | 106 | 328 | 13.15% | No |
| sleep_apnea | 46 | 256 | 10.26% | No |
| sleep_efficiency | 147 | 383 | 15.35% | No |
| sleep_timing | 147 | 385 | 15.43% | No |
| sleepdur | 85 | 317 | 12.71% | No |
| sleepiness | 111 | 351 | 14.07% | No |
| snoring | 95 | 338 | 13.55% | No |

There are 1,670 all-phenotype negative-variance process failures and 242 no-reference-SNP process failures. The latter arise at the same 20 chr6 MHC-spanning loci across traits and one chr9 locus in two traits; the reason for the empty reference intersection has not been established. Unexpected process-failure classifications: 0.

No receipt was rerun or reclassified, and no threshold or scientific result changed. Preserve the scheduled run, but do not make local-sharing or downstream PLACO/fine-mapping/colocalization/molecular claims under the failed gate.

Machine-readable audit: `analysis/lava_full_receipt_integrity_2026-09-26_0540.json`. Auditor: `scripts/53_audit_lava_full_receipts.py`.
