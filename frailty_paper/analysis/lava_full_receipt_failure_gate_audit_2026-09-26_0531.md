# Full-family FI×sleep LAVA receipt and gate audit — 05:31 UTC

Audit window: 2026-09-26T05:31:36+00:00 to 2026-09-26T05:31:40+00:00.

## Receipt integrity

The scan validated **9,913** receipts out of 29,940 expected, with a stable file list and runner count 9,913 at start and 9,913 at end. All matched the locked trait/pair/locus identity and frozen coordinates and had nonempty logs. Structural/identity/log issues: **0**; duplicate receipt identities: **0**; duplicate active claim identities: **0**; claim-file issues: **0**.

Two active worker PIDs were process-verified alive: 78548, 78549. Requested/effective concurrency is 6/2; launch failures 0; stale claim recoveries 0. Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. Prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.

The auditor snapshots only local metadata and reads receipt/log files; it does not mutate the runner or receipts. Receipt-file stream SHA-256 is recorded in the machine-readable JSON.

## Frozen 1% locus-failure gate

The collator flags a locus if process status failed or either univariate FI/sleep status is not `TESTED`. The denominator remains the locked 2,495 loci per trait; because the family is incomplete, these observed failures are lower bounds over the planned family. All 12 traits already exceed the 1% limit.

| Sleep trait | Process failures | Observed failed loci (lower bound) | Fraction of 2,495 | Passes? |
|---|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 1786 | 71.58% | No |
| chronotype | 304 | 1326 | 53.15% | No |
| insomnia | 80 | 296 | 11.86% | No |
| longsleep | 118 | 360 | 14.43% | No |
| napping | 86 | 298 | 11.94% | No |
| shortsleep | 105 | 326 | 13.07% | No |
| sleep_apnea | 46 | 254 | 10.18% | No |
| sleep_efficiency | 146 | 382 | 15.31% | No |
| sleep_timing | 147 | 382 | 15.31% | No |
| sleepdur | 85 | 316 | 12.67% | No |
| sleepiness | 111 | 349 | 13.99% | No |
| snoring | 95 | 336 | 13.47% | No |

Across the validated receipts, 1,668 process failures report negative variance for all phenotypes and 242 report that no specified SNP IDs occur in the reference. The no-SNP outcomes map to loci 950–969 across traits (chr6 MHC-spanning interval) and locus 1484 for two traits (chr9); the cause of the empty reference intersection is not established. Unexpected process-failure classifications: 0.

No QC threshold, receipt, or scientific result was changed. Continue the locked campaign for receipt preservation, but do not run local-sharing inference, PLACO, fine-mapping, colocalization, or molecular interpretation under the failed gate.

Machine-readable snapshot: `analysis/lava_full_receipt_integrity_2026-09-26_0531.json`. Reusable auditor: `scripts/53_audit_lava_full_receipts.py`.
