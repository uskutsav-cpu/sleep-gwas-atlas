# Full LAVA receipt and locked failure-gate audit — 05:17 UTC

Receipt scan completed at 2026-09-26T05:16:47.527882+00:00. Runner state reported 9,866 at scan start/end while a receipt became durable; at 05:16:57.652443 UTC the coordinator caught up to the same 9,867-file count.

## Receipt integrity

All 9,867 receipt JSON files passed the locked `verify_receipt` validator, matched trait/pair/locus identity and frozen locus coordinates, and had a nonempty corresponding log. The receipt-file count reconciled with coordinator state immediately after its next update. Structural issues: **0**; duplicate trait/locus identities: **0**. Process statuses: 7,958 `PROCESSED`, 1,909 `PROCESS_FAILED`. Ordered receipt-stream SHA-256: `02e5b403bfe672af8097977e927388c149488ca4b44d4b52f731eb4292e51a95`.

At coordinator catch-up, two workers remained active on distinct claims: insomnia locus 989 (PID 78549) and sleep_efficiency locus 563 (PID 78548). Launch failures and stale recoveries were zero. The frozen analysis-lock hash remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest hash remains `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`. No workers were interrupted, loci rerun, settings changed, or analysis results overwritten.

## Locked 1% locus-failure gate

`PROCESS_FAILED` receipts alone exceed the frozen 1% limit for all 12 traits. This is a lower bound because the collator also counts non-tested FI/sleep univariate statuses.

| Sleep trait | Process failures | Lower-bound fraction | Negative-variance outcomes | No SNP IDs in reference | Passes 1% gate? |
|---|---:|---:|---:|---:|:---:|
| accel_sleep_duration | 587 | 23.53% | 566 | 21 | No |
| chronotype | 304 | 12.18% | 283 | 21 | No |
| insomnia | 80 | 3.21% | 60 | 20 | No |
| longsleep | 118 | 4.73% | 98 | 20 | No |
| napping | 86 | 3.45% | 66 | 20 | No |
| shortsleep | 105 | 4.21% | 85 | 20 | No |
| sleep_apnea | 46 | 1.84% | 26 | 20 | No |
| sleep_efficiency | 146 | 5.85% | 126 | 20 | No |
| sleep_timing | 146 | 5.85% | 126 | 20 | No |
| sleepdur | 85 | 3.41% | 65 | 20 | No |
| sleepiness | 111 | 4.45% | 91 | 20 | No |
| snoring | 95 | 3.81% | 75 | 20 | No |

Across audited receipts, 1,667 logs contain the recognized all-phenotype negative-variance outcome and 242 report no specified SNP IDs in the reference data. The latter map to loci 950–969 across all traits (chr6:25,684,630–33,864,262, spanning the MHC region) and locus 1484 for two traits (chr9:140,097,760–141,146,682). The audit does not establish why the reference intersection is empty; no outcome was reclassified or repaired.

Continue the locked campaign for receipt preservation, but local-sharing and downstream PLACO/fine-mapping/colocalization/molecular claims remain unjustified. Machine-readable detail is in `lava_receipt_full_audit_2026-09-26_0517.json`.
