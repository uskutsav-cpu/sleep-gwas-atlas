# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 21:20 UTC

Read-only full-family scan: 2026-09-26T21:19:40+00:00 to 2026-09-26T21:20:14+00:00.

- Valid receipts: 17,838 / 29,940; scan-start inventory 17,838, post-scan inventory 17,842, runner state 17,842.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3178, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 786/2,495 | 31.50% | FAIL |
| longsleep | 939/2,495 | 37.64% | FAIL |
| napping | 792/2,495 | 31.74% | FAIL |
| shortsleep | 873/2,495 | 34.99% | FAIL |
| sleep_apnea | 666/2,495 | 26.69% | FAIL |
| sleep_efficiency | 968/2,495 | 38.80% | FAIL |
| sleep_timing | 987/2,495 | 39.56% | FAIL |
| sleepdur | 807/2,495 | 32.34% | FAIL |
| sleepiness | 916/2,495 | 36.71% | FAIL |
| snoring | 880/2,495 | 35.27% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-26_2119_checkpoint.json`.
