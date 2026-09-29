# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 21:25 UTC

Read-only full-family scan: 2026-09-26T21:24:40+00:00 to 2026-09-26T21:25:11+00:00.

- Valid receipts: 17,887 / 29,940; scan-start inventory 17,887, post-scan inventory 17,892, runner state 17,891.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3180, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 788/2,495 | 31.58% | FAIL |
| longsleep | 942/2,495 | 37.76% | FAIL |
| napping | 795/2,495 | 31.86% | FAIL |
| shortsleep | 875/2,495 | 35.07% | FAIL |
| sleep_apnea | 668/2,495 | 26.77% | FAIL |
| sleep_efficiency | 972/2,495 | 38.96% | FAIL |
| sleep_timing | 990/2,495 | 39.68% | FAIL |
| sleepdur | 810/2,495 | 32.46% | FAIL |
| sleepiness | 921/2,495 | 36.91% | FAIL |
| snoring | 881/2,495 | 35.31% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-26_2125_checkpoint.json`.
