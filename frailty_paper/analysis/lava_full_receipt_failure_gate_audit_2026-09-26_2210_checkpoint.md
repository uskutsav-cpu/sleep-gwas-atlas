# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 22:08

Read-only full-family scan: 2026-09-26 22:07:28 UTC to 2026-09-26 22:08:14 UTC.

- Valid receipts: 18,090 / 29,940; scan-start inventory 18,090, post-scan inventory 18,091, runner state 18,091.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3225, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 800/2,495 | 32.06% | FAIL |
| longsleep | 960/2,495 | 38.48% | FAIL |
| napping | 809/2,495 | 32.42% | FAIL |
| shortsleep | 890/2,495 | 35.67% | FAIL |
| sleep_apnea | 680/2,495 | 27.25% | FAIL |
| sleep_efficiency | 985/2,495 | 39.48% | FAIL |
| sleep_timing | 1,007/2,495 | 40.36% | FAIL |
| sleepdur | 826/2,495 | 33.11% | FAIL |
| sleepiness | 935/2,495 | 37.47% | FAIL |
| snoring | 897/2,495 | 35.95% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-26_2210_checkpoint.json`.
