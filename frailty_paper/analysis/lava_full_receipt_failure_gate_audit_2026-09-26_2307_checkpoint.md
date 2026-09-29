# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 23:07

Read-only full-family scan: 2026-09-26T23:07:06+00:00 to 2026-09-26T23:07:58+00:00 UTC.

- Valid receipts: 18,234 / 29,940; scan-start inventory 18,234, post-scan inventory 18,237, runner state 18,237.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3255, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 809/2,495 | 32.42% | FAIL |
| longsleep | 969/2,495 | 38.84% | FAIL |
| napping | 816/2,495 | 32.71% | FAIL |
| shortsleep | 900/2,495 | 36.07% | FAIL |
| sleep_apnea | 687/2,495 | 27.54% | FAIL |
| sleep_efficiency | 995/2,495 | 39.88% | FAIL |
| sleep_timing | 1,018/2,495 | 40.80% | FAIL |
| sleepdur | 834/2,495 | 33.43% | FAIL |
| sleepiness | 946/2,495 | 37.92% | FAIL |
| snoring | 906/2,495 | 36.31% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-26_2307_checkpoint.json`.
