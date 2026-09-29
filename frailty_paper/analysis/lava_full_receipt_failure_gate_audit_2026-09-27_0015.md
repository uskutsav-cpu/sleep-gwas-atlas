# FI×sleep LAVA full-family failure-gate audit — 2026-09-27 00:15 UTC

Read-only full-family scan: 2026-09-27T00:15:54+00:00 to 2026-09-27T00:16:26+00:00 UTC.

- Valid receipts: 18,420 / 29,940; scan-start inventory 18,420, post-scan inventory 18,422, runner state 18,422.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3283, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 818/2,495 | 32.79% | FAIL |
| longsleep | 980/2,495 | 39.28% | FAIL |
| napping | 823/2,495 | 32.99% | FAIL |
| shortsleep | 907/2,495 | 36.35% | FAIL |
| sleep_apnea | 694/2,495 | 27.82% | FAIL |
| sleep_efficiency | 1,006/2,495 | 40.32% | FAIL |
| sleep_timing | 1,034/2,495 | 41.44% | FAIL |
| sleepdur | 841/2,495 | 33.71% | FAIL |
| sleepiness | 952/2,495 | 38.16% | FAIL |
| snoring | 914/2,495 | 36.63% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. The audit is read-only and does not change scientific inputs or thresholds.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-27_0015.json`.
