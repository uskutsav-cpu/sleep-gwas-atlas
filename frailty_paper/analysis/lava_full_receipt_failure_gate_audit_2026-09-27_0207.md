# FI×sleep LAVA full-family failure-gate audit — 2026-09-27 02:08 UTC

Read-only full-family scan: 2026-09-27T02:07:54+00:00 to 2026-09-27T02:08:26+00:00 UTC.

- Valid receipts: 19,116/29,940; scan-start and post-scan inventories: 19,116/19,118; runner state: 19,117.
- Receipt / claim / duplicate identity / duplicate claim issues: 0/0/0/0.
- Active worker PIDs: [72534, 72535]; liveness: `{"72534": true, "72535": true}`; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0/0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3378, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 852/2,495 | 34.15% | FAIL |
| longsleep | 1,033/2,495 | 41.40% | FAIL |
| napping | 857/2,495 | 34.35% | FAIL |
| shortsleep | 947/2,495 | 37.96% | FAIL |
| sleep_apnea | 721/2,495 | 28.90% | FAIL |
| sleep_efficiency | 1,055/2,495 | 42.28% | FAIL |
| sleep_timing | 1,078/2,495 | 43.21% | FAIL |
| sleepdur | 884/2,495 | 35.43% | FAIL |
| sleepiness | 1,000/2,495 | 40.08% | FAIL |
| snoring | 954/2,495 | 38.24% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. The audit is read-only and did not change scientific inputs or thresholds.

The receipt inventory advanced during this scan.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-27_0207.json`.
