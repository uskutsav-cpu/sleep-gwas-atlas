# FI×sleep LAVA full-family failure-gate audit — 2026-09-27 01:35 UTC

Read-only full-family scan: 2026-09-27T01:35:25+00:00 to 2026-09-27T01:35:43+00:00 UTC.

- Valid receipts: 18,899/29,940; scan-start and post-scan inventories: 18,899/18,899; runner state: 18,899.
- Receipt / claim / duplicate identity / duplicate claim issues: 0/0/0/0.
- Active worker PIDs: [72534, 72535]; liveness: `{"72534": true, "72535": true}`; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0/0.
- Frozen lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3350, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 841/2,495 | 33.71% | FAIL |
| longsleep | 1,019/2,495 | 40.84% | FAIL |
| napping | 846/2,495 | 33.91% | FAIL |
| shortsleep | 936/2,495 | 37.52% | FAIL |
| sleep_apnea | 713/2,495 | 28.58% | FAIL |
| sleep_efficiency | 1,039/2,495 | 41.64% | FAIL |
| sleep_timing | 1,070/2,495 | 42.89% | FAIL |
| sleepdur | 871/2,495 | 34.91% | FAIL |
| sleepiness | 983/2,495 | 39.40% | FAIL |
| snoring | 943/2,495 | 37.80% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. The audit is read-only and did not change scientific inputs or thresholds.

The receipt inventory was stable during this scan.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-27_0135.json`.
