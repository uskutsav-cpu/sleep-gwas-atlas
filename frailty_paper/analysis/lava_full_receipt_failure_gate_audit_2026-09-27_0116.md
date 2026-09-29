# FI×sleep LAVA full-family failure-gate audit — 2026-09-27 01:16 UTC

Read-only full-family scan: 2026-09-27T01:16:19+00:00 to 2026-09-27T01:16:29+00:00 UTC.

- Valid receipts: 18,735 / 29,940; scan-start inventory 18,735, post-scan inventory 18,737, runner state 18,736.
- Receipt / claim / duplicate trait-locus / duplicate claim issues: 0 / 0 / 0 / 0.
- Active worker PIDs: [72534, 72535]; liveness: `{"72534": true, "72535": true}`; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3316, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 829/2,495 | 33.23% | FAIL |
| longsleep | 1,006/2,495 | 40.32% | FAIL |
| napping | 837/2,495 | 33.55% | FAIL |
| shortsleep | 924/2,495 | 37.03% | FAIL |
| sleep_apnea | 705/2,495 | 28.26% | FAIL |
| sleep_efficiency | 1,025/2,495 | 41.08% | FAIL |
| sleep_timing | 1,055/2,495 | 42.28% | FAIL |
| sleepdur | 858/2,495 | 34.39% | FAIL |
| sleepiness | 970/2,495 | 38.88% | FAIL |
| snoring | 933/2,495 | 37.39% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. The audit is read-only and did not change scientific inputs or thresholds.

The receipt file list changed during scanning: 2 additional receipt files appeared. The audit validated the 18,735 files present at scan start.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-27_0114.json`.
