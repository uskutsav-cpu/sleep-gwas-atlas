# FI×sleep LAVA full-family failure-gate audit — 2026-09-27 00:46 UTC

Read-only full-family scan: 2026-09-27T00:45:59+00:00 to 2026-09-27T00:46:24+00:00 UTC.

- Valid receipts: 18,575 / 29,940; scan-start inventory 18,575, post-scan inventory 18,577, runner state 18,577.
- Receipt / claim / duplicate trait-locus / duplicate claim issues: 0 / 0 / 0 / 0.
- Active worker PIDs: [72534, 72535]; liveness: `{"72534": true, "72535": true}`; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3302, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 824/2,495 | 33.03% | FAIL |
| longsleep | 991/2,495 | 39.72% | FAIL |
| napping | 830/2,495 | 33.27% | FAIL |
| shortsleep | 914/2,495 | 36.63% | FAIL |
| sleep_apnea | 699/2,495 | 28.02% | FAIL |
| sleep_efficiency | 1,014/2,495 | 40.64% | FAIL |
| sleep_timing | 1,043/2,495 | 41.80% | FAIL |
| sleepdur | 847/2,495 | 33.95% | FAIL |
| sleepiness | 962/2,495 | 38.56% | FAIL |
| snoring | 921/2,495 | 36.91% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. The audit is read-only and did not change scientific inputs or thresholds.

The receipt file list changed during scanning: two additional valid receipts arrived while the audit ran. The audit validated the 18,575 files present at scan start; the runner reported 18,577 by scan end. No issues were found in the validated set.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-27_0045.json`.
