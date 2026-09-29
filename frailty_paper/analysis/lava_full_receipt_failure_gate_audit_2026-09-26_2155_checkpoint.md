# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 21:52

Read-only full-family scan: 2026-09-26 21:51:52 UTC to 2026-09-26 21:52:39 UTC.

- Valid receipts: 18,027 / 29,940; scan-start inventory 18,027, post-scan inventory 18,029, runner state 18,029.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3212, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 797/2,495 | 31.94% | FAIL |
| longsleep | 954/2,495 | 38.24% | FAIL |
| napping | 805/2,495 | 32.26% | FAIL |
| shortsleep | 885/2,495 | 35.47% | FAIL |
| sleep_apnea | 677/2,495 | 27.13% | FAIL |
| sleep_efficiency | 981/2,495 | 39.32% | FAIL |
| sleep_timing | 1,001/2,495 | 40.12% | FAIL |
| sleepdur | 819/2,495 | 32.83% | FAIL |
| sleepiness | 930/2,495 | 37.27% | FAIL |
| snoring | 891/2,495 | 35.71% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-26_2155_checkpoint.json`.
