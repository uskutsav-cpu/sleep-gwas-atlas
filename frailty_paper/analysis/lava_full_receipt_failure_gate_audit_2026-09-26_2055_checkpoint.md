# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 20:56 UTC

Read-only full-family scan: 2026-09-26T20:56:32+00:00 to 2026-09-26T20:56:59+00:00.

- Valid receipts: 17,684 / 29,940; scan-start inventory 17,684, post-scan inventory 17,687, runner state 17,687.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3152, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 777/2,495 | 31.14% | FAIL |
| longsleep | 927/2,495 | 37.15% | FAIL |
| napping | 781/2,495 | 31.30% | FAIL |
| shortsleep | 863/2,495 | 34.59% | FAIL |
| sleep_apnea | 658/2,495 | 26.37% | FAIL |
| sleep_efficiency | 955/2,495 | 38.28% | FAIL |
| sleep_timing | 978/2,495 | 39.20% | FAIL |
| sleepdur | 800/2,495 | 32.06% | FAIL |
| sleepiness | 907/2,495 | 36.35% | FAIL |
| snoring | 871/2,495 | 34.91% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-26_2055_checkpoint.json`.
