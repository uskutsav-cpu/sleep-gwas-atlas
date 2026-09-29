# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 20:16 UTC

Read-only full-family scan: 2026-09-26T20:16:26+00:00 to 2026-09-26T20:16:46+00:00.

- Valid receipts: 17,401 / 29,940; scan-start inventory 17,401, post-scan inventory 17,402, runner state 17,402.
- Receipt, claim, duplicate identity, duplicate claim issues: 0, 0, 0, 0.
- Active worker PIDs: [72534, 72535]; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3086, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 758/2,495 | 30.38% | FAIL |
| longsleep | 906/2,495 | 36.31% | FAIL |
| napping | 762/2,495 | 30.54% | FAIL |
| shortsleep | 839/2,495 | 33.63% | FAIL |
| sleep_apnea | 644/2,495 | 25.81% | FAIL |
| sleep_efficiency | 928/2,495 | 37.19% | FAIL |
| sleep_timing | 956/2,495 | 38.32% | FAIL |
| sleepdur | 781/2,495 | 31.30% | FAIL |
| sleepiness | 887/2,495 | 35.55% | FAIL |
| snoring | 849/2,495 | 34.03% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `analysis/lava_full_receipt_integrity_2026-09-26_2016_checkpoint.json`.
