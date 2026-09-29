# FI×sleep LAVA full-family failure-gate audit — 2026-09-26 19:55 UTC

Read-only full-family scan: 2026-09-26T19:54:09+00:00 to 2026-09-26T19:55:09+00:00.

- Valid receipts: 17,318 / 29,940; scan-start inventory 17,318, post-scan inventory 17,320; runner state at scan end 17,320.
- Receipt, claim, duplicate identity, and duplicate claim issues: 0, 0, 0, and 0.
- Active worker PIDs: 72534, 72535; liveness: {'72534': True, '72535': True}; requested/effective workers: 6/2.
- Launch failures / stale claim recoveries: 0 / 0.
- Frozen analysis lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`; prepared-input manifest SHA-256: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`.
- Process failure categories: `{"all_phenotypes_negative_variance": 3075, "no_specified_snps_in_reference": 252, "other_process_failure": 1}`.

| Sleep trait | Failure lower bound | Failure fraction | Frozen 1% gate |
|---|---:|---:|---|
| accel_sleep_duration | 1,786/2,495 | 71.58% | FAIL |
| chronotype | 1,326/2,495 | 53.15% | FAIL |
| insomnia | 756/2,495 | 30.30% | FAIL |
| longsleep | 903/2,495 | 36.19% | FAIL |
| napping | 759/2,495 | 30.42% | FAIL |
| shortsleep | 836/2,495 | 33.51% | FAIL |
| sleep_apnea | 639/2,495 | 25.61% | FAIL |
| sleep_efficiency | 923/2,495 | 36.99% | FAIL |
| sleep_timing | 952/2,495 | 38.16% | FAIL |
| sleepdur | 778/2,495 | 31.18% | FAIL |
| sleepiness | 882/2,495 | 35.35% | FAIL |
| snoring | 843/2,495 | 33.79% | FAIL |

All 12 frozen trait gates fail. Partial-family local estimates remain inadmissible; this does not establish absence of local sharing. The run remains active for receipt completeness and final collation. This is an integrity/failure-gate audit, not a reinterpretation of the analysis.

Machine-readable evidence: `frailty_paper/analysis/lava_full_receipt_integrity_2026-09-26_1953_checkpoint.json`.
