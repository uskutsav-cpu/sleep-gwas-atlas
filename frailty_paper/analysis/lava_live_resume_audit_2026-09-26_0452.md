# FI×sleep LAVA live state and event-log audit

Snapshot time: 2026-09-26T04:52:42.572056+00:00

## Coordinator state

- Verified receipts: 9,794/29,940; remaining: 20,146. This count is read from the coordinator state file and is not a fresh receipt-file audit.
- Effective workers: 2 of 6 requested. State lists worker PIDs 78548, 78549, 78550, 78551, 78552, 78553.
- Active claims: insomnia locus 560 (PID 78548), napping locus 985 (PID 78549). Claims are distinct.
- Launch failures: 0; stale recoveries: 0.
- Frozen analysis-lock SHA-256: `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.
- Adaptive resource fallback: swap use reached 2772/3072 MiB (≥90% of allocated swap); reducing three workers to two.

## Append-only event log

The full current parallel event log contains 4,977 valid JSON events: 4,922 receipt events, 54 worker-finished events, and 1 scheduler-pause events. Malformed lines: 0. Receipt events map to 4,922 unique `(trait,locus)` pairs; duplicate pair events: 0. Last receipt event: 2026-09-26T04:52:43.312724+00:00 (napping locus 985). This is an event-ledger audit, not a receipt-file checksum audit.

## Measured throughput and resource sample

From the saved 04:41:40 UTC snapshot (9,764 receipts) to this coordinator state (9,794), 30 receipts completed in 11.04 minutes, or approximately 163.0 receipts/hour. At that short-interval rate, 20,146 receipts imply 123.6 hours (5.1 days); this is a provisional observed-rate estimate. The saved concurrency record reports a recent serial estimate of 81 receipts/hour; this interval is about 2.0× that reference. The six-worker trial remains higher at its measured initial rate and was adaptively reduced under swap pressure; no six-worker rate is assumed for this ETA.

A system sample at 04:48 UTC reported 43% system-wide memory free; `/Volumes/Extreme SSD` had 1.6 TiB free, while the internal data volume had 6.0 GiB free at 97% used. Current CPU and live swap values could not be read from this restricted command context (`top`, `iostat`, and `sysctl` telemetry calls were denied); the coordinator still records the prior 2,772/3,072 MiB (90%) swap downshift trigger. No concurrency increase or worker restart was attempted.

## Snapshot checksums

- runner state: `6a879be39fbe93e19d1e04a3c2872807a52b6af7b13d2710a8e0bd9efdbd118f`
- parallel event log: `13bcd1580f64ec0a338214a459596187bdc021be409821305bbe761af8acf250`
- input manifest: `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`
- frozen LAVA configuration (`frailty_paper/config/lava_frailty_sensitivity_v1.yaml`): `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` (matches the runner state's analysis-lock hash)
- live runner lock (`runner.lock`, coordinator PID 78400): `83bc47cfa58a7b0b7b886d32cc25882bc6bcc6448a4662af070722b339de0ca0`
- saved concurrency/resource record (`concurrency.json`): `84a13ab6dc25a3f13d3f9df4951d1df5f9dcac3bdfbe5fbd29c8960a89ec1bc3`

No scientific settings, receipt files, claims, or QC thresholds were changed.
