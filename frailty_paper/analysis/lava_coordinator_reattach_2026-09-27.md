# LAVA coordinator reattachment — 2026-09-27

## Incident and recovery

At 20:54:29 UTC, runner state stopped updating at 21,366/29,940 and showed zero registered workers/active jobs. The coordinator session later exited with `PARALLEL_FINISHED verified=21366/29940 failures=0`. However, the append-only receipt event log continued through 20:58:29, and claim files at 20:58:30 showed jobs assigned to worker PIDs 15810 and 15811. The 20:51 full audit had also confirmed both PIDs live. `runner.lock` was absent after the coordinator exit. The reason the coordinator concluded there were no live workers is not established.

At 21:03:37 UTC, the existing coordinator script was restarted with `--workers 2` after verifying process-group PIDs 15810 and 15811 were still present. It re-registered those existing workers and reported `PARALLEL_LAUNCHED workers=2 pids={'1': 15810, '2': 15811} verified=21412/29940 stale_claims_recovered=0`; it launched no additional workers. Runner state and `runner.lock` then resumed updating. At 21:05:01 UTC it reported 21,420 receipts, two active claims on PIDs 15810/15811, two requested/effective workers and zero stale recoveries.

## Independent verification after reattachment

The read-only full audit ran 21:04:23–21:04:56 UTC and validated 21,416 receipts. Inventory advanced to 21,419 during the scan; runner state reached 21,420 shortly afterward. Receipt, claim, duplicate-receipt and duplicate-claim issues were all zero. Both active worker PIDs were live. Failure categories were 3,726 all-phenotype negative-variance, 252 no-reference-SNP and two other process failures; all 12 frozen 1% trait gates fail. Evidence: `lava_full_receipt_audit_2026-09-27_2104.json` and `lava_runner_state_snapshot_2026-09-27_2104.json`.

The coordinator restart reset its throughput window (`measured_receipts_per_hour: null`, zero elapsed window), so no ETA is reported until a complete post-reattachment window is measured. This restored coordinator monitoring without changing analysis inputs, receipts, locus assignments, failure classifications or locked thresholds. The earlier 20:38 resource downshift to two workers remains the reason for the present two-worker limit, even though a fresh invocation no longer reports that fallback string as active.
