# FI×sleep LAVA coordinator recovery — 2026-09-25 22:07 UTC

## Recovery

The coordinator recorded as PID 21938 had exited while workers 22177 and 22178 remained alive. The external `runner.lock` still named PID 21938; process inspection and a process-group signal check confirmed that PID was absent. Both live worker groups retained distinct atomic claims. No claim was released or worker interrupted.

Before recovery, 27 receipt files newer than the last runner-state update at 22:00:16 UTC passed `45_run_lava_sensitivity.py::verify_receipt`, including identity/status checks and nonempty logs. The replacement coordinator's startup then revalidated 8,852 existing receipts. Prepared inputs passed `43_prepare_lava_sensitivity_inputs.py --verify` with frozen manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2`; the analysis-lock SHA-256 remains `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`.

The replacement coordinator (PID 52234; attached session 78739) started with four requested workers. The resource guard reduced four to three at 22:04:15 UTC when swap reached 3,233/4,096 MiB (the 1,536 MiB absolute trigger). It then reduced three to two after swap crossed the 3,584 MiB deep-pressure trigger. Worker 3 (PID 52535) completed its in-flight locus and exited; workers 1 and 2 remained active. No worker was relaunched over a live claim.

## Verified snapshot

At 22:07:29 UTC, `runner_state.json` reported 8,861/29,940 verified receipts, two effective workers, and two distinct live claims (longsleep locus 522 and sleepiness locus 938); launch failures and stale recoveries were zero. At 22:07:34, the receipt-path count was 8,862 while the state still reported 8,861, consistent with a receipt written between the runner's periodic state update and this read. The coordinator remained live and the next reconciliation was pending. The full event log had 4,027 rows, no malformed lines, and no duplicate receipt keys.

At 22:07:11 UTC, swap was 4,646/6,144 MiB and system-wide free memory was 31%; retain the two-worker rung. An earlier 22:02 sample had 49% free memory but swap reached 3,233/4,096 MiB during coordinator startup, which correctly triggered the 4→3 reduction. The fallback thresholds and scientific configuration were unchanged. No local-sharing inference is supported: the first-pair and chronotype interim QC gates remain failed, and the full LAVA family is incomplete.

This is an incremental operational recovery, not a complete family receipt-manifest rehash or a final QC audit. The campaign and the broader paper remain incomplete; readiness remains NO-GO.


## Follow-up live reconciliation — 2026-09-25 22:11 UTC

At 22:11:12 UTC, the replacement coordinator reported 8,881 verified receipts and the external receipt-file count was also 8,881. Two effective workers had distinct claims on napping locus 523 and sleepdur locus 939. The coordinator continued updating state. Swap was 6,288/7,168 MiB and system-wide free memory 31%; the two-worker deep-pressure rung remains appropriate.
