# FI×sleep LAVA launch deferred — 2026-09-25 09:08 UTC

## Decision

The 29,940-slot FI×sleep LAVA family remains paused. I did not release pause
holds or start either four or six workers: the shared host has exhausted its
configured swap while four Brain6 R jobs and a separate SMFVI analysis are
active. The earlier six-worker interval had also measured below the recent
four-worker interval and triggered the documented fallback condition. Preserve
the current pause until host pressure improves; do not stop unrelated jobs.

## Runner state

Read-only state under
`/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1`
reported `PAUSED_AFTER_WORKERS_EXITED`, 8,176/29,940 receipts, 21,764 pause
holds, zero workers, zero launch failures, zero stale-claim recoveries, and no
`runner.lock`. The pause control remains
`pause_20260925T060720Z_b646` with four workers configured. The latest full
receipt identity/checksum and hold-coverage audit remains 2026-09-25 06:23:41
UTC; this snapshot did not re-audit receipt contents or count claim files.

No receipt, claim, concurrency control, scientific input, or frozen analysis
setting was changed during this check.

## Host snapshot

At 09:08 UTC the eight-core host reported load averages 10.47 / 9.57 / 9.35,
9.32% idle CPU, 7,557 MiB used physical memory, 71 MiB unused, and 4,034 MiB
in the compressor. `vm.swapusage` reported 9,216/9,216 MiB used (zero free);
`memory_pressure -Q` reported 23% system-wide free memory. Four Brain6
`run_power_optimized_sleep_univariate_v1.R` workers and one SMFVI
`run_component_ablations.py` process were active; no LAVA worker or coordinator
was present. The internal volume had 4.8 GiB free and the external analysis SSD
had 1.6 TiB free.

These measurements are host-wide. With swap exhausted and competing analyses
still active, launching LAVA would add memory and I/O load without a safe
headroom margin. The current state therefore remains paused; no new throughput
or completion estimate is available.
