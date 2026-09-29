# Frailty LAVA resource and storage snapshot — 2026-09-26 22:08 UTC

Read-only host observation. No process, scheduler, receipt, claim or lock was changed.

- At 22:08:58 UTC, the repository volume had 2.8 GiB free of 228 GiB; the external analysis volume had 1.5 TiB free of 3.6 TiB. Repository storage remains below the 20 GiB project gate.
- `memory_pressure` reported 40% system-wide free memory. The host has 8 GiB physical memory; see the prior full resource sample at `lava_resource_telemetry_2026-09-26_2146.md` for CPU/load and swap detail.
- CPU utilization and disk throughput were unavailable: `iostat` failed because `sysctl(kern.boottime)` access is restricted. The process-list tool also remained unavailable; independent audit liveness checks passed for both worker PIDs at 22:08:14 UTC.
- The 22:09:16 runner snapshot reported 18,094/29,940 receipts, two effective workers from six requested, zero launches and stale recoveries. Workers 72534 and 72535 held distinct claims for `shortsleep` locus 1635 and `sleep_apnea` locus 1635. The persisted safeguard remains at two workers after swap reached 3,883/5,120 MiB.

This runner counter is not a full integrity scan. The latest independent full scan completed at 22:08:14 UTC and validated 18,090 receipts.
