# Frailty LAVA resource and storage snapshot — 2026-09-26 22:00 UTC

Read-only host and runner observation. No process, scheduler, receipt, claim, or lock was changed.

## Host and storage

- `memory_pressure` at 22:00:15 UTC reported 49% system-wide free memory on an 8 GiB host. A contemporaneous `vm_stat` sample showed 4,107 free 16 KiB pages (about 64 MiB), with memory actively compressed; `memory_pressure` is the system-wide pressure indicator.
- Repository volume (`/System/Volumes/Data`) had 1.9 GiB free of 228 GiB at 22:00:15 UTC. External analysis volume (`/Volumes/Extreme SSD`) had 1.5 TiB free of 3.6 TiB.
- CPU utilization and disk throughput were not available in this sample: `iostat` failed because access to `sysctl(kern.boottime)` is restricted. Process-list inspection also returned “Cannot get process list”; the receipt auditor’s per-PID liveness checks are the process evidence used for the preceding 21:52 audit.

## Runner state

At 22:00:33 UTC the runner state reported 18,058/29,940 receipts, two effective workers from six requested, zero launch failures and zero stale recoveries. Workers 72534 and 72535 held separate claims for `sleepdur` locus 1631 and `sleepiness` locus 1631. The persisted fallback remains two workers after swap reached 3,883/5,120 MiB at 21:52. This runner counter is not a full receipt integrity audit; the latest independent scan validated 18,027 receipts through 21:52:39 UTC.

The repo volume remains below the project’s 20 GiB storage gate. The external volume has ample space for analysis outputs, but the live workspace contains unrelated mixed-worktree data; no cleanup or relocation was attempted during this observation.
