# FI×sleep LAVA resource and throughput checkpoint — 2026-09-27 01:51 UTC

Read-only host and runner-state sampling; no worker or scheduler changes.

- Latest full receipt audit completed 01:50:04 UTC: 19,009/29,940 valid receipts; runner state reached 19,011, leaving 10,929. Two active claims were on distinct traits/loci; the auditor's PID-liveness checks returned true for both. Six workers requested, two effective; zero launch failures and zero stale-claim recoveries.
- The runner's last recorded adaptive safeguard remains swap 3,883/5,120 MiB, which had reduced effective concurrency from three to two. The OS `memory_pressure -Q` sample at 01:50 local reported 64% system-wide memory free. These are distinct telemetry sources and timestamps; direct swap refresh was denied by `sysctl` permissions.
- Host load averages at the sample were 13.24 / 12.78 / 14.12 on an 8-CPU host. Process-level CPU and RAM readings could not be obtained: `ps` and `top` were denied by the sandbox.
- Disk: internal volume 2.8 GiB free (below the 20 GiB write gate); external analysis SSD 1.5 TiB free.
- Since the 01:42:17 UTC state (18,954 receipts), the 01:50:04 UTC state (19,011) advanced by 57 in 7m47s, about 440 receipts/hour. This brief interval is not a stable forecast; do not extrapolate a completion time from it. The earlier recorded serial-era rate is about 81/hour, so this window is directionally faster but is not a controlled comparison.
- Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and prepared-input manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.
