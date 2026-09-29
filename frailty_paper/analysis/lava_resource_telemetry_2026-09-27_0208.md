# FI×sleep LAVA resource and throughput checkpoint — 2026-09-27 02:08 UTC

Read-only host and runner-state sampling; no worker or scheduler changes.

- Full audit completed 02:08:26 UTC: 19,116/29,940 receipts validated; inventory reached 19,118 and runner state 19,117, leaving 10,823 according to the runner. Two active claims on distinct napping/shortsleep locus-1837 jobs passed PID-liveness checks. Six workers requested, two effective; zero launch failures and stale recoveries.
- Runner's last recorded adaptive safeguard remains 3,883/5,120 MiB swap, which reduced concurrency from three to two. Direct swap refresh is denied by sysctl permissions. `memory_pressure -Q` at 02:08 local reported 52% system-wide memory free.
- Load averages were 21.01 / 23.85 / 20.57 on an 8-CPU host. Process-level CPU and RAM readings could not be obtained: `ps` and `top` are denied by the sandbox.
- Disk: internal volume 2.7 GiB free (below the 20 GiB write gate); external analysis SSD 1.5 TiB free.
- From the 02:00:34 UTC state (19,083) to the 02:08:18 UTC runner state (19,117), 34 receipts arrived in 7m44s (~265/hour). This short interval is not a stable forecast. The earlier recorded serial-era rate is about 81/hour; the comparison is not controlled.
- Frozen lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70` and prepared-input manifest SHA-256 `8c9921bd154c8b9bc4333db754d14e19f0af4dcc060c3ecfd91be68545abe8e2` remain unchanged.
