# Storage capacity recheck — 2026-09-26 20:50 UTC

Read-only command: `df -h '/Volumes/Extreme SSD' /Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03`.

| Destination | Size | Used | Available | Capacity | Decision |
|---|---:|---:|---:|---:|---|
| External SSD (`/Volumes/Extreme SSD`) | 3.6 TiB | 2.1 TiB | 1.6 TiB | 57% | Above the 20 GiB analysis gate; retain as the destination for large intermediates. |
| Internal data volume (repository) | 228 GiB | 195 GiB | 4.9 GiB | 98% | Below the 20 GiB analysis gate; do not write large intermediates here. |

The LAVA run is configured on external storage. Capacity was sufficient at this snapshot. Recheck against the actual output destination before each new large analysis.
