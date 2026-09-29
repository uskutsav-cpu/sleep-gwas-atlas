# Storage capacity recheck — 2026-09-26 20:04 UTC

Immediately before any further large derived-data operation, `df -h` reported:

- External analysis workspace `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace`: 1.6 TiB available on a 3.6 TiB volume, above the locked 20 GiB minimum.
- Internal filesystem `/`: 2.0 GiB available on a 228 GiB volume, below the 20 GiB minimum.

Large LAVA and derived-analysis writes remain routed to the external analysis workspace. The live LAVA run is active there; no directories or symlinks were changed by this capacity check. Recheck actual destination capacity immediately before any subsequent large analysis.
