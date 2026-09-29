# Storage capacity recheck — 2026-09-23

Observed at `HEAD` `7d59cdff68dfc1bfaf5f7fba7dcaeacd1366b4ac` at 2026-09-23 20:33 UTC. The same commit was present before and after the read-only checks.

## Capacity

- Internal filesystem (`.`): 228 GiB total, 194 GiB used, 5.5 GiB available (98% capacity).
- External SSD (`/Volumes/Extreme SSD`): 3.6 TiB total, 1.9 TiB used, 1.7 TiB available.
- `frailty_paper/scripts/09_require_storage.sh` passed for each configured analysis destination, reporting 1,775 GiB available at `data/harmonized`, `data/munged`, and `frailty_paper/data`.
- The project gate requires at least 20 GiB on each large-output destination. The external destinations pass; the internal filesystem remains below that minimum.

## Decision

Continue routing raw, harmonized, munged, and large derived data to the external SSD. Set temporary directories and each stage's output path explicitly to that volume, and run the storage gate again immediately before a large stage. Do not place large intermediates on the internal filesystem. No files were created on the external volume as part of this check.

## Commands

```sh
df -h . '/Volumes/Extreme SSD'
bash frailty_paper/scripts/09_require_storage.sh "$PWD" data/harmonized
bash frailty_paper/scripts/09_require_storage.sh "$PWD" data/munged
bash frailty_paper/scripts/09_require_storage.sh "$PWD" frailty_paper/data
```
