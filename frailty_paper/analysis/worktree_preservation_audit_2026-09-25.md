# Frailty worktree storage-path preservation audit — 2026-09-25 06:32 UTC

## Scope

A read-only reconciliation investigated tracked deletions under `frailty_paper/data/` and `frailty_paper/review/pubmed/`. Both directory paths are symlinks into the external SSD (`frailty_paper/data` → `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data`; `frailty_paper/review/pubmed` → `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/pubmed_cache`). Git reports the child paths as deletions in this checkout view; this does not establish that the external files are absent.

## Results

The audit compared each of 21 affected tracked paths with both its resolved external target and the internal preservation copy, using the committed `HEAD` blob as the reference. All 21 internal copies are byte-identical to `HEAD` (99,480 bytes total). Twenty external target files are also byte-identical. The only changed external file is `review/pubmed/search_log.tsv`: its 2026-09-24 `frailty_neighborhood` row records 55,599 retrieved / 59,819 fetched, while the preserved `HEAD` version records the earlier 55,572 / 59,799. The external log is byte-identical to `frailty_paper/review/pubmed_search_log_current_2026-09-24.tsv`; the post-refresh reconciliation reports 55,599 full-query unique IDs, 59,819 segment records, and zero missing/extra segment IDs. The two other query rows remain unchanged.

The machine-readable path, byte-count, and SHA-256 comparison is in [`worktree_preservation_audit_2026-09-25.tsv`](worktree_preservation_audit_2026-09-25.tsv). No files were restored, overwritten, or modified on the external volume. The checksum evidence indicates preservation plus one documented search refresh, rather than an unverified loss of those tracked source records.
