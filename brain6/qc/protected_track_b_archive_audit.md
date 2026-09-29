# Protected Track B archive audit

Audit date: 2026-09-23. This is a read-only provenance audit; the archive was not extracted or modified.

## Archive examined

- Path: `/Volumes/Extreme SSD/Mac Offload 2026-09-04/Research Downloads/sleep_gwas_atlas_track_b_handoff_20260904-101734.zip`
- Size: 185,029,184 bytes.
- SHA256: `3119b03ae8f1d96844e8678317d93f54f17be399b6f80ca4269535c177fefb3a`.
- The archive root identifies a different workspace, `sleep_gwas_atlas`, rather than this `sleep-gwas-atlas-v03` checkout.

## Evidence in the archive

- `results/track_b/00_repository_checkpoint.json` labels itself as a checkpoint before Track B pair selection or local-result access (`checkpoint_utc=2026-09-01T20:41:08Z`). It reports an 8 GiB memory machine, insufficient for several planned downstream stages.
- `results/track_b/CURRENT_STATUS.md` says the run was stopped before any new post-discovery aggregate or conditional checkpoint was published. Its next planned command was to execute PLACO pair analyses.
- The ZIP central directory contains only the directory entry `package/sleep_gwas_atlas/results/track_b/pleiotropy/`; it contains no files in that directory. The archive has policies, scripts, readiness and blocker artifacts, but no protected PLACO result table.

## Conclusion

This archive is not a source for the protected `insomnia__adhd` PLACO output and none of its setup/readiness files are admitted as empirical Brain6 results. The current v3 lock continues to protect that legacy pair from rerun and keeps it in the across-track correction. This audit does not establish that no other copy exists elsewhere; the previously audited current checkout and production scratch still contain no publishable protected-pair output. Do not alter the protected Track B artifacts.


## Current worktree cross-check (2026-09-23)

A read-only tracked-tree and live-worktree check was performed after the archived-handoff audit. The legacy main worktree (`/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas`, HEAD `659d01cf5b10a6debf1841404e7e01ef4750f755`) and the current v03 worktree both contain exactly four files under `results/track_b/pleiotropy/`: `contract.lock.json`, `input_gate.lock.json`, `input_gate.tsv`, and `readiness_gate.tsv`. `git status --short -- results/track_b/pleiotropy` is empty in the main worktree. No PLACO variant results or importable legacy result manifest are present there. This bounded filesystem check supplements, but does not imply a search of unmounted storage.

## Additional bounded artifact search (2026-09-23 19:23 UTC)

A read-only filename/path search checked `/Volumes/Extreme SSD/brain6-work`, `/Users/swethasunilkumar/Documents/Codex`, and `/Volumes/Untitled` for filenames containing `insomnia` or `ADHD` associated with PLACO/pleiotropy/Track B paths, plus directories named `insomnia__adhd`. No protected PLACO output, publishable result table, or receipt was found. The search does not inspect content of arbitrarily named files and does not cover offline/unmounted storage; it does not establish global nonexistence. The protected pair remains excluded from the published v3 master and family completion remains blocked on a source-verified authorized result.

## Additional failed-partial archive candidate check (2026-09-24)

A subsequent read-only check found a separate mounted archive that was not covered by the earlier bounded filename search:

- Archive root: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas`.
- Target pair-B output: `results/track_b/pleiotropy/results/placo/B.full.tsv.gz` (insomnia–ADHD), 274,416,697 bytes, SHA256 `f0eb8b9e8a5b3b4a0ea6151ed75473b0be3d311e463824889a050d83beaaac2e`.
- Sidecar: `B.provenance.json`, SHA256 `fb6565e7a165297a1b733882ba0375fc0668018ddff3887abd3178ac120298d3`.
- Streaming validation read the entire gzip successfully: CRC passed; 5,514,399 rows matched the sidecar; all rows had the declared pair and analysis ID, `TESTED` status, valid P and q values, and zero numerical errors. The output SHA matched the sidecar.
- The archive contract lock SHA matches the current protected contract (`28869d2206c8930c68c253cfa3cc6efb093abc89d666ab6d797abe8e88eb3c95`). The policy SHA and pinned PLACO 0.2.0 source SHA also match the current frozen policy/source. The archived input gate's insomnia and ADHD dense-input hashes, sizes, and row counts match the current pre-results gate.

This is **not an admissible Track B result for the current family**, despite the intact data file and sidecar. The provenance's `input_gate_lock_sha256` (`9e9a40aa47711fdfb8191ca8cdfcf84828ae5d3f2715578b9bf91e939260ac97`) does not match the archived lock file itself (SHA256 `79d2ae0faa0beb6368496617ab4c88a580a34e7f459f4825d5586e31eecb2101`); the aligned input identified by SHA256 `59c21df0bf237e8ede77cfcdbf7b37b54444bc852d39991c88ba74465b2f02e7` is not present in the archive; and the sidecar's runner/materializer source hashes (`18b01861165925227d59f7c55fc4b4549c6a17775100458845dc76ec56d002f9`, `04e9f2f96c81fc1d612ea6c2aa307041594eb865318f87b2a4a469fa1efb2a54`) do not match the available worktree or archived code snapshots. Most decisively, the archive's `continuations/post_lava_terminal_v2/implementation_compatibility_v4/AMENDMENT.md` explicitly records that the old staged full-P result and shards were retained historically because full-scan monitor evidence was unsuccessful; those shards were rejected for adoption under the corrected V4 fingerprint, which required a new computation.

Therefore the output is preserved as a historical failed-attempt artifact only. It was not copied into the current repository, the protected Track B lock or four-pair published master was not changed, and no Track B scientific result is promoted. A fresh authorized execution is not initiated from this repair task; availability of this archive does not override the protected-pair no-rerun rule. This finding narrows the prior absence statement: a candidate file exists in this separate failed-partial archive, but no source-verified, terminal-valid, importable Track B result is available.
