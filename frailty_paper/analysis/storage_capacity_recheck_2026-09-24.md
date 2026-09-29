# Storage capacity recheck

Audit time: 2026-09-24 00:13 UTC.

- Repository/internal filesystem: 3.4 GiB available (below the configured 20 GiB minimum).
- External SSD `/Volumes/Extreme SSD`: 1.7 TiB available.
- The configured raw, harmonized, munged and frailty-data destinations on the external SSD all exist and are resolvable. The registered source-data directory remains on the external volume.
- No large download, harmonization, or analysis output was started during this check. Continue routing large inputs, temporary files and derived outputs to the external volume and enforce the 20 GiB gate against each actual destination before a large stage.

## External mount spot-check (2026-09-24 00:44 UTC)

The external SSD remains mounted at `/Volumes/Extreme SSD` with 1.7 TiB available; the internal filesystem reports 3.3 GiB available, still below the 20 GiB gate. `frailty_paper/data` resolves to `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data`. The three files reported as deleted in the Git worktree are present through that symlink and match the acquisition manifest checksums: `eqtl_catalogue_dataset_metadata_r7.tsv` (75,969 bytes; `70b0b6eac8a607568d11b3f2d8ad6cf6e9bb430ffafa4883627797d8967785ea`), `brainscope_eqtl_links.json` (2,811 bytes; `eabbfd27c3952fa33028a2790867c8af0fbadd65ea6ae3bc3864272fb3b2d262`), and `brainscope_key_resource_links.json` (3,363 bytes; `b20cff40752554ee386c3325d636f2a141e07a7b1dbea2a366cadecfeb9e35d4`). This is a spot-check of those three files, not a new full 347-file manifest verification. No output-producing stage was run.

## Full manifest verification (2026-09-24 00:50 UTC)

With the external volume mounted, reran `frailty_paper/.venv/bin/python frailty_paper/scripts/16_verify_acquired_resource_manifest.py --repo . --external-storage-root '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'`. Result: `RESOURCE_MANIFEST_OK rows=347 files_verified=347` (exit 0). The 00:44 section above records the earlier three-file spot-check; this full run supersedes its statement that no full verification had been run. This is a read-only integrity check and no output-producing analysis stage was run. At recheck, internal free space was 3.3 GiB and the external SSD had 1.7 TiB available.
