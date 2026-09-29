# Track B recovery inventory

Historical run: `f8a3de8642698a84cddaf59aa2d94e518a5994139fe8fc2bb5ff530ccdaf5066`. Source: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/FAILED-PARTIAL-2026-09-01-sleep_gwas_atlas`. Backup source: `/Volumes/Extreme SSD/Codex Archive/2026-09-19/codex-2026-09-01-sleep_gwas_atlas.tar`.

| Component | Observed state |
| --- | --- |
| `B.full.tsv.gz` | Retained; 5,514,399 tested rows; SHA `f0eb8b9e8a5b3b4a0ea6151ed75473b0be3d311e463824889a050d83beaaac2e` |
| `B.provenance.json` | Retained; SHA `fb6565e7a165297a1b733882ba0375fc0668018ddff3887abd3178ac120298d3` |
| V1 contract/input gate | Retained; SHA `28869d2206c8930c68f253cfa3cc6efb093abc89d666ab6d797abe8e88eb3c95` / `79d2ae0faa0beb6368496617ab4c88a580a34e7f459f4825d5586e31eecb2101` |
| V2 terminal gate | Retained; SHA `9e9a40aa47711fdfb8191ca8cdfcf84828ae5d3f2715578b9bf91e939260ac97`; matches legacy sidecar field |
| Exact materializer/runner/coordinator | Present in backup tar; all 5/5 expected code hashes match |
| Original aligned input | Intentionally removed; declared SHA `59c21df0bf237e8ede77cfcdbf7b37b54444bc852d39991c88ba74465b2f02e7`; 5,514,399 rows; exact bytes not independently regenerated |
| Nuisance checkpoint | Original bytes embedded in cleanup plan; SHA `f51d16d39329b5d739f10a56809e418e448f02e664db137fcb59d888015ed243` |
| Shards | 276 executed; original shard RDS/individual receipts intentionally removed after attested cleanup |
| B full-scan monitor | Exit 0; peak 1,889,976,320 bytes; output hash matches |
| B full-scan log | 276/276 shard-completion markers and terminal success line |
| Frozen insomnia/ADHD inputs | Archived dense-source hashes match current source cards: {'insomnia': True, 'adhd': True} |
| B-to-v3 row schema | Different: legacy 17-field full ledger versus v3 10-field pair output; requires explicit checked adapter |
| Cleanup receipt | 562 listed work artifacts removed; completion receipt bound to plan: True |
| Current Brain6 v3 admission | No frozen importer acceptance; protected slot remains empty |

Search scope: current v03 and main worktrees, mounted Extreme SSD filename inventory, failed-partial archive, 2026-09-19 backup tar, older 2026-09-04 tar and handoff ZIP, relevant Git paths/history, and shell-history pointers. Unmounted/offline media were not inspected.
