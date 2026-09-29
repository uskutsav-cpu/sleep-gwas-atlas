# Repository, collection and review audit — 2026-09-23

## Repository state

- Repository: `/Users/swethasunilkumar/Documents/Codex/2026-09-19/can/work/sleep-gwas-atlas-v03`.
- At audit start, current branch was `frailty-paper-v1`, `HEAD` was `49ce9a2`, and `HEAD` was attached. The local `main` branch is in a separate worktree; it was not changed.
- The working tree already contained extensive other work: 19 modified paths, six deletions, and 1,024 untracked paths at the initial status snapshot, including Brain6 chunk/checkpoint files. Those changes were left untouched.
- Inventory includes the established `config/`, `data/`, `results/`, `scripts/`, `extensions/brain6/`, and `frailty_paper/` trees. The existing repository provides Track B replication, LAVA, PLACO+, SuSiE/coloc, QTL and Brain6 workflows; the frailty paper area contains a review pipeline, source manifests, frozen analysis plan, FI and latent-factor LDSC outputs, and provisional figures/manuscript.

## Capacity and storage

At the audit snapshot, `df -h` reported about 8.4 GiB free on the internal filesystem and 1.8 TiB free on `/Volumes/Extreme SSD`. The external `analysis-workspace` occupied about 39 GiB. The locked storage gate passed with 1,797 GiB free at the external workspace. Raw GWAS and derived harmonized/munged paths resolve under the configured external root. Recheck capacity and route `TMPDIR` before each large-output stage; the internal volume remains below the 20 GiB workspace minimum.

## Collection inventory and integrity

`frailty_paper/manifests/all_acquired_resources.tsv` contains 347 unique resource IDs and 347 unique paths:

| Resource class | Files | Manifest bytes | Current verification |
|---|---:|---:|---|
| GWAS summary statistics | 32 | 13,447,043,914 | Every listed path, size and SHA-256 passed the manifest verifier using the explicit external-storage root. The 20 non-sleep GWAS sources are structurally scanned in `gwas_source_scan.tsv`; the 12 locked sleep sources are separately scanned in `gwas_source_scan_sleep_panel.tsv`. Structural validity does not imply every source is scientifically usable. |
| PubMed XML batches | 311 | 1,320,887,288 | Every manifest size/hash passed. A fresh streaming XML parse found 60,477 `PubmedArticle` plus 512 `PubmedBookArticle` nodes (60,989 total); all 311 files contained at least one record and parsed successfully. |
| Public-resource index/metadata files | 4 | 107,521 | Manifest size/hash passed. These are indexes/metadata, not acquired eQTL or BrainSCOPE QTL summary statistics. |

The PubMed builder's current deduplicated queue has 56,092 unique records and links 4,897 duplicate occurrences. Review validation reports zero title/abstract decisions, zero full-text records/decisions, zero exclusions and zero included studies. Licensed Embase, Scopus, Web of Science and optional PsycINFO exports are absent. The review is not complete.

Acquired data include the 12 locked sleep GWAS, primary FI, six selected aging/context GWAS, five physical-component candidate GWAS, and eight Catalog accessions GCST90624046–GCST90624053 (seven mapped frailty factors plus pneumonia). The eight latent files have source and duplicate-key audits; the pneumonia file retains 466,894 non-finite/missing beta or SE rows. Physical-component candidate files pass structural scans but remain ineligible because deposit-level build/effect/model provenance and exact phenotype equivalence are unresolved. The exact UKB Fried Frailty Score item has been identified but no verified full file is acquired. The HFRS custom endpoint is unverified; the standard FinnGen route requires the official access form/email, and R18_SENILITY is not an equivalent proxy.

The collection log contains a truncated PubMed network response and a later `--config: command not found` message from `07_run_all_collection.sh`. The current script passes `bash -n`; the logged message does not reproduce from a syntax-only check. `launchctl` returned no matching collection job, while process-table access through `ps`/`pgrep` was denied in this environment. No live collector handle was available to poll and no collection/download was restarted. The current 347-file inventory passed verification, so no redownload was needed for Phase 1.

## Existing analyses and validation

- Analysis-plan lock: `ANALYSIS_PLAN_LOCK_OK`, 12 sleep traits and inherited multiplicity family 396.
- Source-level overlap ledger: valid, 34 rows; exact participant intersections remain unknown.
- Primary sleep × FI: 12 pairs, 9/12 inherited all-396 BH q≤0.05; same-family Bonferroni sensitivity agrees for those nine.
- Secondary sleep × latent frailty: all 84 pairs completed; 48/84 BH q≤0.05 in the locked secondary family. Results remain sensitivity-only because exact cohort overlap is unknown.
- Aging-context extraction: 72 read-only frozen-atlas pairs, not replication.
- All 58 tests in `make -C frailty_paper test` passed. Plan, overlap and review validators passed; the latter correctly reports that no records have been screened.
- No eligible independent pairwise sleep–frailty replication is established. Consequently, frailty-specific local/LAVA, PLACO/shared-locus, fine-mapping, coloc and QTL/cell-type stages remain gated or not justified.

## Current execution order

1. Keep the verified acquisition manifest and raw sources immutable; retain the scientific source/provenance failures in the per-source audits.
2. Resolve only genuine source-level blockers that could unlock a predeclared phenotype (exact Fried file access, physical-component provenance, or HFRS endpoint/access). Do not substitute non-equivalent outcomes or pass an access form without user authorization.
3. Complete manual database searches through authorized institutional access, import untouched exports, then begin dual screening and extraction. Do not infer completion from the PubMed-only queue.
4. Proceed with additional harmonization/h²/global rg only for datasets that pass the frozen source, phenotype, ancestry/build, overlap and QC gates.
5. Seek paired, independent replication for headline pairs before any frailty-specific local-sharing or mechanistic work. Keep downstream stages gated if no valid pair/locus exists.
6. Build only supplementary tables and manuscript claims supported by frozen, validated outputs; finish reporting checklists, reviewer audit, frozen package and final readiness report after the evidence package is mature.
