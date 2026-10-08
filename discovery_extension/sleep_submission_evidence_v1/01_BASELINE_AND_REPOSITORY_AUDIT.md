# Baseline and repository audit

Execution date: 2026-10-08 (America/Chicago). Scope: frozen sleep phenome extension only; no manuscript work.

The research branch starts at authoritative extension main `303c905c229104b7bf623df6de7431abb87e0727`. The historical core checkpoint is `22df92c8d6894f93c50561ae8a86acf2f48f50cb`; the frozen extension records 12 sleep traits and 100 selected external phenotypes. These are separate from the immutable 45-trait core and its 396 comparisons. Main is not merged or modified by this package.

Historical source contracts, result tables, selection locks, provenance JSON, acquisition receipts and replication LDSC logs were examined. Source identities must come from `config/analysis_panel.tsv`, not the older `config/traits.tsv`: the older file retains obsolete aliases and discrepant sleep-apnea counts. The new metadata table preserves the locked panel values.

The source-free numerical audit generates actual counts elsewhere in this package. The full historical result universes are 1,200 discovery and 217 replication candidates. The frozen report's 603 significant discovery comparisons and 23 qualifying validation pairs are hypotheses to check, not input constants.

`git worktree list --porcelain` identifies two earlier local sleep-project worktrees. Only sleep namespaces were eligible for restoration; Brain6 and Frailty scientific data were excluded. Git LFS lists only Brain6 files. Git history contains no tracked versions at the ten historically missing core result paths; the checkpoint tree itself stores `.gitkeep` files at those output directories. Search details are in `sources/provenance_recovery_scan.json`.

`/Volumes` contains only `Untitled -> /`; the expected Extreme SSD is absent. The filesystem has approximately 1 GiB free at task start. A 246.91 GiB dense mirror is infeasible and was not started. A bounded remote-access pilot consumes at most 4,096 bytes from each of three pinned objects and retains metadata only; its result is in `sources/provenance_bounded_access.json`.

Decision: proceed with validated frozen processed results, historical logs and receipt checks, retaining explicit native-rerun and source-restoration blockers. No success of a software test demonstrates a new GWAS run.

The initial access pilot used system curl with verified TLS for three HTTP 206 range reads (12,288 bytes total), retaining no raw bodies. Python HTTPS initially failed local certificate trust; SSL verification was not disabled. Exact Pan-UKB version and FinnGen generation headers match historical ETag/length. MVP ETag/Last-Modified metadata differ, so subsequent complete-body checks tested whether scientific compressed bytes changed. Initial metadata evidence remains separate in `sources/provenance_bounded_access.json`.

Current whole-body verification completed for both MVP sources using streamed SHA-256 and MD5 with verified TLS and no raw retention. GCST90479148 abdominal pain: 575,939,664 bytes in 1,150.383 seconds; GCST90479330 allergy: 612,834,482 bytes in 1,312.911 seconds. Both digests match the historical full-stream receipts exactly (`PASS_SOURCE_BODY_HASH` only; native LDSC remains `BLOCKED`). Therefore these two current compressed objects retain original byte identity despite altered HTTP metadata. The full bodies plus the earlier pilot and range reads consumed 1,208,754,180 bytes, below the 2 GiB campaign cap. Current receipts are `logs/provenance_remote_hash_gwas_catalog_GCST90479148.json` and `logs/provenance_remote_hash_gwas_catalog_GCST90479330.json`; exact digests are also in `tables/source_integrity_ledger.tsv`. This establishes current remote source integrity for two objects, with zero raw sources recovered and zero native GWAS/LDSC reruns.
