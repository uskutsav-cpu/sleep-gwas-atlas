# Review workflow

1. The exact PubMed search files and counts are in `pubmed/`; acquisition integrity and index drift are in `pubmed_reconciliation.md`.
2. Build the record-level audit and blank screening queues with:

   ```bash
   frailty_paper/.venv/bin/python frailty_paper/scripts/08_build_review_records.py
   ```

   The builder now refuses to overwrite any nonblank reviewer decision, full-text location, exclusion record, or included-study extraction. Archive and reconcile screening work explicitly before rebuilding records; do not clear the guard to force regeneration.

3. Reconcile saved XML against current official PubMed search IDs (network required; this is read-only):

   ```bash
   frailty_paper/.venv/bin/python frailty_paper/scripts/09_reconcile_pubmed_retrieval.py
   ```

4. Run licensed searches and place untouched exports in `manual_exports/` using `manual_search_instructions.md`. Imports must be added without overwriting PubMed source records and must retain original database identifiers.
   Normalize RIS/CSV exports and capture their source-file hashes with:

   ```bash
   make -C frailty_paper import-review-exports
   ```

   Then rebuild the combined record/dedup/screening tables using step 2, before any reviewer decisions have been entered. Each source record remains in `all_records.tsv`; `duplicate_sources` and `retained_record_id` preserve the cross-database deduplication audit.
   The 2026-09-24 targeted publisher audit originally logged PMID 41692760 as a potentially post-cutoff record. A later source-level reconciliation confirmed it is already present in the frozen 2026-09-22 PubMed corpus and main queue; its earlier unresolved indexing note is superseded. The record remains unscreened and must pass the registered dual-review workflow. Do not add a supplemental PRISMA row or change counts for it.

   **Supplemental later PubMed window (kept outside the frozen protocol queue):** a reconciled 2026-09-23–25 window contains 17 article occurrences and 17 window-unique records; two PMIDs already occur in the frozen snapshot and 15 are new. All 17 remain unscreened. To rebuild its separate queue and provenance:

   ```bash
   python3 frailty_paper/scripts/55_build_pubmed_window_review_records.py \
     --increment-dir '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/pubmed_cache_increments/2026-09-23_2026-09-25' \
     --reconciliation frailty_paper/review/reconciliation_2026-09-25_window/pubmed_retrieval_reconciliation.json
   ```

   Review `reconciliation_2026-09-25_window/pubmed_window_unique_review_queue.tsv` separately. The builder verifies the frozen all-records hash and reconciliation inputs, preserves the 17 occurrence rows, and creates no reviewer decisions. It does not modify `all_records.tsv`, the main screening queue, the frozen cutoff or PRISMA counts. Include these later-indexed records in a future review only through a documented protocol/search amendment; never merge them silently.
5. The current screening queue is intentionally blank. Two independent title/abstract and full-text decisions plus adjudication are supported by the columns. Do not report review completion before supplied database exports and screening are resolved.

   **Independent reviewer packets:** To prevent one reviewer from seeing the other's in-progress decisions, export separate, hashed batches from the frozen blank queue:

   ```bash
   python frailty_paper/scripts/build_independent_screening_packets.py \
     --output-dir frailty_paper/review/screening/reviewer_packets_YYYYMMDDTHHMMZ \
     --batch-size 500
   ```

   The exporter refuses a queue with any existing decision or an output path that already exists. It writes one directory per reviewer slot, 500 records per TSV by default, and a manifest binding the source queue, frozen protocol, packet row identities and immutable source fields. Packets contain no other reviewer's decision fields and the exporter makes no decisions. Keep each directory with its assigned reviewer; do not share the other packet or a partially completed canonical queue. When one reviewer's packet set is returned, validate and preview it first, then merge only that slot:

   ```bash
   python frailty_paper/scripts/import_independent_screening_packets.py \
     --packet-dir frailty_paper/review/screening/reviewer_packets_YYYYMMDDTHHMMZ \
     --slot reviewer_1
   ```

   The default is a dry run. Add `--apply` only after reviewing the reported counts. Import validates every packet, the manifest, protocol hash, source rows, record coverage and decisions; it refuses conflicting existing decisions and creates a hash-named pre-import backup before an atomic queue update. Repeat with `--slot reviewer_2` when that reviewer's independent decisions are ready. Adjudication remains a separate human step. Generated packet directories are ignored by Git and must be stored and transferred securely by the review team. This workflow does not replace reviewer assignment, human judgment, adjudication, full-text screening, or licensed database searches.

   The frozen main queue contains 56,117 records; the 17 later-window records above are additional, separate and unscreened. The frozen-queue audit found 935 empty-abstract PMIDs and rechecked every source XML occurrence; none has abstract text. All main-queue titles are present. Missing abstract text is not an exclusion reason: use the available title, preserve an unclear/awaiting-classification decision where the criteria cannot be resolved, and seek full text as indicated. See `analysis/title_abstract_queue_metadata_audit_2026-09-24.md` and `analysis/missing_abstract_source_audit_2026-09-24_post_pubmed_refresh.json`.

   When the PubMed cache is on a mounted external volume, keep the logical `pubmed/` path pointed at that cache. The source audit accepts `--external-storage-root` only for the explicitly mounted root, for example `--external-storage-root '/Volumes/Extreme SSD'`.

   A supplemental Europe PMC lookup matched all 935 blank-abstract PMIDs; none had abstract body text (one contained only a section heading). The result is provenance-tracked separately and the screening queue remains unchanged. See `analysis/europe_pmc_missing_abstract_audit_2026-09-24.md`.

   An optional deterministic keyword-priority list is available at `screening/title_abstract_triage_priority.tsv`; reproduce it with:

   ```bash
   python frailty_paper/scripts/51_prioritize_screening_queue.py
   ```

   The accompanying `screening/title_abstract_triage_priority.provenance.json` records the source queue, frozen protocol, script and output hashes, scoring groups, and record counts. The ranking uses title/abstract keyword matches only to order human attention. It makes no screening decisions, removes no records, and does not replace independent protocol-based decisions. Keyword matches are not eligibility judgments, and the ranked order may bias reviewer attention; screen every record and feel free to ignore the ranking. Rebuild only from the frozen queue and protocol, and do not treat this file as the decision queue.

6. For included observational studies, use the design/version routing and blank dual-review templates in `risk_of_bias/`. The analytical cross-sectional (2026), cohort (2025) and current case-control JBI DOCX tools are archived under `risk_of_bias/tools/` with SHA-256 values in `risk_of_bias/tool_registry.tsv`. Recheck the registry and record the exact file/checksum in each study-register row before appraisal. Genetic association studies are handled separately through STREGA/STROBE and source-specific genetic QC rather than forcing a nonmatching JBI checklist. No studies have yet been appraised.

Large XML, extracted record tables, the master CSV, queues and raw data are locally regenerable and ignored by git. Their byte counts and SHA-256 hashes are retained in the small source/output manifests.

## Independent screening packet preparation — 2026-09-26 08:35 UTC

To support the protocol's independent human decisions without exposing one reviewer's in-progress choices to the other, added separate packet export/import tools and generated 113 batches for each reviewer (56,117 records per reviewer; batch size 500). Manifest: `frailty_paper/review/screening/reviewer_packets_20260926T0835Z/packet_manifest.json`, SHA-256 `4cdce0bb868e9b3cdac33988c093100891810e02b774a9e3822959e9edf63810`; it binds the unchanged queue SHA-256 `3ae3fb3b4765c85536c3ccf64d20297a632af892e2c66e622ee665db1ea43b2e` and frozen protocol hash. Full dry-run import validated complete record coverage for both reviewer slots with zero decisions applied. The canonical queue remains unchanged: 56,117 records and zero decisions. The reviewer validator passes structurally; human decisions, licensed exports and adjudication remain outstanding. Four focused packet-tool unit tests pass.
