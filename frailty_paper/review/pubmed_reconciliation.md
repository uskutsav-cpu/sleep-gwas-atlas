# PubMed acquisition reconciliation

**Publication-date cutoff:** 2026-09-22. The original search log is preserved as `pubmed_search_log_original_2026-09-22.tsv`; `pubmed/search_log.tsv` records the resumed retrieval run from 2026-09-24. Exact query strings and segment date ranges remain unchanged. The retrieval script's `fetched` value is the sum of ESearch segment counts, which can exceed unique IDs because PubMed can return the same record under more than one publication-date bucket.

## Current retrieval — 2026-09-24

- The resumed searches returned current full-query counts of 1,152 primary, 36 genetic, and 55,599 frailty-neighborhood records. The neighborhood ESearch segment counts sum to 59,819 occurrences; the union of the 10 segments contains 55,599 unique PMIDs. Two separately retained index-update records bring the active XML source total to 61,009 occurrences across 311 XML files.
- `scripts/09_reconcile_pubmed_retrieval.py` was run after the refresh at 2026-09-24 22:34 UTC. All 12 saved date segments matched current ESearch ID sets exactly: zero missing PMIDs, zero extra PMIDs, and zero query errors. Machine-readable results and per-segment expected-ID hashes are in `reconciliation_2026-09-24_post_refresh/`.
- The later ID-aware collector pass re-queried ESearch and compared exact PMID sets for every cached/fetched batch, not only record totals. It reused all 311 XML files and wrote 12 current `esearch_ids.txt` snapshots. The external `pubmed_cache.sha256.tsv` manifest covers 342 files (311 XML batches, 12 ID snapshots and acquisition/query metadata); no XML batch required another download.
- The later ID-aware collector pass re-queried ESearch and compared exact PMID sets for every cached/fetched batch, not only record totals. It reused all 311 XML files and wrote 12 current `esearch_ids.txt` snapshots. The external `pubmed_cache.sha256.tsv` manifest covers 342 files (311 XML batches, 12 ID snapshots, query/acquisition metadata and logs); no XML batch required another download.
- Before the refresh, the final neighborhood segment contained 7,532 unique PMIDs against 7,558 current ESearch IDs: 29 current IDs were missing and three cached IDs were no longer returned. The prior segment remains preserved at the external cache's `pubmed_cache_preserved_2026-09-24/` path with its own SHA-256 manifest. The active replacement segment was downloaded atomically and rechecked; its reconciliation is exact.
- The corrected builder parsed 61,009 XML record occurrences (60,503 journal articles and 506 book articles), retained 56,117 unique review records, and linked 4,892 duplicate occurrences. There are zero manual-database imports. The blank title/abstract queue has 56,117 rows and zero screening decisions; no full texts have been assessed and no studies have been included.
- The all-acquired-resources manifest now contains 354 rows and passes the 22-column metadata audit (rechecked 2026-09-25 04:53 UTC). The recorded full-file verification covered the original 351 resources; the three subsequently registered GEO metadata assets were individually checked for path, size and SHA-256 in their source audits. The PubMed XML batches have per-file SHA-256 records in `pubmed_source_files.tsv`; the full external cache has a separate checksum manifest. The active cache resides on the external SSD and is exposed at the logical `pubmed/` path; the original internal cache is preserved as `pubmed.internal-preserved_2026-09-24/`.
- The current missing-abstract source audit passes: all 935 empty-abstract PMIDs are present in source XML, and no duplicate occurrence supplies abstract text. Evidence: `analysis/missing_abstract_source_audit_2026-09-24_post_pubmed_refresh.json`.

## Initial 2026-09-22 acquisition checkpoint (historical)

- The original snapshot has 310 non-empty XML batches. A separately stored two-record NCBI index update brings the source total to **311 XML files**, all individually SHA-256 hashed in `pubmed_source_files.tsv`.
- The retained XML contains **60,989 PubMed records**: 60,477 `PubmedArticle` nodes and 512 `PubmedBookArticle` nodes. The original extraction loop counted only `PubmedArticle`; that caused the apparent shortfall in the first audit. The corrected parser includes both record types and now falls back to `BookTitle` for book records with no `ArticleTitle`.
- The original `search_log.tsv` reports fetched counts of 1,152 primary, 36 genetic, and 59,799 neighborhood records (60,987 in the original search snapshot). Two newly indexed neighborhood-search records were fetched separately and recorded in `pubmed_update_log.tsv`; their source XML remains separate from the original 310 batches.
- The initial builder output contained **60,989 PubMed records and zero manual-database records**, retaining **56,092 records** with **4,897 duplicate occurrences**. This checkpoint is superseded by the 2026-09-24 refreshed retrieval above. Deduplication remains PMID-first, then DOI when titles agree or a unique titled candidate exists, then normalized title only when identifiers do not conflict; otherwise records remain separate. Licensed database exports are still required before the combined search is complete.
- `pubmed_master.csv` came from the earlier script that omitted book-record parsing. Use the corrected `all_records.tsv` and `deduplicated_records.tsv` for screening; the original file is preserved unchanged.
- The acquisition-record builder's `record_build_manifest.json`, `record_build_outputs.tsv` and `pubmed_source_files.tsv` record script/source/output checksums. All five non-empty registered GWAS files listed in `manifests/all_acquired_resources.tsv` passed gzip integrity checks and have non-zero sizes and SHA-256 hashes.

## First current-index check (2026-09-22; superseded)

A read-only NCBI E-utilities reconciliation was refreshed at **2026-09-22 22:30 UTC**. All 12 saved date-segment queries were checked with no query errors. The primary and genetic segment records match current ESearch IDs. The neighborhood segment has 7,539 current ESearch IDs versus 7,538 archived XML PMIDs; the two missing IDs (42771750 and 42771753) were fetched into `pubmed/frailty_neighborhood/index_update_20260922/batch_000000.xml`, separately from the original snapshot. One archived PMID falls outside its current publication-date segment, reflecting index/date-segment drift; the original XML is retained unchanged. The exact new-file checksum, PMID list and retrieval source are in that folder's `acquisition.json`; update summary is in `pubmed_update_log.tsv`.

The full machine-readable query/segment counts, expected-ID hashes and differences are in `pubmed_retrieval_reconciliation.tsv`, `pubmed_retrieval_reconciliation_by_query.tsv` and `pubmed_retrieval_reconciliation.json`.

## Screening status

The current deduplicated title/abstract queue contains 56,117 rows. **Zero records have been screened**, no full texts have been assessed, and no studies have been included. Manual licensed-database exports are still required. PRISMA values beyond identification and deduplication are currently zero because screening has not started; they are not completion claims.

## Targeted external follow-up

On 2026-09-23, a targeted publisher search revisited Huang et al. (2026),
PMID 41692760, DOI `10.1186/s12890-026-04178-2`, a bidirectional OSA–frailty
MR study. A source-level audit on 2026-09-24 confirms that the record was
already present in the frozen 2026-09-22 PubMed acquisition: it occurs in the
primary sleep–frailty, genetic sleep–frailty, and neighborhood query XML
sources. All three files match their recorded byte counts and SHA-256 values
in `pubmed_source_files.tsv`; the source query batches are dated 2026-09-22.
The PubMed history in each XML occurrence gives Entrez date 2026-02-15 and
PubMed date 2026-02-16, within the registered publication-date limits ending
2026-09-22. `all_records.tsv` contains the three source occurrences;
`deduplicated_records.tsv` and `screening/title_abstract_queue.tsv` each contain
one row, `primary_sleep_frailty:41692760:11f58fecd8d0`. The queue's two reviewer
decisions and adjudication remain blank. The record is therefore counted in
the existing PubMed snapshot and needs no supplemental PRISMA row or count
change. It remains unscreened and must pass the registered independent dual
title/abstract workflow. Automated priority labels are ordering aids only.

Evidence details and the correction to earlier absent-from-corpus notes are in
`analysis/pubmed_record_41692760_snapshot_verification_2026-09-24.md`.
The remaining licensed searches and dual screening are still outstanding.

## Supplemental PubMed window — 2026-09-23 through 2026-09-25

This later-indexed window is outside the frozen publication-date cutoff of
2026-09-22. Its three query summaries and two saved date segments reconcile
exactly to NCBI ESearch: 17 PubMed article nodes, zero book nodes, zero query
errors, and zero missing or extra segment-ID occurrences. The window has 17
unique records, no duplicate occurrences within the window, two PMID matches
to the frozen snapshot, and 15 PMIDs absent from that snapshot. The external
source increment remains separate and unmerged.

The reproducible builder
`scripts/55_build_pubmed_window_review_records.py` created a separate
17-row title/abstract queue and preserved all 17 source occurrences in
`reconciliation_2026-09-25_window/`. Its provenance records 17 unscreened
records, zero reviewer decisions, no mutation of the frozen main queue, and
no PRISMA-count change. The two snapshot matches and 15 later-only records
all remain unscreened; do not combine them with the frozen review or change
its PRISMA flow without a documented protocol/search amendment.

Key checksums:

- NCBI reconciliation report: `b1b550fa52863a2e3fabcb7b4c786f31f1dc31fec0de7cbbc9086834b093c9ff`.
- External window manifest: `eb565315eb117eecfd7694a2d78e5f849200fdd8f4ba26a167cc58b9a4c45e71`.
- Frozen `all_records.tsv` verified unchanged: `56b8b04079dd89afe1e19c8773f8830c0e9db74ca07ed8650870284170906887`.
- Supplemental output-manifest: `4c6d5f36f96188e84a5e6020ac5367e898c93d3aabd0f52fd0c77b2aa56511dd`.
- Supplemental occurrence table: `212c06eb73a0c4943f3d2056982076d0973df6b1e26817f266102504cba49c26`.
- Supplemental unique queue: `f938b5b694935088ae4ec7a3237624ea82e187506ff939a6b19bb643137f42b7`.

Run instructions and reviewer handling are in `README.md`. The unchanged
frozen queue remains at 56,117 records with zero decisions; the 17-row
supplemental queue is not counted as screened, excluded or included.

## Early reconciliation retry (2026-09-24 00:38 UTC; superseded)

The read-only `scripts/09_reconcile_pubmed_retrieval.py` retry could not resolve `eutils.ncbi.nlm.nih.gov` from this environment. All 15 ESearch requests failed with DNS `NameResolutionError`; no query or segment response was obtained, and failed output files were restored. At that time the retry left PMID 41692760's snapshot status unresolved. The later source-level verification documented above resolves that question from the archived acquisition files: the PMID is in the original query snapshot and already included once in the deduplicated queue. The failed live query adds no evidence beyond this archived record and does not require a count change.
