# PubMed date-window increment — 2026-09-25

## Acquisition

An incremental official PubMed E-utilities search covered publication dates **2026-09-23 through 2026-09-25 inclusive**. It used the exact three configured query strings and wrote to a separate resumable cache, leaving the full-history cache untouched.

| Query ID | ESearch/XML records |
|---|---:|
| `primary_sleep_frailty` | 1 |
| `genetic_sleep_frailty` | 0 |
| `frailty_neighborhood` | 16 |
| **Total query occurrences** | **17** |

The raw cache is `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/pubmed_cache_increments/2026-09-23_2026-09-25`. Its `window_manifest.json` preserves search terms, date bounds, per-query counts, config and collector hashes, and the explicit `SEPARATE_CACHE_NOT_MERGED` state. The source XML and ESearch ID lists remain there. SHA-256 for the manifest is `eb565315eb117eecfd7694a2d78e5f849200fdd8f4ba26a167cc58b9a4c45e71`.

## Validation and overlap

The existing reconciliation script re-queried NCBI ESearch at 2026-09-25 01:36:58 UTC. It checked the two non-empty date-window segments across all three query summaries: zero query errors, zero missing segment IDs, zero XML PMIDs outside ESearch, 17 PubMed article nodes, and zero book-article nodes. Detailed output is in `review/reconciliation_2026-09-25_window/`.

All 17 PubMed IDs in the increment match an XML article node, with no repeated ID within a query window and no overlap between the three query windows. Two IDs already appear in the prior full-history cache (`42675876`, `42781174`); the increment therefore contributes **15 previously unseen unique PMIDs** while retaining all 17 search occurrences for an auditable duplicate trail. PMID `42781174` has only year/month (`2026-Sep`) in its PubDate, but NCBI returned it for the exact requested publication-date window; this metadata limitation is retained, not inferred away.

## Screening and integration state

At the pre-increment checkpoint, the existing review builder had 61,009 PubMed source occurrences and a 56,117-record deduplicated title/abstract queue. A current read-only scan confirmed zero nonblank reviewer decisions, and `guard_review_outputs` passed. The 17 increment occurrences have not yet been merged into those cumulative outputs, so the existing queue counts remain the baseline and the new records have **no screening decisions**. The increment remains in its separate raw cache until it can be merged and the record builder safely regenerated.

The merge/rebuild is deferred while the active four-worker LAVA run is using the same external SSD and the host is under substantial memory pressure. At 2026-09-25 01:34 UTC, host-wide monitoring showed 76 MB unused physical RAM and 5,170 MB of 5,632 MB swap used. The collector's separate-cache guard and exact-window tests passed in the 122-test frailty suite. No search string, prior XML, reviewer output, or screening decision was overwritten.

## Candidate titles (not screening decisions)

The one primary-query hit is PMID `42776128`, “Physical Performance and Daytime Sleepiness as Correlates of Life-Space Mobility in Community-Dwelling Older Adults.” The 16 broad-neighborhood hits are `42781174`, `42779354`, `42779207`, `42779137`, `42778408`, `42778327`, `42778253`, `42778246`, `42778091`, `42776469`, `42776389`, `42775994`, `42775568`, `42773252`, `42675876`, and `41158073`. They remain unreviewed; query matching alone does not establish relevance or eligibility.
