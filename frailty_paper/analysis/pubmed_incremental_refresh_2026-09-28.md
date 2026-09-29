# PubMed incremental refresh — 2026-09-28

## Result

Ran the repository's frozen three-query PubMed configuration for publication dates 2026-09-28 through 2026-09-28. At 2026-09-28 19:53 UTC, NCBI ESearch returned zero records for each query: primary sleep–frailty, genetic sleep–frailty, and frailty-neighborhood. The collector retained the exact query strings and returned counts in the separate source cache. The subsequent ESearch reconciliation completed with zero query errors, zero missing/extra IDs, and zero PubMed article/book nodes. There were no XML source records, no unique supplemental review records, and no reviewer packets to prepare.

The main 56,117-record snapshot and its title/abstract queue remain untouched. The canonical queue SHA-256 is `3ae3fb3b4765c85536c3ccf64d20297a632af892e2c66e622ee665db1ea43b2e`; the builder provenance confirms `main_screening_queue_mutated=false`, `prisma_counts_changed=false`, and zero decisions created. This is a date-bounded query snapshot as of the retrieval time; it does not establish that no additional records will later be indexed for this publication date.

## Reproducibility

- Query configuration: `config/pubmed_queries.json`, SHA-256 `057a46bbff6af9c720a35bcb285622bbfda043477667920293629aafed5a368d`.
- Collector: `scripts/01_pubmed_search.py`, SHA-256 `e0a9d6f314a68cfe79921a9fcbaad50b80848b953ec89c03c8160f03907dd477`.
- Bounded-window command: `frailty_paper/scripts/54_collect_pubmed_window.py --start-date 2026-09-28 --end-date 2026-09-28 --outdir frailty_paper/review/incremental_windows/2026-09-28_2026-09-28/source`.
- Record builder: `frailty_paper/scripts/55_build_pubmed_window_review_records.py`; it checked the frozen record-build manifest and generated empty, header-only supplemental tables.
- Separate source files: `review/incremental_windows/2026-09-28_2026-09-28/source/`; window-manifest SHA-256 `335e8155a5c89e5ea1b6f0992b4dfaa3510dd82d7eabed97ddd38e9b72042943`.
- ESearch reconciliation: `review/incremental_windows/2026-09-28_2026-09-28/reconciliation_2026-09-28/`; summary SHA-256 `476e5e7d9d5a6282542eb147cf46b779b8588997e02b40746aac8ba227bd06dc`.
- Supplemental review build: `review/incremental_windows/2026-09-28_2026-09-28/review/`; provenance SHA-256 `3803094353b4d9cf8adc1376143be1822c1a3c60353dc3fe975e13035e0af2b5`.

No source access controls were bypassed. The initial shell attempt lacked the `requests` dependency and the first retry encountered DNS failure; after using the repository venv and obtaining turn-scoped network access, the official NCBI collector and reconciliation completed successfully.
