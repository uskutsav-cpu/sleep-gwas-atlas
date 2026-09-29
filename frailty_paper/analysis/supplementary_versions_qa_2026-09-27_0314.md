# Supplementary software/version table QA — 2026-09-27 03:14 UTC

- Rebuilt the source-linked supplementary package with pinned Python 3.11.11 using `scripts/39_build_supplementary_tables.py --repo . --replace-outputs` after refreshing Figure 4 provenance.
- S16 contains 52 rows and records `rsvg-convert` 2.62.3 for Figures 1, 2, and 4. The Figure 4 provenance is a direct builder input.
- The focused supplementary-table suite passed 7/7 tests (`analysis/test_supplementary_tables_py311_2026-09-27_0310.log`). The complete pinned Python 3.11.11 Frailty suite passed 169/169 (`analysis/test_suite_frailty_py311_2026-09-27_0314.log`).
- Provenance verification passed for all 29 recorded inputs and all nine output hashes. The S1 GWAS-manifest scope remains 33 rows; its row-scoped digest matches the builder provenance.
- S16 remains `PARTIAL_VERSION_INVENTORY`; this records a current validation run and corrects renderer usage evidence but does not claim a full scientific replay or complete tool-invocation coverage.
