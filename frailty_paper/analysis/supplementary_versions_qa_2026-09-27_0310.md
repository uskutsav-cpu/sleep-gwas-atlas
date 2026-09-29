# Supplementary software/version table QA — 2026-09-27 03:10 UTC

- Rebuilt the source-linked supplementary package with pinned Python 3.11.11 using `scripts/39_build_supplementary_tables.py --repo . --replace-outputs`.
- S16 remains 52 rows and explicitly records `rsvg-convert` 2.62.3 for Figures 1, 2, and 4. Figure 4 renderer provenance is now an input to the builder rather than inferred from the earlier figure records alone.
- The focused supplementary-table suite passed 7/7 tests (`analysis/test_supplementary_tables_py311_2026-09-27_0310.log`).
- Provenance verification passed for all 29 recorded inputs and all nine output hashes. The S1 GWAS-manifest scope remains 33 rows; its row-scoped digest matches the builder provenance.
- S16 remains `PARTIAL_VERSION_INVENTORY`; this corrects renderer usage evidence but does not claim full replay or complete tool-invocation coverage.
