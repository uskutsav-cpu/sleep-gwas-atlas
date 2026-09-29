# Supplementary software/version table QA — 2026-09-27 03:24 UTC

- Rebuilt the source-linked supplementary package with pinned Python 3.11.11 using `scripts/39_build_supplementary_tables.py --repo . --replace-outputs` after the 03:22 Figure 4 refresh.
- S16 contains 52 rows and records `rsvg-convert` 2.62.3 for Figures 1, 2, and 4. It also cites the 03:20 full 169/169 suite log.
- The full Frailty suite passed 169/169 at 03:20 UTC (`analysis/test_suite_frailty_py311_2026-09-27_0320.log`); the focused supplementary-table suite passed 7/7 after the final S16 rebuild (`analysis/test_supplementary_tables_py311_2026-09-27_0324.log`).
- Provenance validation passed for all 29 recorded inputs and all nine outputs. The S1 GWAS-manifest subset remains 33 rows and its row-scoped digest matches the builder provenance.
- S16 remains `PARTIAL_VERSION_INVENTORY`; this is software and artifact validation, not a completed scientific replay or full tool-invocation inventory.
