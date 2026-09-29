# Supplementary software/version table QA — 2026-09-27 02:39 UTC

- Rebuilt with pinned Python 3.11.11 using `scripts/39_build_supplementary_tables.py --repo . --replace-outputs`.
- S16 contains 52 rows, including the latest 169/169 package-suite record and current manuscript audit records. It remains `PARTIAL_VERSION_INVENTORY`; it does not claim every tool invocation or a complete scientific replay.
- Seven focused supplementary-table tests pass. The full pinned Frailty suite passes 169/169 (`test_suite_frailty_py311_2026-09-27_0239.log`), including the PubMed ESearch transient-error handling tests.
- Supplementary provenance records 28 inputs and nine generated output hashes. All non-scoped input/output hashes match. The GWAS resource-manifest input is intentionally scoped to the 33 rows consumed by S1 and its scoped hash/count match the recorded values.
- Current manuscript audits pass: bibliography 29 citation uses/27 entries, quantitative claims 28/28, and reporting locators 92 rows/93 references.
- Limits remain: the inventory is partial, these validation audits cover defined checks only, and the scientific analysis/review is unfinished.
