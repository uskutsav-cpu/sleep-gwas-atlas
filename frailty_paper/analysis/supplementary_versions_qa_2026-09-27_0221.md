# Supplementary software/version table QA — 2026-09-27 02:21 UTC

- Rebuilt with pinned Python 3.11.11 using `scripts/39_build_supplementary_tables.py --repo . --replace-outputs`.
- S16 contains 52 rows, including the latest 167/167 package-suite record and current manuscript audit records. Its status remains `PARTIAL_VERSION_INVENTORY`: this table does not claim an exhaustive record of every tool invocation or a complete scientific replay.
- The seven focused supplementary-table tests pass after updating the latest suite evidence assertion.
- Supplementary provenance records 28 inputs and nine generated output hashes. The GWAS resource-manifest input is intentionally scoped to the 33 GWAS rows used by S1; its recorded hash is the canonical hash over those consumed rows, not the whole mixed-worktree manifest. The other recorded inputs and all outputs match their recorded SHA-256 values.
- Full pinned Python 3.11.11 Frailty suite: 167/167 passed; see `test_suite_frailty_py311_2026-09-27_0226.log`.
- Limits remain: the inventory is partial, manuscript audits cover their defined checks only, and the scientific analysis/review is unfinished.
