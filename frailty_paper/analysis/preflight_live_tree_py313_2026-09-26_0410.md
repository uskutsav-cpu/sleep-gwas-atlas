# Integrated frailty preflight — Python 3.13 compatibility runtime

**Run window:** approximately 2026-09-26 04:10–04:15 UTC.
**Command:** `PYTHONPATH="$PWD/frailty_paper/.venv/lib/python3.13/site-packages" make -C frailty_paper PYTHON=python3 preflight FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'`
**Python:** 3.13.5 (recorded by the Q_SNP audit provenance).

All preflight stages completed through the final unittest target:

- Resource manifest: 354/354 files verified; 22/22 required metadata columns.
- Frozen plan: 12 sleep traits; multiplicity family 396.
- Cohort-overlap ledger: 34 rows structurally valid; exact intersections remain unknown.
- Review queue: 56,117 records; zero decisions, full-text records, or exclusions.
- Bibliography: 29 citation uses, 27 entries; no missing, duplicate, or uncited keys.
- Missing abstracts: 935/935 PMIDs found in retained verified XML; zero contain `AbstractText`; 53 duplicate source occurrences.
- Quantitative manuscript claims: 28/28 passed.
- Q_SNP: 32,035,589 rows across seven factors; zero qualifying variants and zero exclusions at the locked inclusive threshold.
- Reporting locators: 92 rows, 93 references, zero errors.
- Final unittest target: 151 tests passed (`Ran 151 tests ... OK`). A separate pytest collection recorded in `analysis/test_suite_py313_compat_2026-09-26_0404.md` passed 154 tests and 92 subtests.

The shell session completed after the final `OK` test output. This mixed live-tree check is not a clean-checkout attestation and does not resolve human screening, licensed-source access, exact cohort overlap, source eligibility, LAVA completion/QC, or downstream gated analyses. The standard pinned Python 3.11 environment still lacks pytest; this run used the installed Python 3.13 runtime and project dependency directory without altering either environment.
