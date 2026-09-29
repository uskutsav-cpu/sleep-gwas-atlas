# Integrated frailty preflight — 2026-09-26 03:49–03:52 UTC

Command: `make -C frailty_paper PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python preflight FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'` (pinned Python 3.11.11).

The run validated 354/354 manifest files and all 22 required metadata columns; locked plan (12 traits, multiplicity 396); 34-row overlap-ledger structure (exact intersections remain unknown); review structure (56,117 records, zero decisions and zero full-text records); bibliography (29 uses/27 entries, no key errors); missing-abstract sources (935/935 PMIDs found in verified XML, none with AbstractText; 53 duplicate source occurrences); quantitative manuscript claims (28/28); seven-factor Q_SNP audit (32,035,589 rows, zero qualifying variants/exclusions at the inclusive locked boundary); and reporting locators (92 rows, 93 locators, zero errors).

The first attempt exposed that the Makefile target did not forward an explicitly supplied `FRAILTY_EXTERNAL_STORAGE_ROOT`. The target now forwards it to the audit script; the rerun passed. The root allowlist behavior in the audit script is unchanged.

The suite ran 152 tests and ended with one import error: the untracked Brain6-only `test_audit_brain6_sleep_power_pilot` module imports `pytest`, unavailable in this pinned environment. The other 151 passed. This full preflight therefore exits nonzero; no Brain6 test or implementation was changed.

Scope: package/input validation only. It does not establish completion of human screening, licensed-search acquisition, sample-overlap resolution, LAVA, independent replication, downstream locus analyses, or manuscript readiness.
