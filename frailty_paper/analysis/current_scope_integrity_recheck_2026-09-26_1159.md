# Current frailty scope and input integrity recheck — 2026-09-26 11:59 UTC

Read-only validations from the repository root:

- `make -C frailty_paper PYTHON=.venv/bin/python verify-manifest FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'` — passed: 354/354 files matched manifest sizes and SHA-256.
- `make -C frailty_paper PYTHON=.venv/bin/python audit-manifest-metadata` — passed: 354 rows, all 22 required metadata columns.
- `make -C frailty_paper PYTHON=.venv/bin/python validate-plan` — passed: locked plan v1, 12 sleep traits, 396 tests.
- `make -C frailty_paper PYTHON=.venv/bin/python validate-overlap` — passed: 34-row cohort ledger and frozen 12-trait/source IDs match; exact participant intersections remain unknown.

These are structural, provenance, and contract checks. They do not resolve source eligibility, access restrictions, human screening, or cohort-level participant overlap. No inputs, thresholds, analysis settings, or review decisions were changed.
