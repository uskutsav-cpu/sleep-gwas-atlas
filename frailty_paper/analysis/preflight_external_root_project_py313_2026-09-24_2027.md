# Integrated frailty preflight with external storage root

**Run window:** 2026-09-24 20:27–20:30 UTC
**Command:** `FRAILTY_EXTERNAL_STORAGE_ROOT="/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1" make -C frailty_paper PYTHON=.venv/bin/python preflight`
**Interpreter:** `frailty_paper/.venv/bin/python`, Python 3.13.13 (project virtual environment).
**Source revision at start:** `174e14ea5813ec00da9a54e5b21c0b71737cc525`.

The integrated preflight completed its manifest and structural gates: 351/351 acquired files verified; all 22 required manifest columns present; locked plan valid (12 sleep traits, 396 primary pairs); 34-row cohort ledger valid; review queue valid with 56,092 records and zero decisions; bibliography audit passed (29 citation uses, 27 entries); all 935 missing-abstract records mapped to source XML; manuscript numeric-claim audit passed 28/28; latent Q_SNP audit scanned 32,035,589 rows with zero significant variants and zero excluded rows; reporting locators passed 92 rows / 93 line references with zero errors; package suite passed 112/112 tests.

A preceding invocation without `FRAILTY_EXTERNAL_STORAGE_ROOT` failed the repository-root path check on 39 externally mounted paths. Supplying the configured external root allowed the verifier to check those same manifest entries by their declared size and SHA-256. No resource, frozen analysis input, screening decision, or scientific threshold was changed. The 112 passing tests are listed in the preflight output; the persistent preflight audit JSON outputs are under `analysis/`.

This validates the current frailty package against the mounted external resources. It does not validate unrelated Brain6/root-workflow changes, constitute a clean-checkout replay, resolve unresolved study eligibility or sample overlap, or complete manual database searches and screening.

Current hashes of persistent audit outputs:

- `manuscript_bibliography_audit_2026-09-23.json` — SHA-256 `1379683bada51cae999f7daf94e86085f16da15e82c4cff5121b37dc4fa8ed2c`.
- `latent_factor_qsnp_pruning_audit.provenance.json` — SHA-256 `2134124f0c26fb61c7dbc172a645210b2918ac43560d28e738bb647e10f99c34`.
- `manuscript_quantitative_claims_audit_2026-09-24.json` — SHA-256 `56838ca8bc8b3c04ffbe33587f0a340dd30afb4a0ab0de73c0c2124974bc3656`.
- `missing_abstract_source_audit_2026-09-24.json` — SHA-256 `cecc67de93ae4513f96b0d10a0ae07327969e69b7a205940e2c39959276a25e4`.
- `reporting_checklist_locator_audit_2026-09-23.json` — SHA-256 `9a2314b2e34e22e377a9ee34d58cf4cd500cbaa98d8287d9c9232f388b6c22ab`.
