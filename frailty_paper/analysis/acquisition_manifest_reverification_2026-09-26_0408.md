# Acquisition manifest integrity re-verification — 2026-09-26 04:08 UTC

Commands:

- `make -C frailty_paper PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python verify-manifest FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'`
- `make -C frailty_paper PYTHON=../work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python audit-manifest-metadata`

Observed results: `RESOURCE_MANIFEST_OK rows=354 files_verified=354` and `MANIFEST_METADATA_OK rows=354 required_columns=22`.

The manifest TSV is 292,025 bytes with SHA-256 `dcab0d9c1b67fb81e35cb15b6dc0f951e7aa86662f850383290fb1277148555b`. The configured external root supplied all file bytes. This verifies registered path, size and checksum integrity and metadata structure; it does not establish source eligibility, licensing, cohort independence, or analysis readiness.

At 04:08:26 UTC the external LAVA runner state reported 9,663/29,940 receipts, two effective workers with distinct active claims (sleep_apnea locus 553; sleep_timing locus 978), zero launch failures and stale recoveries, and unchanged analysis-lock SHA-256 `74d0df6e897ff75ca7128876438cfc9bcb3495f4033cb35f77bef4cee4a00e70`. This is a runner-state snapshot, not a receipt-file re-audit.
