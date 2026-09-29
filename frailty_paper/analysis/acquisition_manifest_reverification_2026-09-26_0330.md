# Acquired-resource manifest re-verification — 2026-09-26 03:29–03:30 UTC

The read-only resource integrity check passed against the explicit external data root.

- Manifest: `frailty_paper/manifests/all_acquired_resources.tsv`
- Manifest SHA-256: `dcab0d9c1b67fb81e35cb15b6dc0f951e7aa86662f850383290fb1277148555b`
- Command: `make -C frailty_paper verify-manifest PYTHON=python3 FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'`
- Result: `RESOURCE_MANIFEST_OK rows=354 files_verified=354` (exit status 0)
- Metadata command: `python3 frailty_paper/scripts/45_audit_acquisition_manifest_metadata.py --manifest frailty_paper/manifests/all_acquired_resources.tsv`
- Metadata result: `MANIFEST_METADATA_OK rows=354 required_columns=22` (exit status 0)
- Verifier SHA-256: `4be937ef1aca90c877d0aa93eb76649d5692b9f8545a0cdf58dc5875a84d51ed`
- Metadata-auditor SHA-256: `4c03f2c50098bd3e7507806800aa15c82a242ceaeb94387c0416ce9f0adeedb1`

The verification read and hashed the 354 files without modifying them. It establishes size/checksum agreement with the acquisition manifest and structural completeness of the 22 metadata columns only; it does not establish source eligibility, correct phenotype interpretation, cohort independence, or analysis readiness.
