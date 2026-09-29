# Phase 1 acquisition manifest recheck — 2026-09-26 08:55 UTC

Revalidated the current 354-row `manifests/all_acquired_resources.tsv` against the external storage root.

- Metadata audit: `MANIFEST_METADATA_OK rows=354 required_columns=22`.
- File audit: `RESOURCE_MANIFEST_OK rows=354 files_verified=354`; every manifest file was present and matched its recorded byte size and SHA-256.
- Command: `python3 frailty_paper/scripts/45_audit_acquisition_manifest_metadata.py --manifest frailty_paper/manifests/all_acquired_resources.tsv`, followed by `python3 frailty_paper/scripts/16_verify_acquired_resource_manifest.py --repo . --manifest frailty_paper/manifests/all_acquired_resources.tsv --external-storage-root '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'`.
- Disk at recheck: internal volume 6.6 GiB free; external SSD 1.6 TiB free. The internal volume remains below the 20 GiB analysis-space gate. Keep derived outputs and temporary work on the external SSD and rerun the destination-specific capacity gate before large jobs.

This verifies manifest integrity and file identity only. It does not establish that every acquisition has sufficient reuse licensing, source provenance, build metadata, or analysis eligibility. Those limitations remain recorded in the resource manifest, source audits, and `BLOCKERS.md`.
