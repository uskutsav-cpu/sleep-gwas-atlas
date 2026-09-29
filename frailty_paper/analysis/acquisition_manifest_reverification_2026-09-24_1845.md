# Acquisition manifest reverification — 2026-09-24 18:45 UTC

This read-only checkpoint revalidated the current Phase 1 resource manifest
against the mounted external storage root. No source files were downloaded or
modified.

## Results

- Resource manifest: 351 rows; every listed file exists and matches its
  recorded byte count and SHA-256 (`RESOURCE_MANIFEST_OK rows=351
  files_verified=351`).
- Metadata audit: all 351 rows pass all 22 required metadata columns
  (`MANIFEST_METADATA_OK rows=351 required_columns=22`).
- Manifest SHA-256: `cb7e26c1192ae73c52acf45ce69f79a8a8aff928ab365fa919ecc8cb71894685`.
- Internal filesystem: 5.2 GiB available; external SSD: 1.7 TiB available.

## Reproduction

From the repository root:

```sh
python3 frailty_paper/scripts/16_verify_acquired_resource_manifest.py \
  --repo . \
  --external-storage-root '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'
python3 frailty_paper/scripts/45_audit_acquisition_manifest_metadata.py \
  --manifest frailty_paper/manifests/all_acquired_resources.tsv
```

Both scripts use the Python standard library and are read-only. The mounted
resource path resolves to `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data`.

## Scope limits

This establishes current byte integrity and metadata structure only. It does
not resolve phenotype eligibility, genome-build/effect semantics, study
provenance gaps, cohort overlap, licenses, or manual-access blockers. The
working tree still reports the established `frailty_paper/data/` symlink as
untracked and three Git-tracked files beneath that path as deleted; the
read-only verifier resolved the mounted external data path and verified all
351 manifest entries. These Git path entries are not included as manifest
resource rows.
