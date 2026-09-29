# Phase 1 acquisition-manifest recheck — 2026-09-25 08:40 UTC

## Current integrity result

Rechecked the current registered acquisition inventory without restarting collection.

- `45_audit_acquisition_manifest_metadata.py`: PASS, 354 rows and all 22 required columns.
- `16_verify_acquired_resource_manifest.py --external-storage-root '/Volumes/Extreme SSD'`: PASS, all 354 registered file paths, recorded sizes and SHA-256 values verified; zero missing files or size/hash mismatches.
- Manifest SHA-256: `dcab0d9c1b67fb81e35cb15b6dc0f951e7aa86662f850383290fb1277148555b`.
- Registered file bytes total 14,770,326,676 (13.76 GiB): 32 GWAS files, 311 PubMed XML batches, five public resource indexes, one donor-metadata record, one library-metadata record, one publisher hit-table workbook, one sample-metadata workbook and two GEO MINiML archives.
- All 354 paths were present and had the registered byte counts before hashing. `bash -n frailty_paper/scripts/07_run_all_collection.sh` passed.

The external data directory remains mounted at `/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/frailty_paper/data`; the PubMed cache is mounted under the external analysis workspace. At the 08:38 UTC disk sample, the internal filesystem had 2.9 GiB free and the external SSD had 1.6 TiB free. Large analysis outputs and scratch must use the external workspace; the internal volume remains below the 20 GiB analysis gate.

No full-collection or PubMed acquisition process was present in process inspection. `logs/full_collection.log` is still a historical 2026-09-22 file; the current collection script is syntactically valid, so that old logged error does not establish a current script defect. No downloads were restarted.

## Scope and unresolved Phase 1 eligibility

This recheck establishes byte integrity against the current manifest, not source eligibility, licensing, genome-build truth, cohort overlap or harmonization readiness. Existing source-specific audits remain authoritative for those questions. The Fried score file is still unavailable through its ordinary challenge-blocked route; the five physical-component GWAS lack source-linked build/effect/provenance; the HFRS supplementary hit table is not a full summary-statistics resource; licensed database exports and human screening remain outstanding; and functional resource access/identity blockers remain documented in `BLOCKERS.md`. The acquired 354-row manifest must not be read as evidence that every desired resource was acquired or is eligible for analysis.
