# Current acquisition-manifest and GWAS scan coverage audit

Date: 2026-09-23

Branch: `frailty-paper-v1`

## File integrity

Ran the existing read-only acquired-resource verifier from the repository root with the configured external storage location explicitly allowed:

```text
python frailty_paper/scripts/16_verify_acquired_resource_manifest.py --repo . --external-storage-root "/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1"
RESOURCE_MANIFEST_OK rows=347 files_verified=347
```

This verifies the manifest's required identity/path/byte-count/SHA-256 fields against each currently reachable file. It confirms the previously registered sizes and hashes; it does not re-establish scientific validity, licenses, phenotype equivalence, or source provenance.

## Source-scan coverage

A read-only path-coverage comparison of `all_acquired_resources.tsv` against the source scan indexes found:

- 32 unique manifested GWAS-summary-statistic paths.
- All 32 paths are represented in source scan reports: 20 in `gwas_source_scan.tsv` and 12 in `gwas_source_scan_sleep_panel.tsv`.
- Specialized reports overlap those 32 paths and document latent-factor, Parkinson, and newly acquired-source checks; none revealed an additional unmanifested path in this comparison.
- No manifested GWAS path was missing from the union of scan reports.

The check compared path coverage only; it did not reparse every source row or independently reproduce each scan's QC. The latent-factor scan's separate basic, duplicate-key, allele-orientation, and HapMap3 audits remain bounded to the scopes documented in those reports.

## Manifest schema completeness

The acquisition manifest has all 22 requested columns. Every one of the 347 rows has a nonempty value in each column. Non-applicable or unresolved metadata are represented explicitly (for example, `NA` or `UNKNOWN`) rather than as empty cells. Explicit unknown values are still unresolved scientific metadata; populated cells do not mean those fields have been verified.

## Current interpretation

The acquisition manifest remains complete at the registered file/hash level. Acquisition is not equivalent to analysis eligibility: Fried full summary statistics remain unavailable; physical-component source build/effect/model provenance remains unresolved; the exact HFRS endpoint is unverified; exact cohort intersections remain unknown; downstream frailty local and molecular analyses remain gated. These blockers are tracked in `BLOCKERS.md` and the source-specific audits.
