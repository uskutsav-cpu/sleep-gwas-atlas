# Phase 1 acquisition-manifest re-verification

**Rechecked:** 2026-09-23 22:18 UTC
**Repository:** `frailty-paper-v1`, live `HEAD` `9025a20e6925bcf0460d9adb75fc248d233de850`

## Current evidence

From the repository root, the pinned Python 3.11.11 verifier was run with the documented external storage root:

```sh
FRAILTY_EXTERNAL_STORAGE_ROOT='/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1' \
  work/conda-envs/frailty-paper-py311-clean-2026-09-23/bin/python \
  frailty_paper/scripts/16_verify_acquired_resource_manifest.py \
  --repo . --external-storage-root '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'
```

Observed result: `RESOURCE_MANIFEST_OK rows=347 files_verified=347`. The manifest contains 22 columns and 347 rows (32 GWAS summary-statistic paths, 311 PubMed XML batches, and four public resource-index records); every manifest cell is populated, with unresolved and non-applicable values represented explicitly. The acquired-resource scan-coverage audit reports all 32 manifested GWAS paths in the source-scan report union. Its scope is path coverage; it does not re-run every row-level QC or establish scientific eligibility.

The locked plan was independently revalidated with `scripts/13_validate_analysis_plan.py`: `ANALYSIS_PLAN_LOCK_OK version=1 sleep_traits=12 multiplicity=396`.

## Collection-log reconciliation

Historical `logs/full_collection.log` and `logs/finish_collection.log` record interrupted transfers and an earlier weight-loss download attempt that failed expected size/MD5 checks. The current registered weight-loss file is present at the expected 242,305,000 bytes, and the complete manifest verifier passed its recorded SHA-256. No partial file is currently present beside it. This resolves the current file-integrity question for the registered copy; it does not erase the failed-attempt history or resolve the physical-component source/provenance blocker. No collection job was restarted.

## Capacity and remaining limits

At this recheck, internal free space was 4.8 GiB and the configured external SSD had 1.7 TiB free. Large acquisition, harmonization and temporary outputs must continue to use the external workspace. Fried full summary statistics, the exact custom HFRS endpoint, manual database exports, and several source-eligibility facts remain unresolved as documented in `BLOCKERS.md`; the 347-file manifest is not evidence that those absent resources were acquired.
