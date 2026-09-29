# HFRS FinnGen public-source availability recheck — 2026-09-26

**Question:** Can the exact FinnGen HFRS full summary-statistics file used by Mak et al. be identified and downloaded with independently verifiable bytes?

## Findings

- The primary article is Mak et al., *Nature Aging* 5, 1589–1600 (2025), DOI [10.1038/s43587-025-00925-y](https://doi.org/10.1038/s43587-025-00925-y). It reports FinnGen HFRS GWAS N=500,737, a dementia-excluded HFRS sensitivity, and UK Biobank replication. Its data-availability section directs readers to general FinnGen access instructions; it does not give a HFRS-specific GWAS Catalog accession, release endpoint code, direct object path, or checksum.
- The official FinnGen [Access results page](https://www.finngen.fi/en/access_results), checked 2026-09-26, lists Data Freeze 13 as publicly released on 2026-06-02. For bulk summary-statistic downloads it instructs users to submit the official access form and receive download instructions by email. The page describes public availability of results generally but does not identify the paper's custom HFRS and HFRS-without-dementia outputs.
- The official [GWAS results format documentation](https://docs.finngen.fi/finngen-data-specifics/green-library-data-aggregate-data/core-analysis-results-files/gwas-results-format) describes the schema of standard FinnGen endpoint summary files, but does not establish that the custom continuous HFRS analyses are included in those endpoint files.

## Decision

The public release status has changed since earlier access checks, but **the exact HFRS summary-statistics object remains unverified**. This is not evidence that no file exists: an opaque-named or separately delivered custom result may exist. No new HFRS summary-statistic bytes were acquired, no checksum could be recorded, and no file was admitted to analysis. The published hit tables remain thresholded lead-variant evidence only. Do not substitute a diagnosis endpoint or infer that a generic FinnGen release contains the custom score results.

## Next verification step

After the authorized FinnGen download workflow provides its manifest/instructions, match the paper's exact phenotype (full HFRS or HFRS without dementia), cohort/release, and analytic model to the candidate object; then preserve the source URL/path, acquisition date, byte count, SHA-256, header/allele semantics, and independent file validation before considering harmonization.

**Access boundary:** This check used public pages only. No form was submitted and no access control was bypassed.
