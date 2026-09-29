# HFRS summary-statistics access recheck — 2026-09-29

## Scope

Recheck whether the official publication and FinnGen access documentation identify a directly downloadable, full summary-statistics object for the custom Hospital Frailty Risk Score (HFRS) and HFRS-without-dementia GWAS reported by Mak et al.

## Primary-source findings

- The publisher's 2025 article describes a custom continuous HFRS GWAS in FinnGen (N=500,737), UK Biobank replication (N=407,463), and a sensitivity analysis removing dementia weights. It reports 53 independent lead variants for the main HFRS analysis. The article's Data availability section says FinnGen summary statistics are made available after a one-year consortium embargo through FinnGen's results-access route. It does not give an HFRS-specific endpoint code, accession, direct full-statistics URL, object name, or checksum.
- FinnGen's current official access-results page lists public R13 browsing and explicitly says bulk summary statistics require an online form, after which download instructions are emailed. The official R4 download documentation describes the same access mechanism. Public browsing and licensed/registered bulk download are distinct routes.
- The one-year embargo has elapsed as of this recheck. That fact alone does not establish that the exact custom HFRS files are in the standard DF13/R13 endpoint bundle. The prior DF13 endpoint scan and 109-code crosswalk found no literal HFRS/frailty endpoint label; matching 99 component ICD-10 categories is not evidence that the custom weighted score is the same phenotype.

## Decision

The exact HFRS and HFRS-without-dementia summary-statistics objects remain unidentified and ineligible for harmonization or genetic analysis. Do not substitute an ordinary FinnGen endpoint, the paper's thresholded lead-variant tables, or the UK Biobank replication results. No form was submitted and no data were downloaded during this recheck.

## Required next evidence

Obtain the official FinnGen bulk manifest/instructions through the authorized registration route, then establish the exact phenotype definition, endpoint/object identifier, cohort and analytic sample, file checksum, schema, build, effect allele convention, and complete genome-wide coverage before considering use. If neither custom object is listed, request the exact full summary statistics from the study authors through an approved channel.

## Sources

- Mak et al. (2025), [Nature Aging article](https://www.nature.com/articles/s43587-025-00925-y), including sample sizes, lead-variant counts, and Data availability statement.
- [FinnGen Access results](https://www.finngen.fi/en/access_results).
- [FinnGen R4 data download documentation](https://finngen.gitbook.io/documentation/r4/data-download).
- Earlier route evidence: `hfrs_official_route_recheck_2026-09-27.md`, `hfrs_df13_public_endpoint_index_scan_2026-09-27.md`, and `hfrs_df13_109code_crosswalk_2026-09-27.md`.
