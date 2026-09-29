# HFRS public-results route audit — 2026-09-25

## Question

Has the one-year FinnGen embargo expired, and does the currently documented public route establish that the exact Mak et al. R12 HFRS and HFRS-without-dementia summary-statistics files are available for immediate download?

## Findings

Mak et al. was published online on 2025-08-05. The article reports FinnGen HFRS discovery GWAS of 500,737 participants and states that FinnGen summary statistics become available to the scientific community after a one-year embargo. Its FinnGen methods identify the HFRS and HFRS-without-dementia analyses as custom SAIGE v0.35.8.8 analyses, adjusted for birth year, sex and ten PCs. As of this audit date, the stated embargo period has elapsed.

The current FinnGen results page lists DF13 as the latest public core release (released 2026-06-02) and R12 as the earlier release. It says to request bulk summary-statistic access using FinnGen's online form and that download instructions arrive by email. It also lists the R13 core PheWeb browser. The earlier search of the standard R12 core manifest returned 2,469 endpoint rows and no metadata labels matching frailty, HFRS or hospital frailty; that is a bounded core-manifest name search and does not search separately delivered custom analyses.

FinnGen's custom-GWAS fine-mapping documentation says results from its unmodifiable pipeline are placed in the userresults browser and endpoint-specific green-library folder (`sandbox_custom_gwas/<endpoint>/finemap`). That describes fine-mapping outputs for endpoints run through that pipeline. It does not identify an HFRS endpoint or prove that Mak et al.'s custom SAIGE files were placed there. The accessible article and provider pages reviewed here still do not identify an exact endpoint code, accession, direct full-statistics URL, or checksum for either HFRS definition.

## Decision

The one-year embargo is no longer the timing blocker. Exact-file identity and access remain unresolved. Keep the HFRS full-statistics and HFRS-without-dementia rows ineligible for harmonization, LDSC or LAVA; the publisher workbook remains thresholded reported-hit tables only. Do not treat standard core R12 endpoint results or later DF13 results as substitutes for the custom SAIGE analyses. The remaining official bulk-download path is the provider's results-access form followed by emailed instructions. No form was submitted, and no author or provider message was sent.

## Sources checked

- Mak et al., *Nature Aging* (published 2025-08-05), [primary article](https://www.nature.com/articles/s43587-025-00925-y), including Methods and Data availability.
- FinnGen, [Access results](https://www.finngen.fi/en/access_results), current DF13/R12 release and bulk-download instructions.
- FinnGen Handbook, [Fine-mapping of custom GWAS analyses](https://docs.finngen.fi/working-in-the-sandbox/which-tools-are-available/untitled/finemapping-of-custom-gwas-analyses), describing the public browser and endpoint-specific green-library route for eligible unmodifiable-pipeline fine-mapping results.
- Existing bounded standard-release scan: `analysis/hfrs_r12_public_manifest_recheck_2026-09-24.md`.
