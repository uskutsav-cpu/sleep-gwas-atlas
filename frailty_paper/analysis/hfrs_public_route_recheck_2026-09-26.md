# HFRS public-results route recheck — 2026-09-26

## Question

Does the current FinnGen public-release update identify or directly expose the exact custom HFRS and HFRS-without-dementia summary-statistics files used by Mak et al.?

## Findings

The current FinnGen access page says Data Freeze 13 (DF13) was publicly released on 2026-06-02, includes 2,755 disease endpoints and can be browsed through R13. Bulk summary-statistic downloads remain available through an access form; FinnGen emails the download instructions afterward. The page describes the results as free to download once that route is completed. It does not provide an HFRS-specific accession, endpoint code or direct file URL.

The current FinnGen clinical-endpoints page links a DF13 endpoint-definition workbook and lists DF14 as the current definitions page. The linked DF13 workbook URL returned HTTP 406 to a direct request from this environment. The workbook was not acquired, and that response does not show that the endpoint or statistics are absent. No access form was submitted and no access-control workaround was attempted.

A separate official FinnGen Handbook route documents user-generated custom GWAS outputs: users can search the `userresults.finngen.fi` browser by the analysis title assigned in the Sandbox, while custom result files are stored under a release-specific `sandbox_custom_gwas/<phenotype_name>` location in the green library. The browser workflow requires signing in with a FinnGen user account. This is a potentially useful exact-name route if the paper authors disclose their analysis title/phenotype name, but it is not evidence that the paper's SAIGE outputs were published there; the paper and accessible public pages still do not give that name or a file identifier.

Mak et al. was published online 2025-08-05, so the stated one-year FinnGen embargo has elapsed. The article describes custom SAIGE analyses of continuous HFRS and HFRS without dementia in FinnGen, with N=500,737. That article-specific sample and analysis identity do not match automatically to a standard DF13 core endpoint or to a later release with a different total participant count. The publisher workbook remains thresholded reported-hit tables, not genome-wide summary statistics. The prior standard R12 manifest scan remains a bounded label search and does not exclude an opaque-named or separately delivered custom file.

## Decision and next step

The embargo is no longer the timing blocker. Exact custom-file identity, endpoint, checksum and usable all-variant content are still unverified. Keep both HFRS files ineligible for harmonization, LDSC and LAVA; do not substitute a standard DF13 endpoint. Continue only through the official results-access form/email route or source-specific documentation that identifies the exact custom files. A targeted R13 browser lookup remains deferred until an eligible locus exists under the frozen analysis plan.

## Official sources

- FinnGen, [Access results](https://www.finngen.fi/en/access_results): DF13 release date, sample, endpoint count, R13 browser and bulk access procedure.
- FinnGen, [Clinical endpoints](https://www.finngen.fi/en/researchers/clinical-endpoints): DF13/DF14 endpoint-definition resources; linked DF13 workbook URL: `https://www.finngen.fi/sites/default/files/inline-files/FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx`.
- FinnGen Handbook, [Custom GWAS GUI tool](https://docs.finngen.fi/working-in-the-sandbox/which-tools-are-available/untitled/custom-gwas-tool): analysis-title lookup in `userresults.finngen.fi`, sign-in requirement and release-specific custom-GWAS green-library location.
- Mak et al., [Nature Aging article](https://www.nature.com/articles/s43587-025-00925-y): publication date, custom HFRS analyses, discovery sample and data availability.
- Previous bounded HFRS checks: `hfrs_r12_public_manifest_recheck_2026-09-24.md`, `hfrs_df13_release_route_recheck_2026-09-24.md`, and `hfrs_public_route_audit_2026-09-25.md`.
