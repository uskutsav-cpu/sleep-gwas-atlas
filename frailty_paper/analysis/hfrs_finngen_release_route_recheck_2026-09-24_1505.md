# HFRS FinnGen release-route recheck — 2026-09-24 15:05 UTC

## Question

Does the current official FinnGen documentation identify a public, downloadable summary-statistics file for the custom R12 HFRS and HFRS-without-dementia analyses reported by Mak et al.?

## Evidence

The primary article describes a continuous HFRS constructed from 109 weighted ICD-10 codes, a FinnGen R12 analysis, a reported FinnGen sample size of 500,737, and a separate UK Biobank replication. The article says FinnGen results have a one-year embargo and that summary statistics are subsequently released to the research community; it does not provide an HFRS-specific accession, endpoint identifier, direct full-summary-statistics URL, or checksum ([Mak et al., *Nature Aging*, 2025](https://www.nature.com/articles/s43587-025-00925-y), Data availability and Methods).

The official FinnGen release handbook now lists `DF14/R13` with 500,186 participants in core GWAS, 4,662 total endpoints, 2,466 core endpoints, and public results expected around Q2 2026. It lists R12 results as public since Q2 2024. The handbook also says core GWAS summary statistics are released to FinnGen's green Google Cloud bucket about three months after a release ([FinnGen data freezes and releases](https://docs.finngen.fi/finngen-data-specifics/finngen-data-freezes-and-releases)).

The prior 2026-09-23 access-page recheck recorded that FinnGen's access page calls its 2026-06-02 public release `DF13`, describes 2,755 standard endpoints and approximately 500,186 participants, and directs bulk-summary users through an online form followed by emailed download instructions (`analysis/hfrs_finngen_access_recheck_2026-09-23.md`). The release naming and endpoint totals differ between that access page and the current handbook table; these labels cannot be treated as interchangeable without FinnGen's clarification. Neither source associates a standard endpoint with the custom weighted HFRS used in the article.

## Acquisition decision

- Exact R12 HFRS and HFRS-without-dementia full summary-statistics files remain **NOT ACQUIRED / NOT VERIFIED**.
- The release-wide documentation establishes when and where standard FinnGen core results are published; it does not establish that the paper's custom HFRS files are present, nor does it identify their filenames or phenotype IDs.
- No authenticated portal, form, email, author contact, or access-controlled endpoint was used in this recheck. No bulk listing was treated as available when the official web tool could not open the storage-listing API.
- Do not substitute the newer standard release or diagnosis endpoints for the R12 weighted-score analysis. Keep the registered HFRS analysis blocked until an ordinary official route or author-provided copy yields identifiable complete files and the build, alleles, effect scale, N, and checksum can be verified.

## Next permitted step

Use the existing unsent inquiry draft (`review/author_queries/hfrs_summary_statistics_request.md`) only when the project owner authorizes sending it, or obtain the files through the ordinary FinnGen access process. Ask FinnGen to clarify the `DF13` versus `DF14/R13` release labels and whether either HFRS file is included in the public core-results bucket. A release-wide public announcement alone does not answer that study-specific question.
