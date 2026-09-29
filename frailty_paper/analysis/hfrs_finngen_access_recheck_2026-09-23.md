# HFRS and current FinnGen release access recheck

Audit date: 2026-09-23.

## What changed in the current access documentation

The HFRS GWAS paper reports its custom 109-weighted-ICD-10-code phenotype in FinnGen R12 (N=500,737), and identifies FinnGen's results route for its summary statistics. FinnGen's current access page now identifies Data Freeze 13 (DF13) as publicly released on 2026-06-02, with 2,755 standard endpoints and approximately 500,186 participants. The page says summary-statistic downloads are obtained by filling the online form; download instructions are sent by email. This clarifies that the current public FinnGen release is newer than the R12 release originally examined.

## HFRS availability decision

This release update does not establish that the paper's custom HFRS GWAS or its dementia-excluded sensitivity GWAS is one of the standard DF13 endpoints. The source paper analyzes an HFRS score it constructs from 109 weighted ICD-10 codes; no matching DF13 endpoint identifier, accession, or direct summary-statistics path was verified in this audit. Therefore:

- Keep the exact HFRS and HFRS-without-dementia results **NOT ACQUIRED / NOT VERIFIED**.
- Do not replace the R12 discovery results with DF13 or a similarly named diagnosis endpoint; release, phenotype construction, sample size and effect scale would need to match or be handled as a separately justified analysis.
- The previously identified `R18_SENILITY` endpoint is an ICD-10 diagnosis endpoint, not evidence of equivalence to the custom weighted score.
- The public summary-statistics route remains subject to an online access form and emailed download instructions. No form was submitted and no request was sent.

## Required next evidence

A source that can unblock the registered HFRS family must identify the exact HFRS (and, if used, HFRS-without-dementia) endpoint/file, release, phenotype definition, sample size/cases/controls, ancestry/build, effect and SE fields, and source checksum. Obtain it through FinnGen's ordinary authorized route or an author-provided release; do not infer equivalence from the general DF13 availability notice.

## Primary sources

- [Mak et al., Nature Aging (2025)](https://www.nature.com/articles/s43587-025-00925-y): describes the R12 custom HFRS analysis and data route.
- [FinnGen Access Results](https://www.finngen.fi/en/access_results): states DF13 public release date, current endpoint count, browser and form/email download procedure.
- [FinnGen DF12 release announcement](https://www.finngen.fi/en/results-based-full-finngen-cohort-500000-participants-released): documents the public R12 release.
