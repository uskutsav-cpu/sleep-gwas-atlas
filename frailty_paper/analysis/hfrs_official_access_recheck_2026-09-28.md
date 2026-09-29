# HFRS and FinnGen official-access recheck — 2026-09-28

## Question

Does the current primary-source documentation identify the full FinnGen R12 summary-statistics object for the custom Hospital Frailty Risk Score (HFRS) and HFRS-without-dementia GWAS, and can the public standard-release catalog substitute for it?

## Findings

- Mak et al. report a continuous score formed by summing 109 weighted ICD-10 codes, FinnGen R12 discovery at N=500,737, UK Biobank replication at N=407,463, and 53 independent lead variants for the primary HFRS analysis. The article's data-availability statement directs readers to FinnGen's results-access page and describes the one-year embargo/release cadence; it does not give a custom endpoint identifier, direct full-summary-statistics URL, object name, checksum, or complete genome-wide file.
- FinnGen's official Access Results page, checked 2026-09-28, states that DF13 was publicly released on 2026-06-02 with 2,755 endpoints and 500,186 participants. It says summary statistics are freely downloadable from cloud storage after completing the online access form; download instructions are then sent by email. The public endpoint browser supports discovery, but its release-wide availability does not identify the custom HFRS score.
- The previously acquired publisher supplement remains thresholded results, not full summary statistics: it includes reported FinnGen variants meeting P<5e-8 and UK Biobank/meta-analysis fields. It is insufficient for the planned genome-wide LDSC/LAVA analyses.
- The existing literal ICD-10 crosswalk to DF13 found 99 of the 109 HFRS code rows as exact code tokens in standard endpoint definitions and 10 without literal matches. Neither those component endpoints nor their aggregate reconstructs the paper's weighted person-level HFRS phenotype. No proxy endpoint is justified.

## Decision and action

The exact custom FinnGen HFRS and HFRS-without-dementia full-statistics objects remain **unverified**. Keep HFRS unavailable for the locked pairwise analysis until the object identity, release, phenotype definition, analytic sample, model, genome build, effect scale, checksum, overlap and reuse terms are established. Do not submit the FinnGen access form or treat standard ICD endpoints as substitutes. No form was submitted and no file was downloaded in this recheck.

## Source records

- Mak et al., *Nature Aging* (2025), [article](https://www.nature.com/articles/s43587-025-00925-y), including methods and data-availability statements.
- FinnGen, [Access Results](https://www.finngen.fi/en/access_results), including current DF13 release facts and bulk-download instructions.
- FinnGen, [DF13 release announcement](https://www.finngen.fi/en/finngen-data-freeze-13-results-now-available).
- Local exact-code crosswalk: `analysis/hfrs_df13_109code_crosswalk_2026-09-27.md` and its 109-row TSV.
- Publisher supplementary workbook/source audit: `analysis/hfrs_published_supplement_audit_2026-09-24.md`.

## Provenance

Read-only web verification on 2026-09-28. No network download, access-form submission, data transformation, analysis, or endpoint substitution was performed. The claim about the custom file's absence is bounded to the inspected article/supplement and official public access/catalog descriptions; opaque or separately delivered files are not ruled out.
