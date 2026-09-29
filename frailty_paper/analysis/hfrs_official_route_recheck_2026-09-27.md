# HFRS full-statistics route recheck — 2026-09-27

## Question

Has the official publisher or FinnGen route now exposed an identifiable, directly downloadable full summary-statistics object for the custom HFRS and HFRS-without-dementia analyses in Mak et al.?

## Primary-source findings

- The [Mak et al. Nature Aging article](https://www.nature.com/articles/s43587-025-00925-y) reports continuous custom HFRS GWAS in FinnGen (N=500,737), UK Biobank variant replication (N=407,463), and an HFRS-without-dementia sensitivity analysis. Its Data availability section links generally to FinnGen access instructions; it does not provide a HFRS-specific endpoint code, accession, direct file URL, object name, or checksum.
- The current [FinnGen Access results page](https://www.finngen.fi/en/access_results) states that DF13 was publicly released on 2026-06-02, lists 2,755 standard endpoints and 500,186 participants, and says bulk summary-statistics access starts with an online form; download instructions are then emailed. The public browser can be searched without that download step. The general DF13 listing does not identify either custom HFRS definition.
- The linked [official FinnGen download registration form](https://elomake.helsinki.fi/lomakkeet/124935/lomake.html) confirms that the form asks for investigator identity, organization, email and a use/compliance certification, then emails download directions. The form was inspected but not submitted.
- The FinnGen clinical-endpoint list linked from the results page timed out during this check. That failed fetch is not treated as evidence that no endpoint exists.

## Decision

The exact full HFRS and HFRS-without-dementia files remain unidentified and ineligible for harmonization, LDSC or LAVA. DF13 availability does not establish that the custom HFRS GWAS is included among the standard endpoints. The publisher workbook remains thresholded reported-hit tables only. No substitute endpoint was used, no external form was submitted, and no new summary-statistic bytes were acquired.

## Next step

Use the official FinnGen form with investigator/contact details and certification, then use the emailed instructions and manifest to search by the exact paper phenotype, definition, cohort and analysis. Preserve object identity, byte count, SHA-256, schema, build, allele semantics and cohort details before considering analysis.
