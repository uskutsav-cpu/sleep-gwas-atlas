# HFRS summary-statistics route recheck after DF13 release

**Checked:** 2026-09-24 20:03 UTC. **Result:** the one-year embargo for the 2025 paper has elapsed, and FinnGen DF13 is publicly released, but the exact custom HFRS summary-statistics files remain unverified and were not downloaded.

## Evidence

- Mak et al. published the HFRS study on 2025-08-05. It used FinnGen R12 and a continuous score assembled from 109 weighted ICD-10 codes, with a separate HFRS-without-dementia sensitivity definition. The paper reports N=500,737 for the FinnGen analysis. Its Data Availability section says FinnGen summary statistics become public after a one-year embargo and points to FinnGen's results page; it does not provide an HFRS-specific accession, file URL, or checksum.
- FinnGen's current [Access results page](https://www.finngen.fi/en/access_results) identifies DF13 as publicly released on 2026-06-02, with 2,755 standard disease endpoints and 500,186 participants. The results are browsable in the [R13 browser](https://r13.finngen.fi/). The same access page still requires an online form to receive emailed instructions for summary-statistic downloads.
- The prior official R12 manifest scan searched all 2,469 manifest rows and found no HFRS/frailty label match. It was a standard-manifest metadata scan, not an exhaustive scan of every result file or a search for opaque file names (`hfrs_r12_public_manifest_recheck_2026-09-24.md`).

## Interpretation and next step

The elapsed embargo and release of DF13 do not establish that the study's custom R12 composite is included in a standard release. DF13 disease-endpoint results are not interchangeable with a continuous weighted sum of 109 codes, nor with the paper's separately defined no-dementia score. The article's N also differs from the R12 core-release sample count, so do not infer matching participants, phenotype construction, effect scale, build, or alleles from the release label alone.

No access form was submitted, no email or author message was sent, and no file was downloaded. The exact FinnGen R12 HFRS and HFRS-without-dementia genome-wide files remain **NOT ACQUIRED / NOT VERIFIED**. If needed, the next permitted route is the official FinnGen form or an author inquiry using the existing unsent draft; confirm file identity and full metadata before considering analysis. Do not substitute a DF13 disease endpoint, thresholded supplementary hits, or R12/R13 senility as the custom HFRS GWAS.

## Sources

- Mak et al. (2025), [Nature Aging article](https://www.nature.com/articles/s43587-025-00925-y), especially Methods and Data Availability.
- FinnGen, [Access results](https://www.finngen.fi/en/access_results), DF13 release and download instructions.
- FinnGen, [DF13 browser](https://r13.finngen.fi/).
- Internal bounded check: [`hfrs_r12_public_manifest_recheck_2026-09-24.md`](hfrs_r12_public_manifest_recheck_2026-09-24.md).
