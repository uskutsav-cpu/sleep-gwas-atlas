# Physical-component GWAS supplement review

**Review date:** 2026-09-23. **Source:** Chen et al., *GeroScience*, DOI `10.1007/s11357-025-01734-2`, publisher-linked Supplementary file 1 (DOCX, 4.57 MB), downloaded from the visible Springer Nature supplementary-material link. The downloaded temporary source was identified as a Microsoft Word DOCX, 4,801,766 bytes; SHA-256 `3d556268324387129565220d2e59094fd4651425827963cf58085baa89f41c9c`. It was not copied into the repository; the source URL is the publisher asset linked from the article's “Supplementary Information” section: <https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs11357-025-01734-2/MediaObjects/11357_2025_1734_MOESM1_ESM.docx>.

## Evidence extracted

Table S3 explicitly gives UK Biobank field IDs and binary/missingness definitions for the same five physical-frailty component labels listed in the Zenodo record 14011550 README: weight loss (2306), exhaustion (2080), low physical activity (6164 and 1011), slow walking speed (924), and low grip strength (31, 21001, 46 and 47). The per-component operational coding is summarized in `physical_component_supplement_crosswalk.tsv`; it keeps the source-specific low-grip BMI/sex cutoff wording in that table rather than replacing it with a simplified threshold.

The article abstract describes its individual-level analysis of 49,530 UK Biobank participants with cardiovascular disease. Its publisher Data Availability statement says physical-frailty GWAS summary statistics were obtained from a previous study (reference 24 in the article, Ye et al. 2023). That statement does not identify the five separate Zenodo files as those same statistics. The Zenodo release's author, component names, and component definitions are concordant with the article supplement, which improves phenotype-definition evidence, but does not resolve the file-level summary-statistic lineage or demonstrate the GWAS sample/model behind the much larger per-file case/control totals.

## Decision

This confirms what the component phenotypes mean and which UK Biobank fields/response rules they use. It does **not** establish the five-file release's genome build, ancestry and variant/sample QC, signed effect interpretation/model, covariates, or exact generating analysis. Therefore it does not pass the frozen source-eligibility gate and does not authorize coordinate harmonization or LDSC. The Neale field-level GWAS remains a distinct proxy and is not substituted. No thresholds were changed and no result was calculated.
