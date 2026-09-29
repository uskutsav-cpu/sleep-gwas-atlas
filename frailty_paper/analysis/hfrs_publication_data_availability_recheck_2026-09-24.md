# HFRS publication data-availability recheck

Audit time: 2026-09-24 02:37 UTC.

## Evidence from the primary study

Mak et al. report a custom continuous Hospital Frailty Risk Score (HFRS) GWAS in FinnGen Release 12 (N=500,737), UK Biobank replication (N=407,463), and 53 independent lead variants for the main HFRS definition. Their Methods define the score from 109 weighted ICD-10 codes and analyze the score continuously; a separate sensitivity analysis removes dementia weights. These are the study's reported design facts, not validation of an available file.

The article's Data availability statement says individual-level data cannot be made publicly available, and points to the FinnGen results route for summary statistics after a one-year embargo. The article does not name a study-specific HFRS endpoint/accession, direct HFRS summary-statistics URL, or file checksum. It also does not identify a full-statistics link alongside its supplementary-material references.

## Current official FinnGen route

FinnGen's current access page says Data Freeze 13 (DF13) was publicly released on 2026-06-02 and provides 2,755 standard disease endpoints for 500,186 participants. The same page says summary-statistic downloads are initiated by completing an online form, after which download instructions are emailed. Its interactive browser is public. Neither the release announcement nor the general access page identifies the paper's custom R12 HFRS or HFRS-without-dementia file as a standard endpoint.

## R12 sample-count reconciliation

The HFRS article's reported FinnGen N=500,737 differs from the official R12 core-analysis count of 500,348 by 389. FinnGen's release documentation explains the two counts as different sample-processing stages: 500,737 nonduplicate Finnish-ancestry inliers had covariates projected, then 355 samples lacking minimum phenotype data and 34 failing sex checks were excluded, leaving 500,348 for core analyses. Thus, the count difference alone does not establish an error in the HFRS paper; it leaves the exact HFRS analytic inclusion set and its relation to the R12 core set unresolved. Require the full-file metadata or author/FinnGen confirmation before treating the reported N, controls, or sample-overlap assumptions as reconciled.

Sources for this reconciliation: FinnGen's [R12 release announcement](https://www.finngen.fi/en/results-based-full-finngen-cohort-500000-participants-released) reports 500,348 participants; the [FinnGen release documentation](https://docs.finngen.fi/finngen-data-specifics/finngen-data-freezes-and-releases) distinguishes total samples from core-GWAS samples; the [FinnGen QC documentation](https://finngen.gitbook.io/documentation/methods/phewas/quality-checks) describes the 500,737 and 500,348 processing stages. The [HFRS primary article](https://www.nature.com/articles/s43587-025-00925-y) reports N=500,737.

## Access decision

- Exact FinnGen R12 HFRS and HFRS-without-dementia summary files remain **NOT ACQUIRED / NOT VERIFIED**.
- The paper's N=500,737 matches FinnGen's pre-core projected sample count, while the official R12 core-GWAS count is 500,348; the HFRS-specific analytic sample's QC stage has not been verified from a full file or study-specific metadata.
- A newer standard FinnGen release or a similarly named diagnosis endpoint is not evidence of equivalence to the published custom score.
- No web form was submitted, no email or author inquiry was sent, and no access restriction was bypassed.
- The existing unsent author/FinnGen inquiry draft remains the authorized next route if exact-file access is needed; before use, confirm release, phenotype construction, N, build, effect scale, alleles and checksum.

## Primary sources

- Mak et al. (2025), *Nature Aging*, [study article](https://www.nature.com/articles/s43587-025-00925-y), PMID 40764432. The article reports the R12 custom-score analysis, defines the weighted-code HFRS, and directs readers to FinnGen for post-embargo summary statistics.
- FinnGen, [Access results](https://www.finngen.fi/en/access_results). The current page records the DF13 release, participant/endpoint counts, interactive browser, and form/email download steps.

## Official R12 manifest follow-up (2026-09-24)

A later read-only streaming audit of the official public R12 summary-statistics manifest found 2,469 listed rows and no `HFRS`/frailty label matches. The check searched all manifest values, not the contents of every result file. FinnGen’s standard R12 core GWAS documentation also describes a Regenie analysis of 2,502 core endpoints/500,348 samples, distinct from the custom HFRS SAIGE analysis reported in the paper (N=500,737). This narrows the standard-manifest route but does not establish that no custom file was separately delivered or stored under an opaque name. Exact full custom HFRS files remain unverified; see `hfrs_r12_public_manifest_recheck_2026-09-24.md` for the command, checksum, result and sources.

## DF13 route recheck (2026-09-24 20:03 UTC)

The one-year publication embargo has elapsed, and FinnGen's current access page identifies DF13 as publicly released on 2026-06-02 (2,755 standard endpoints; 500,186 participants). However, its bulk-download instructions still require an online access form followed by emailed instructions. The article used a custom continuous 109-code HFRS score in R12; release of newer standard disease endpoints does not verify availability or equivalence of that custom phenotype. The exact full HFRS files remain not acquired/not verified. Details: `hfrs_df13_release_route_recheck_2026-09-24.md`.
