# FinnGen R9 replication candidates for preselected Brain6 pairs

**Audit date:** 2026-09-25. **Status:** metadata feasibility only; no R9 GWAS summary-statistic bytes were acquired or analyzed.

## Scope and candidate selection

This follow-up targets only the four already selected, original-family-significant Brain6 pairs without an admitted replication estimate: insomnia–MDD, long sleep–bipolar, long sleep–schizophrenia, and long sleep–Parkinson's disease. FinnGen Data Freeze 9 (R9) is a separate Finnish population cohort and predates the currently archived R13 candidates. It is a plausible external-cohort directional replication source, subject to differences in ancestry, diagnostic ascertainment, phenotype definition, and exact participant overlap. It does not make the insomnia exposure independent because the Brain6 insomnia GWAS remains the same UK Biobank exposure.

The official R9 manifest was retrieved from `https://storage.googleapis.com/finngen-public-data-r9/summary_stats/R9_manifest.tsv` (661,463 bytes; SHA-256 `13f3866a224e35471abc0290cba5ed151aad7a4e270ae70b3eacf3a92cf124a2`). Four relevant endpoint objects are listed there. Anonymous HTTP HEAD requests returned 200 and exposed object generation, byte count, ETag, and modification time; these are metadata observations, not data acquisition receipts. Per-endpoint details are in `finngen_r9_replication_candidates.tsv` and its provenance JSON.

| Pair | R9 endpoint | Cases / controls | Approximate effective N* | Interpretation |
|---|---|---:|---:|---|
| Insomnia–MDD | `F5_DEPRESSIO` | 43,280 / 329,192 | 153,004 | Register-defined depression (ICD-10 F32/F33); related to MDD, but narrower than the UKB help-seeking plus PGC discovery phenotype. |
| Long sleep–bipolar | `F5_BIPO` | 7,006 / 329,192 | 27,440 | Bipolar affective disorder from Finnish hospital-discharge/death registers; diagnosis and ascertainment differ from PGC3. |
| Long sleep–SCZ | `F5_SCHZPHR` | 6,515 / 364,160 | 25,602 | Schizophrenia from Finnish hospital-discharge/death registers; narrower than some PGC case definitions. |
| Long sleep–Parkinson's | `G6_PARKINSON` | 4,235 / 373,042 | 16,750 | Register-defined Parkinson's disease; case count may limit precision for the discovery effect. |

\*Computed from manifest case/control counts as `4 × cases × controls / (cases + controls)`; this is a summary effective-N approximation, not the per-variant N from GWAS files.

The R9 source describes its summary-statistic coordinate/effect conventions in the official format documentation; a candidate file still needs post-access header, build, effect-allele, beta/SE/P/N, variant-count, checksum, and transformation-QC verification before analysis. Finnish ancestry is European-related but founder structure limits direct portability claims. MDD, bipolar, and schizophrenia are clinically related to the source outcomes, not exact phenotype matches. No sample-overlap-zero assertion is made; exact individual overlap with source GWAS cohorts was not verified.

## Additional Parkinson's source screen: FinnGen R4 fallback

The 2024 Parkinson's multi-ancestry meta-analysis reports its European Nalls et al. 2019 stratum separately from FinnGen R4, and reports that the FinnGen data did not include UK Biobank participants. Its FinnGen R4 `G6_PARKINSON_EXMORE` row contains 1,587 cases, zero proxy cases, and 94,096 controls (approximate effective N 6,243). This is a source-specific candidate for a comparable Finnish-register PD outcome against the locked Nalls 2019 discovery, but only the standalone R4 stratum would be eligible; the combined 2024 meta-analysis is not independent because it contains the Nalls European results. The R4 download instructions also require an online form and emailed instructions. No R4 summary-statistic bytes, release manifest entry, file schema/build, or source checksum were acquired, so it is not admitted.

The larger R9 `G6_PARKINSON` candidate (4,235 cases; approximate effective N 16,750) remains the first-choice metadata candidate if authorized FinnGen access is obtained. R4 is retained as a lower-power fallback, not as an additional estimate or a reason to change the long-sleep–PD pair selection. Machine-readable source screening is in `finngen_pd_r4_fallback_candidate_20260925.tsv` and its provenance JSON.

## Access and admission decision

The official FinnGen R9 download instructions require completion of the online access form followed by emailed download instructions. That workflow has not been completed. No GWAS file was downloaded, no access gate was bypassed, and no replication estimate or classification was changed. Consequently these are **candidate sources only**, not admitted results.

After authorized access, retain the original bytes in the research archive and bind each to the release/object generation, byte count, cryptographic hash, acquisition event, file header, and an explicit raw-to-analysis transformation receipt. Harmonize only then. Preselect the four endpoint mappings above before calculating rg; evaluate direction and uncertainty against the locked Brain6 pair estimates and retain Finnish-founder and phenotype-definition caveats. Do not call the reused insomnia exposure a fully independent two-trait replication.

## Sources

- FinnGen R9 manifest: https://storage.googleapis.com/finngen-public-data-r9/summary_stats/R9_manifest.tsv
- Official R9 download/access instructions: https://finngen.gitbook.io/documentation/r9/data-download
- FinnGen R9 release background: https://finngen.fi/en/public-release-finngen-data-freeze-9-results-and-summary-statistics
- R9 endpoint definitions: https://r9.risteys.finngen.fi/endpoints/F5_DEPRESSIO ; https://r9.risteys.finngen.fi/endpoints/F5_BIPO ; https://r9.risteys.finngen.fi/endpoints/F5_SCHZPHR ; https://r9.risteys.finngen.fi/endpoints/G6_PARKINSON
- Kim et al. 2024 PD multi-ancestry meta-analysis, including source-specific FinnGen R4 and Nalls European strata: https://doi.org/10.1038/s41588-023-01584-8 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10786718/
- Official FinnGen R4 download instructions: https://finngen.gitbook.io/documentation/r4/data-download
