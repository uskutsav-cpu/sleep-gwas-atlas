# Result-free MVP insomnia source review

Review date: 2026-10-08, America/Chicago. This is a source-admissibility review, not an analysis result. No newly estimated insomnia–disease correlation has been inspected. All association-body probes retained the header only.

## Source decision

**CONDITIONALLY FEASIBLE as independent related-phenotype validation; NOT an exact replacement of the frozen UKB insomnia GWAS.** MVP GCST90475826 measures clinically recorded insomnia (PheCode 327.4); the original UKB source measures frequent self-reported trouble falling asleep or waking at night. Both belong to insomnia, but matching the trait name does not establish equivalence of their genetic estimands. Before any pair outcomes, the new protocol must explicitly identify clinical-to-symptom phenotype transport and restrict its conclusions accordingly. A positive test cannot establish exact replication of UKB field-1200 liability or explain a clinical mechanism.

## Exact identity and composition

The current [Catalog v2 accession](https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90475826) explicitly identifies MVP, EUR only, 78,566 cases, 329,572 controls, N=408,138 and downloadable full summary statistics. The [primary dbGaP phenotype manifest](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/document.cgi?phd=8759&study_id=phs002453.v1.p1) independently maps these counts to `Phe_327_4`, type binary, category PheCodes/Neurological, and public/non-sensitive status. Its complete four-population row has N=575,155 (116,662 cases /458,493 controls): AFR 24,423/83,343; AMR 12,531/40,549; EAS 1,142/5,029; EUR 78,566/329,572. The selected file is the EUR stratum, not that four-population meta-analysis or the narrower organic/persistent-insomnia endpoint.

The [authors' primary methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10327290/) classify PheCode cases using at least two mapped ICD-9-CM or ICD-10-CM instances and controls using zero instances. Their published [2024 resource](https://pmc.ncbi.nlm.nih.gov/articles/PMC12857194/) also describes collated clinical diagnosis codes and public full results. The resource covers MVP participants linked to the VA EHR; it does not contain UKB or FinnGen in this selected GWAS. Exact ICD-to-PheCode mapping version and case-exclusion implementation should be recorded from the primary supplement or public phenotype definitions if needed for stricter phenotype harmonization; neither age/sex composition nor absence of comorbid conditions can be inferred from the accession totals.

Official body: [GCST90475826.tsv.gz](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475826/GCST90475826.tsv.gz).

| Identity field | Verified metadata |
|---|---|
| Compressed bytes | 575,504,178 (548.84 MiB) |
| Upstream MD5 | `2e9c12624b653aac527444fc426e037b` |
| HTTP ETag | `224d7f32-630b35ab64241` |
| Object last modified | 2025-03-19 15:05:39 GMT |
| Metadata modified | 2025-04-28 |
| Genome assembly | GRCh38; 1-based coordinates |
| Upstream format | GWAS-SSF v1.0, not harmonised, not sorted |
| EUR effective N from manifest counts | 253,768.615 (4×cases×controls/N) |

Build and format come from the exact [source YAML](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475826/GCST90475826.tsv.gz-meta.yaml); MD5 comes from that YAML and [official checksum file](https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475826/md5sum.txt). Publication date and file modification date are distinct. Whole-body MD5/SHA-256 are **not yet independently verified** in this review.

## Schema and direction

A verified-TLS, 65,536-byte range probe returned the following header, retained in `provenance_v1_header_evidence.json`:

```text
chromosome base_pair_location effect_allele other_allele odds_ratio standard_error effect_allele_frequency p_value rsid ci_upper ci_lower alt n case_af num_cases control_af num_controls r2 q_pval i2 direction
```

Use `effect_allele` as the beta allele and transform positive `odds_ratio` to log(OR). Higher values indicate increased clinically recorded insomnia liability. Do not infer the effect allele from `alt` alone. The historical MVP outcome materializer derives log-OR SE from log confidence bounds, `(log(ci_upper)-log(ci_lower))/(2×1.959963984540054)`, and validates beta/SE against P. The same checks are necessary here because the header alone does not establish the scale of `standard_error`. Record row-specific case/control N if it varies; the constant manifest effective N cannot silently replace substantial per-variant variation. Liftover or a directly verified GRCh38-compatible mapping is required for the frozen GRCh37/HapMap3/LDSC path. This schema evidence does not certify all rows or numerical QC.

## Rights and access

This accession specifies [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/), rather than an accession-level CC0 statement. Those terms impose no restrictions beyond original-owner terms and expect attribution. The published resource states public browsing/download of these non-sensitive results through dbGaP `phs002453`. Public full-summary scientific reuse is supported; protected individual-level data access is unnecessary for this file. This review does not establish an unrestricted redistribution sublicense. Keep raw data on the SSD, exclude them from Git, and attribute MVP/the primary publication and Catalog in later human-authored work.

## Disjoint-source design and remaining gates

Discovery is Jansen UKB-only insomnia × Pan-UKB EUR outcome. The proposed validation is MVP EUR clinical insomnia × the existing exact FinnGen R13 outcome. Thus neither validation source contains the UKB discovery cohort, based on their reported study composition. This is a stronger design than reuse of the UKB insomnia GWAS. Individual-level cross-enrollment is not directly audited: the justified description is **cohort-distinct by study design; exact participant intersection unknown**, not mathematically proven zero overlap. A nonsignificant LDSC cross-trait intercept cannot prove independence.

The complete 217-row result-free admission table is `provenance_v1_mvp_insomnia_217_admission.tsv`. There are 40 original insomnia candidates; 10 have existing FinnGen R13 external sources, two have MVP sources, and 28 have no eligible external source. The proposed pilot preserves every row and uses only the 10 FinnGen-source rows. Keep the frozen family denominator of 217 and alpha=0.05/217=0.0002304147465437788. Do not relabel unavailable rows as negative tests or reduce the denominator to 10.

The 10 metadata candidates are reflux disease, diaphragmatic hernia, noninfective gastroenteritis/colitis, unspecified chronic bronchitis, smoking dependency, COPD, gonarthrosis, arthrosis, intervertebral disc disorders and cholelithiasis. Six have historically exact outcome matches; four have documented definition differences. This preserves the original source selection. **Eligibility is not established by a previous positive outcome**: every selected row must undergo the same current source/hash/harmonization/h2/intercept/SNP-overlap gates, including sources that previously failed QC.

Before outcomes, freeze: related clinical insomnia estimand; exact source URL/MD5/ETag/size; verified full-body SHA-256 and harmonized-output hash; effect and N policy; source-specific filters; matched LD ancestry/build; phenotype differences; both sleep and outcome h2 Z≥4/intercept≤1.2 gates inherited from the contract; a declared power criterion and result-independent stop; direction and 0.05/217 success rule; all noneligible classifications. Do not substitute another source after results. The statistical reviewer must decide whether any effect-difference calculation is justified under independently sourced cohorts and whether uncertainty from unresolved individual overlap remains material. Only after those gates can this be called newly executed, cohort-distinct two-source validation of a related insomnia phenotype. There are currently **zero new validation results** from this review.

## Bounded acquisition plan

A single exact compressed body is approximately 549 MiB. Acquire once onto a new SSD output/data path, preserving original files; budget ≤600 MiB compressed plus ≤200 MiB streaming/HapMap3 outputs, and allow no whole-genome decompressed retention. Stream decompression and filtering while consuming the complete compressed body; retain compressed MD5/SHA-256, byte count, source headers and processing receipts. Verify available SSD capacity, native library/runtime and network suitability before execution. This review acquired no complete GWAS body.
