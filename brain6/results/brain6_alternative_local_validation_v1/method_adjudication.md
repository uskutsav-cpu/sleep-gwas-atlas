# Brain6 secondary local genetic validation: method choice before outcomes

This branch is separate from the frozen seven-trait LAVA v3 endpoint. Its output
cannot turn the LAVA family into a pass or promote any protected candidate under
the LAVA rule. The choice below was made without reading any SUPERGNOVA local
covariance result or testing a PLACO candidate region.

| Method | Original methods and implementation | Brain6 suitability decision |
|---|---|---|
| **SUPERGNOVA** | [Zhang et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8422619/) derive observed-scale local covariance from paired signed GWAS Z scores and LD. Their case-control derivation preserves the significance and sign of covariance on conversion to liability scale. They estimate the sample-overlap/phenotypic covariance intercept from genome-wide data rather than requiring exact shared-person counts. Their [implementation](https://github.com/qlu-lab/SUPERGNOVA/tree/319e84e114a4f954005a4592756c56cfee083667) accepts one declared study sample size per GWAS and PLINK EUR reference genotypes, calculates a two-sided local covariance p value, and leaves local correlation undefined when either local h² is negative. | **Selected** for five prespecified Brain6 sleep–disorder pairs, conditional on exact input hashes, source-specific eligibility, successful common-SNP QC and complete family execution. Covariance, not correlation, is the primary secondary statistic. This method needs a scalar analyzed-study N and calibrated signed Z; it does not solve a genuinely unidentified GWAS model simply by accepting a number. |
| **ρ-HESS** | The [original HESS instructions](https://huwenboshi.github.io/hess/local_rhog/) specify a shared-sample count and phenotypic correlation to correct overlapping cohorts. The same page notes that the phenotypic correlation inferred from a cross-trait LDSC intercept also depends on that shared count. | **Not selected**: the relevant UK Biobank/consortium participant-level overlaps are not established. Supplying zero overlap or invented shared counts would be unjustified. |
| **HDL-L** | The [current HDL.L implementation](https://github.com/zhenin/HDL/blob/master/HDL/R/HDL.L.R) requires SNP-level `N` in each input and eigen-decomposed LD reference files. The [HDL-L methods](https://www.nature.com/articles/s41588-025-02123-3) describe likelihood-based local inference. | **Not selected for the five-pair family**: the exact Dashti long-sleep public source lacks per-SNP N, and the requisite source-matched HDL LD eigen-reference has not been established. Assigning the study total to every SNP would recreate the unsupported missingness assumption; the method's published UK Biobank LD eigen-reference may itself overlap UK Biobank discovery samples. |

## Input eligibility fixed before secondary outcomes

The paired association inputs are the archived GRCh37 EUR LDSC-format munged
files named in `munged_input_sha256.txt`, whose Z values derive from the exact
Brain6 GWAS sources. A new adapter removes missing/invalid Z, non-rs/duplicate
IDs and strand-ambiguous SNPs; it never changes a Z sign except by explicit
allele alignment to the EUR reference in the original SUPERGNOVA code. The
method uses the study-wide **total** N published for each source (not the
historical effective N written in some LDSC files). The N column in the munged
files is intentionally ignored by SUPERGNOVA when `--N1/--N2` are supplied.
This is a documented method-specific scalar-N approximation for meta-analyses
with variant-varying sample contribution, not a claim that every variant had
that exact N or an approved LAVA `N` reconstruction.

| Trait | Total N | Source/model eligibility |
|---|---:|---|
| Insomnia | 386,533 | Jansen 2019 exact UKB-only file has native per-SNP N; source README defines N and effect-allele OR, and signed Z is available. Eligible secondary. |
| ADHD | 225,534 | Demontis 2023 EUR meta-analysis supplies per-variant case/control participation, but a scalar study total approximates variant-varying N. Eligible secondary with N-heterogeneity limitation. |
| MDD | 500,199 | Howard 2019 public UKB+PGC meta-analysis has no per-SNP N; published public-release aggregate count is a scalar approximation. Eligible secondary with N-heterogeneity and mixed-phenotype limitation. |
| Schizophrenia | 130,644 | PGC3 EUR meta-analysis has per-variant effective case/control counts, but observed-scale association Z and public-release aggregate counts support secondary scalar-N analysis; eligibility limited by varying cohort contribution. |
| Bipolar disorder | 413,466 | PGC3 meta-analysis as above; eligibility limited by varying cohort contribution. |
| Parkinson disease | 482,730 | Nalls 2019 public EUR meta-analysis includes direct and UKB proxy cases, described as such throughout; eligibility limited by source phenotype mixture and varying cohort participation. |
| Long sleep | 339,926 | Dashti 2019 ≥9 h versus 7–8 h UKB contrast. Primary paper and coauthored follow-up identify linear BOLT-LMM; a pre-outcome genome-wide SE-scale audit supports the linear scale, but the release README does not explicitly bind columns to a model and lacks per-SNP analyzed N. SUPERGNOVA uses signed Z and scalar total N without LAVA's binary logistic reconstruction. **Eligible only as secondary, source-limited evidence**; no LAVA repair or exact per-SNP N claim follows. |

The [pre-outcome Dashti model audit](../lava_confirmatory_source_triage_v1/longsleep_original_model_20260927.md)
supplies the underlying 70,150-variant SE-scale evidence and exact-source
links. Published source definitions and counts are in
`brain6/manifests/gwas_external_availability.tsv`. Every pair uses its original
phenotype labels, including Parkinson/proxy and Howard broad depression. The
five pair IDs are `insomnia__adhd`, `insomnia__mdd`,
`longsleep__scz`, `longsleep__bipolar`, and `longsleep__parkinson`.

## Precommitted family and interpretation

Use all non-MHC Berisa–Pickrell EUR GRCh37 LDetect blocks and all five pairs.
Exclude reference MAF below 5%, non-rs or duplicated IDs, non-biallelic SNPs,
and A/T or C/G ambiguous alleles before analysis. Use a GRCh37 1000 Genomes
Phase 3 EUR (503-person) PLINK LD panel with HapMap II GRCh37 genetic-map cM
interpolation; the downloaded BIM has zero cM and is invalid unmodified for
SUPERGNOVA's genetic windows. The panel is split by chromosome without
changing genotype records, since the implementation requires sorted cM within
each LD-score pass and genetic cM resets at chromosome boundaries. A combined
22-chromosome file must not be passed to its `getBlockLefts` routine. Report
skipped blocks below the method's 120
shared SNP minimum as `NOT_TESTED_LOW_SNP` and all other missing blocks as
errors, never as null results.

Use two-sided Bonferroni `0.05 / (5 × number of frozen non-MHC blocks)` for
the covariance p value across the entire five-pair family, regardless of
whether a block returns an estimate. No candidate-only correction is permitted.
The local correlation estimate is descriptive only if both local h² are
nonnegative and the software emits a finite correlation. A protected candidate
gets `SUPPORTED_FWER` only if one or more predeclared LD blocks overlapping
its frozen coordinate interval has a finite covariance p value below this
family threshold and the full pair run passes input/runtime completeness QC.
The sign is the sign of local covariance and is compared descriptively with
the lead-variant cross-trait direction; discordance is never concealed. No
significance choice, block boundary, source, N arm, or method changes after
seeing local results. A nonsignificant candidate is `NOT_SUPPORTED`, not a
disproof of shared biology; an unestimated block is `NOT_TESTED_INPUT_QC`.

The final paper must present this as **secondary local covariance evidence**
alongside the failed prespecified LAVA endpoint. It is not independent cohort
replication: these are the same GWAS associations analyzed by another method.
