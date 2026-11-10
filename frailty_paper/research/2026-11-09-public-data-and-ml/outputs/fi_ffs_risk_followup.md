# FI–FFS genetic correlation and published frailty-locus follow-up

Research update, 27 October 2026. This is a source audit and reconstruction from saved project outputs and a published supplement. It is not a new GWAS, a new joint-locus test, or an individual risk prediction.

## Direct FI–FFS correlation

Two saved auxiliary LDSC runs directly compared the Atkins FI GWAS with the Ye FFS GWAS. They used SNP sets restricted to the support of two insomnia-source analyses; the labels identify which set was used, not the phenotype being correlated in this table.

| Matched SNP set | FI–FFS genetic correlation | Block-jackknife SE | Approximate 95% CI | SNPs |
|---|---:|---:|---:|---:|
| FinnGen-insomnia matched set | 0.7570 | 0.0212 | 0.7154–0.7986 | 1,138,262 |
| MVP-insomnia matched set | 0.7548 | 0.0200 | 0.7155–0.7940 | 1,151,928 |

I checked each saved exact file against its recorded SHA-256, recomputed each point estimate as genetic covariance divided by the square root of the two SNP-heritability estimates, and independently recalculated its SE from the 200 delete-block estimates. The results matched the saved estimates to numerical precision; the 95% intervals are normal approximations. Full traceability appears in `direct_fi_ffs_genetic_correlation.csv` alongside this note.

These two estimates use the same FI and FFS GWAS, so they are alternative-SNP-support sensitivities, not independent replications. UK Biobank contributes to both sources; exact person overlap remains unknown. The results support substantial genome-wide genetic overlap between these particular measures; they are not a correlation of measured clinical scores and do not imply that 75% of variants, patients, or disease risk are shared. They do not explain why insomnia correlates differently with the two measures. The prior note's statement that the insomnia contrast was not the direct FI–FFS correlation remains true, but this additional direct estimate was subsequently located.

## Published FFS locus inventory

I examined the [Ye et al. 2023 supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1007%2Fs11357-023-00771-z/MediaObjects/11357_2023_771_MOESM2_ESM.xlsx), sheets `ST1. COJO independent loci` and `ST3. GWAS Catalog`. The downloaded file's SHA-256 is `cc95a22e0c8731473d29c72951ee3b887bc9344f87bbb01f73dc4ee8709282e9`. The accompanying `published_ffs_37_locus_inventory.csv` transcribes all 37 COJO loci from ST1, including the source's nearest-gene annotation, marginal and conditional statistics, and counts of previously catalogued sleep-related entries in ST3. These are **FFS-associated loci**, not 37 established insomnia–frailty shared-risk loci. The source's nearest gene is an annotation, not a demonstrated effector gene.

All 37 source rows have conditional p-values below 5×10⁻⁸, but one—the MAP3K3 annotated indel `17:61736108_TTTTG_T`—has a marginal p of 6.87×10⁻⁸ and conditional p of 3.85×10⁻⁸. Both columns are retained rather than describing every marginal test as genome-wide significant. `refA` and betas are transcribed as reported; precise build/allele semantics should be authenticated before harmonization or scoring. The estimated conditional effects are not prevalidated weights for an insomnia–frailty risk score. The population is the FFS study population of European-ancestry UK Biobank participants; overlap with sleep GWAS sources must be evaluated.

Thirteen previously catalogued sleep or circadian association entries occur in eight of the source's 37 locus groups. The companion `published_ffs_loci_sleep_catalog_links.csv` retains each sleep trait and variant, the FFS lead, reported gene mappings and position difference. This reflects the source authors' catalog grouping; proximity alone neither establishes LD nor demonstrates that the traits share a causal variant. It does not prove same-direction effects, independent replication, or a causal gene. The sleep entries include chronotype and sleep duration as well as insomnia and should not be pooled as one trait.

### Insomnia-specific example

For locus 14, the FFS COJO lead is **rs7703746** (reported chromosome 5 position 103,907,337), with marginal FFS beta 0.01315 and p=6.98×10⁻¹³. ST1 reports nearest gene `RP11-6N13.1`. The same source locus grouping includes **rs1592757**, reported at 103,889,998 in the insomnia-symptom GWAS of [Lane et al.](https://pubmed.ncbi.nlm.nih.gov/30804566/), p=2×10⁻¹⁰; ST3 reports `NUDT12/RAB9P1` and maps it to `AC099520.1`. The two reported positions differ by 17,339 bp. The variants are distinct and the original studies may share UKB participants. A subsequent 1000 Genomes European-reference pairwise lookup returned r²=0.6433, but the exact study allele contrasts remain unresolved for these multiallelic rsIDs (see the LD section below). Genome-build and allele coding must be verified on the actual input releases before any direct positional or signed-effect comparison. This is a **candidate region**, not a validated causal variant, confirmed gene, or shared-risk score.

A different prior study reported two additional candidate variants, rs34290943/DALRD3 and rs10865954/MIR191, from an insomnia–FI analysis at 3p21.31. These are separate from the chromosome 5 FFS candidate and are not established FFS findings. [Song et al. 2024](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2024.1358996/full).

## What the inventories can support

The 37-row FFS list and 13 sleep-catalog links give a concrete, auditable starting set for future risk investigation. They cannot calculate an individual's risk or identify 15–30 verified shared variants. The project’s frozen FI×sleep LAVA analysis completed all slots but failed all 12 prespecified trait QC gates. A new region-based screen would need a separately justified, prespecified design, full dense source data, source/allele/build and overlap checks, suitable LD and independent confirmation; failed results must remain failed. Identifying a shared signal is followed by credible-set/colocalization and gene mapping, not guaranteed by nearness.

For a real polygenic score, the endpoint, weights, effect alleles, and development population must be defined and its performance evaluated in independent genotyped participants. Existing [PGS005229](https://www.pgscatalog.org/score/PGS005229/) is an FFS-derived comparator score, not an insomnia–frailty combined score. No usable shared score was calculated here.

## Published FI and FFS polygenic-score comparison

[Foote et al.](https://www.nature.com/articles/s41588-025-02269-0) evaluated an FI-derived score ([PGS005228](https://www.pgscatalog.org/score/PGS005228/)) and an FFS-derived score ([PGS005229](https://www.pgscatalog.org/score/PGS005229/)) against observed FI in three older-adult cohorts. Both used SBayesR and European-ancestry discovery data; the FI score contains 2,379,429 variants and the FFS score 2,761,868. The catalog reports the following R² values; the paper describes the analysis as **incremental variance explained above age, sex and ancestry principal components**:

| Observed outcome and evaluation cohort | N | FI-derived PGS R² | FFS-derived PGS R² |
|---|---:|---:|---:|
| FI in ELSA | 7,181 | 0.0200 | 0.0138 |
| FI in PISA | 3,265 | 0.0210 | 0.0300 |
| FI in LBC1936 | 1,005 | 0.0050 | 0.0190 |

**Source discrepancy:** A later audit found that several PGS Catalog values in this table conflict with the paper’s Supplementary Table 44, including differences beyond rounding. The cause is unresolved; do not use these numbers as verified Table 44 values or rank scores on this basis. See `fi_ffs_pgs_published_validation_audit.md` for the side-by-side source values and physical-frailty results. These are cohort-specific reported model statistics, not lifetime risk, probabilities of becoming frail, or a comparison of insomnia–frailty shared scores. They cannot establish which score is generally superior without a paired statistical comparison, nor do they guarantee calibration in a different population. Both discovery and evaluation ancestry were reported as European, so transferability is unproven. The evaluation FI differs in its number of component deficits across the three cohorts. The catalog's reported OR and beta are supplied alongside linear-model results in the paper, so the R² values here should not be used to infer absolute event risks or the scale of an individual risk estimate.

## Published insomnia–FI cross-trait candidates: 124 variants

The full supplemental file for [Song et al.](https://doi.org/10.3389/fnagi.2024.1358996) was recovered through [Europe PMC's supplementary-files service](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10903740/supplementaryFiles). Its `Data_Sheet_1.docx` has SHA-256 `6a1483cdae3997fdc46420e182f49cf0e7b1591b587806bbb7a4dd46fd27ee28`. Table S4 contains 124 distinct listed variant IDs. I transcribed all rows to `published_insomnia_fi_124_cross_trait_candidates.csv` and made `published_insomnia_fi_top30_candidates.csv` by sorting those rows by the reported CPASSOC p-value. All 124 rows have reported MTAG p-values for each trait and CPASSOC p-values below 5×10⁻⁸.

**Interpretation is crucial:** These are the paper's cross-trait statistical candidates. A CPASSOC association can be driven by one trait, while MTAG results share information across traits; neither alone proves that each variant independently influences both traits in the same direction. The table contains no effect estimates, confidence intervals, verified variant-to-gene links, or risk-score weights. Its genome build and reported effect-allele interpretation have not yet been independently authenticated. The original FI GWAS and insomnia source also contain UK Biobank, so exact overlap must be assessed. Reproducing this published table is not an independent replication or a new project association.

The top-30 list is an **illustrative ranking**, not 30 independent risk loci or a predictive score. Fifteen of its 30 entries lie within chromosome 6 positions 25–34 million, an approximate MHC flag whose exact boundaries and genome build require verification. Nearby variants can be correlated, particularly in the MHC; LD independence was not established. The table includes all 124 candidates so that the ranking and its concentration are visible. It also preserves the two 3p21.31 variants named in the paper’s main text: rs10865954 ranks 29th and rs34290943 ranks 95th by the reported CPASSOC p-value. These ranks do not establish relative causality.

A scientifically valid sparse score would require signed and allele-aligned weights from appropriate discovery statistics, a frozen selection procedure with LD handling, and independent genotype-level evaluation. The paper's S4 p-values alone cannot provide those weights. No such score has been claimed here.

## Descriptive lookup in the project's three GWAS inputs

On 27 October, I looked up **all 124** published Table S4 rsIDs in the existing harmonized Jansen UKB-only insomnia input, the harmonized Atkins FI input and the article-linked Ye FFS file. This was a post-publication descriptive lookup, not a new preregistered discovery family, replication test, local-sharing analysis or test of causality. The project’s failed LAVA gate remains unchanged. The complete one-row-per-candidate results, including absent IDs, original per-input alleles, beta, SE and p, are in `published_candidate_marginal_lookup.csv`.

| Source input | Candidates found by rsID | Matching reported chromosome and position |
|---|---:|---:|
| Jansen insomnia, harmonized | 99/124 | 99/99 |
| Atkins FI, harmonized | 103/124 | 103/103 |
| Ye FFS, raw article-linked file | 110/124 | 110/110 |

All 99 insomnia matches also appeared in the FI and FFS files. Among these 99, **none** had both insomnia and FI *marginal-source* p-values below 5×10⁻⁸. This is a descriptive check of these particular sources, not a test that there is no shared risk: the published multi-trait statistics answer a different question, sample sizes and phenotype definitions differ, and lack of significance is not absence of effect. In particular, no multiple-testing or selection-bias adjustment has been performed for this lookup. Across the 99 non-palindromic, position- and directly allele-matched FI–insomnia rows, all aligned beta signs were the same; that observation is conditional on published candidate selection and shared or overlapping source cohorts, and is not evidence of independent replication. For FFS, 87 aligned signs matched insomnia and 12 differed, assuming the documented-format inference that Ye's A1 is the effect allele; individual signs may be noisy and these counts are not a test of phenotype differences.

The project’s harmonization excluded the approximate chromosome 6 MHC, strand-ambiguous variants and other QC failures. Of 25 insomnia IDs not found, 21 have source-reported chromosome 6 positions in the 25–34 Mb filter interval. The remainder can have multiple exclusion reasons; absence is **not a zero effect**. All 110 FFS matched positions agree with the published table, but this does not by itself resolve exact-build, strand and source-specific effect conventions.

For example, rs113851554 is ranked first by Song’s CPASSOC p-value (4.07×10⁻⁴⁶). In this lookup its insomnia marginal p was 2.347×10⁻²¹, FI marginal p 0.1987, and FFS marginal p 0.19. Those p-values come from non-identical analyses and do not prove the FI or FFS effects are absent. They illustrate why a strong cross-trait rank is not enough to establish a significant effect on each individual trait.

Input provenance: Jansen insomnia harmonized file SHA-256 `81c71d838b19f0297e09ec027d00382732fee1bfd112f44cfd2feca44a4f8bf5`; Atkins FI harmonized file `241ed294ccc367b22787a0c051f8d377c020e16089e82920a555286f3148e6de`; Ye raw file `912e290b2064999c7147cbe3597f87cbcb14b25d387a8b6db2f25862e7274e55` (matches the project's registered exact file). Files were scanned completely for rsID and checked for duplicate matches. The saved QC reports use GRCh37/hg19 for the harmonized data; the exact FFS build is supported by coordinate concordance and its A1/BETA convention remains an inference rather than an author-confirmed sidecar. Direct/swap matching was used without strand complementation, with A/T and C/G pairs excluded from directional comparisons. Matching rsIDs/positions cannot by themselves verify original alleles, causality or independence.

## An unvalidated 15-variant score prototype

To make the requested score concrete, `unvalidated_15_variant_score_prototype.csv` contains **15 existing GWAS-weight rows**, not a validated shared-risk score. The rule is transparent: take the 30 published Song rows with the smallest reported CPASSOC p-values and retain those present in all three supplied inputs after the project's existing QC. Exactly 15 remain. This availability-driven, post-result choice is not optimized or prespecified independently of discovery results; removal of MHC is inherited from the project filters, and the 15 have **not been LD-pruned**.

The table aligns each trait’s signed beta to the same counted allele. For a person with observed dosage `d_i` of that allele, a researcher could calculate three separate *raw research indices*, one each for FI, FFS and insomnia:

`index_trait = sum(d_i × beta_trait_i)` for i=1…15.

These are **three distinct uncalibrated weighted sums**, not a single shared insomnia–frailty score, and they are not probabilities or clinical risk estimates. FI weights come from the FI source; FFS weights come from Ye and rely on inferred A1 semantics; insomnia weights come from the Jansen source. Trait-specific effect scales differ, so the three columns cannot be added together. No genotype data were applied, missing-genotype policy was evaluated, ancestry transfer checked, score tested in an independent cohort, or calibration established. The weights include variants with weak FI or FFS marginal evidence and selection was not LD-aware. This prototype should not be used to rank people, make health decisions, or claim a predictive result.

A legitimate shared-risk model would need a defined outcome and time horizon, a prespecified LD-aware feature selection and weight-estimation method, robust input conventions, and independent individual-level validation with calibration and incremental performance against existing FI, FFS and insomnia PGS baselines. The prototype records the missing work rather than filling it with invented validation.

## Clinical insomnia sources and pairwise LD checks

The existing project includes two additional clinical insomnia datasets: FinnGen R13 (51,643 cases, 446,273 controls) and European-ancestry MVP (78,566 cases, 329,572 controls). I checked the full 124 published candidate IDs against the project's existing HapMap3-filtered files, whose SHA-256 values matched the saved run manifests. Both filtered files contained 31 candidate IDs; 30 could be compared directionally with the harmonized Jansen input by rsID, chromosome and non-palindromic allele matching. The filtered files do not contain coordinates, so the match is not an independent position check. Direct, swapped and unambiguous strand-complement matches were considered, with effect signs changed for swapped orientation. Full source values and non-matches are in `published_candidates_external_insomnia_lookup.csv`.

Among the **15 prototype variants**, only seven were present in both external HapMap3 files with a comparable Jansen allele. All seven had the same aligned sign in FinnGen; five matched and two differed in MVP. These counts are **post-selection descriptive summaries**, not formal tests: effects have uncertainty, missing HapMap3 markers are not null results, the data sets use different insomnia definitions, and exact participant overlap and relatedness are unresolved. The FinnGen and MVP inputs change the insomnia source but reuse the same FI and FFS GWAS elsewhere; they do not constitute independent frailty validation. MVP's filtered source lacks an SE column; no SE was invented. No multiple-testing correction or model comparison was performed for these lookups.

To assess LD without asserting score validity, I queried [Ensembl's GRCh37 pairwise LD endpoint](https://rest.ensembl.org/documentation/info/ld_pairwise_get) for 1000 Genomes phase 3 European-reference estimates. Three pairs within 500 kb among the 15 prototype markers returned r² values 0.0949, 0.1708 and 0.1369. These are only the three nearby pairs queried; all other pairs and population-specific LD remain untested. A different pair, the FFS lead **rs7703746** and insomnia-catalog **rs1592757**, returned r²=0.6433 and D′=0.9942 in that reference. This supports following up the chromosome 5 region but does not establish colocalization, a shared causal allele, or the same direction of effect. The two Song chromosome 3 candidates returned r²=0.0612. Exact responses, source URLs, pairs and hashes are recorded in `selected_variant_pair_ld_reference.csv`.

Ensembl GRCh37 variation records agree with the reported chromosome and positions for these four candidate variants, but list multiple alleles for rs7703746, rs1592757 and rs10865954. The pairwise LD response does not specify which alleles were used. Therefore the reported pairwise r² must not be assumed to describe the exact study-effect alleles until the reference's multi-allelic representation is verified. One reference estimate also cannot replace ancestry-matched locus-level LD, fine-mapping, score LD handling or replication.

### Variant identity audit of the score prototype

The Ensembl GRCh37 variation endpoint returned the expected chromosome, position and listed alleles for all 15 prototype IDs. However, **nine of the 15 rsIDs have more than two alleles listed in the reference**. The remaining six have two listed alleles. This does not prove the GWAS used an invalid biallelic contrast, but it requires the exact reference/alternate allele representation, decomposition and study effect coding to be checked before scores or LD are interpreted. The prototype now marks these variants explicitly; the source responses and hashes are documented in `prototype_variant_identity_reference.csv`. Coordinate and allele-list agreement is not verification of allele-specific frequency, LD, weights or score validity. No individual scoring was performed.

## Full-release follow-up of all 15 prototype markers

A subsequent full scan recovered **all 15 variants** in both the raw FinnGen R13 insomnia release and the raw MVP European insomnia release. This supersedes the earlier seven-marker *HapMap3-subset coverage* as a coverage statement; the earlier seven-marker results remain correct for those filtered files. The raw inputs have 21,327,010 and 19,703,815 source rows, respectively. Each marker had exactly one rsID-matching row in each file. The scanned compressed-file hashes matched the saved acquisition/QC manifests: FinnGen `353129dc8252461a2a97087ce6de39548d0b0bf1e488a6ef36da75e9f4690049`; MVP `46b62349ba0d98b22bf302be973bb54a80707cfb02ce953426f8224910761835`.

FinnGen uses GRCh38, and the MVP source's saved GWAS-SSF sidecar states GRCh38, one-based coordinates. All 15 project GRCh37 positions received a single Ensembl GRCh37→GRCh38 mapping on the forward strand; the resulting coordinates matched both raw releases and the Ensembl GRCh38 rsID mappings. Each source's reported allele pair directly matched or swapped the prototype allele pair. Results and provisional orientations are in `prototype_15_variants_full_external_insomnia_lookup.csv`; mapping URLs and reference-response hashes are in `prototype_variants_build_reference.csv`.

These checks establish coordinate and allele-string agreement, not full variant-level QC or proof that multi-allelic contrasts and sample independence are resolved. Nine markers have more than two reference alleles, requiring source-by-source contrast verification. FinnGen raw beta refers to ALT; MVP signed effects were calculated as log of the reported odds ratio, which is rounded in the source, and all 15 MVP SE values are `#NA`; none was invented. Source-specific INFO/imputation quality for these raw, non-HapMap3 markers was not requalified. No tests were adjusted for multiple comparisons, and the 15 were selected after the original study's results were known.

| Comparison with Jansen insomnia effect direction | Same | Opposite | Not interpretable |
|---|---:|---:|---:|
| FinnGen R13, all 15 | 14 | 1 | 0 |
| MVP, all 15 | 13 | 2 | 0 |

Direction agreement alone does not establish a statistically reliable effect; raw effects and p-values are retained. Exact participant overlap and relatedness are unknown. Both datasets evaluate **insomnia**, not frailty; even an independently established insomnia association would not validate a shared insomnia–frailty score.

One concrete example is rs113851554: the original Jansen insomnia marginal p=2.347×10⁻²¹, FinnGen insomnia p=1.162×10⁻¹⁰ and MVP insomnia p=4.216×10⁻¹⁵, with matching provisional effect direction. But its FI and FFS marginal p-values in the project's other inputs are 0.1987 and 0.19, respectively. This is a promising **insomnia-associated marker**, not evidence that the marker is independently associated with both frailty outcomes. Non-significant p-values do not establish zero effects, and no causal or shared-locus conclusion follows from the comparison.

Sources: [FinnGen genome-build documentation](https://docs.finngen.fi/finngen-data-specifics/red-library-data-individual-level-data/genotype-data/human-genome-build-we-use-and-how-do-we-report-genotypes), the saved GWAS Catalog metadata sidecar for GCST90475826, and the Ensembl mapping and variation records. This was an exploratory source lookup, not a new preregistered association study, a validation of the 15-marker score, or a replacement for the failed project locus gate.

## Fifteen-marker comparison across FI, FFS and three insomnia sources

`prototype_variant_five_trait_evidence.csv` brings the 15 exploratory markers into one table. Effects are aligned to the listed counted allele; FI and FFS remain in their own original phenotype scales, while insomnia effects reflect the individual study scales. The table adds approximate **unadjusted** Wald 95% intervals where source SEs were available: beta ± 1.95996 × SE. MVP's supplied SEs are `#NA`, so no MVP interval was calculated. Ye’s FFS effect-allele interpretation is a documented format inference, and nine reference IDs are multiallelic. These confidence intervals are conditional on the source estimates; they do not adjust for prior variant selection, cohort overlap, variant uncertainty or multiple testing.

Eight of the 15 selected markers have source-reported marginal p<0.05 in both FI and FFS, with the same aligned sign. The 0.05 threshold here is a **descriptive screen**, not a corrected familywise or discovery claim. The candidates were chosen from prior results, FI and FFS share UK Biobank participants, the variant identities require further work, and the FI phenotype itself includes sleep-related items. These observations do not prove eight shared causal variants or independent replication. A more stringent threshold is not a substitute for replication or fine-mapping, and a non-significant result is not proof of absence.

Two examples, selected for discussion after looking at the table, illustrate why the five-source evidence is nuanced:

* **rs323509**, counted allele A: FI beta 0.0169 (approximate 95% CI 0.0100–0.0238); FFS beta 0.0120 (0.0082–0.0158). The source p-values were 1.802×10⁻⁶ for FI, 8×10⁻¹⁰ for FFS, 6.826×10⁻⁹ for Jansen insomnia, 1.5771×10⁻⁴ for FinnGen insomnia, and 0.05089 for MVP insomnia. Ensembl lists multiple alleles for this rsID. The original FI and FFS measurement scales differ, so beta magnitudes are not directly comparable.
* **rs12536335**, counted allele G: FI beta −0.0158 (approximate 95% CI −0.0225 to −0.0091); FFS beta −0.00579 (−0.00944 to −0.00214). The source p-values were 2.633×10⁻⁶ for FI, 0.0017 for FFS, 7.651×10⁻⁷ for Jansen insomnia, 0.00625446 for FinnGen insomnia and 4.57×10⁻⁸ for MVP insomnia. The negative effects are relative to the specified counted allele, not estimates of the effect of changing sleep on frailty or a clinical intervention.

These examples are hypotheses for a later, source-validated locus analysis. Their choice was informed by observed results; the FI and FFS records are not independent confirmation, and the insomnia sources test sleep ascertainment rather than a joint causal mechanism. No common-variant colocalization, independently validated joint score, or per-person risk has been established.

## Five-source descriptive sensitivity screen

I checked the same fixed **15-row prototype** against all five existing marginal source results. This is a post-selection description, not a new hypothesis family or a validation test. Exactly two markers have p<0.05 in every source and the same aligned sign across all five:

| Variant, counted allele | FI p | FFS p | Jansen insomnia p | FinnGen insomnia p | MVP insomnia p |
|---|---:|---:|---:|---:|---:|
| rs12536335, G | 2.633×10⁻⁶ | 0.0017 | 7.651×10⁻⁷ | 0.00625446 | 4.57×10⁻⁸ |
| rs7734045, C | 0.005453 | 0.015 | 1.594×10⁻⁶ | 0.0352176 | 0.0184 |

For perspective, none of the 15 has p<0.05/15 in **every** one of the five sources. This is merely a stringent numerical benchmark, not a valid post-selection adjusted p-value or proof against shared association. The choice of markers used previously observed cross-trait statistics; FI and FFS samples overlap; exact sleep cohort overlap and source QC remain unresolved; SNPs can be correlated. Even the same sign in five studies does not identify a shared causal variant or rule out confounding. The two highlighted rows are candidates for a separately specified analysis and should not be labeled confirmed shared-risk variants. Their underlying betas, CIs and p-values are preserved in `prototype_variant_five_trait_evidence.csv`; the full 15-row numerical screen is in `prototype_five_source_descriptive_triage.csv`.

## Published FI and FFS polygenic-score file audit

I downloaded the complete original PGS Catalog scoring files and verified their published MD5 checksums before inspecting every data row. This is a check of the publicly distributed scoring files, **not** a re-analysis of Foote et al.’s original evaluation data.

| Catalog file | Rows | Effect and other allele both present and identical | Other allele blank |
|---|---:|---:|---:|
| [PGS005228, FI](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS005228/ScoringFiles/PGS005228.txt.gz) | 2,379,429 | 2,379,429 | 0 |
| [PGS005229, FFS](https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS005229/ScoringFiles/PGS005229.txt.gz) | 2,761,868 | 0 | 388,036 |

In the FI file, the `other_allele` field repeats the `effect_allele` on every row. That does **not** by itself show that the effect allele or weights are wrong, or that the published cohort results are invalid: a scoring system may operate on effect-allele dosage without an other-allele field. It does mean this file's two reported fields cannot be used as an ordinary two-distinct-allele harmonization pair. Resolve the source coding and allele matching independently before using it with a genotype dataset; do not silently infer the other allele. The FFS file has a blank other allele for 388,036 rows; a missing other allele can be permitted in score formats, but those rows require appropriate matching checks rather than an automatic allele-pair match.

Eight of our 15 prototype rsIDs occur in each score file, with all eight matching the prototype's reported GRCh37 coordinates. The other seven were not found by rsID; that is not a demonstration that their score weights are zero, because alternate identifiers or representation may be present. The matching rows and source weights are in `published_fi_ffs_pgs_candidate_lookup.csv`, and full file-count and checksum details are in `published_fi_ffs_pgs_source_audit.json`. These are reference lookups; no 15-marker subset was evaluated against outcomes, no published whole PGS was applied to individual genotypes, and the published whole-score performance cannot be attributed to the subset. A later technical computation of the separate 15-marker prototypes in reference genotypes is described in `reference_15_score_diagnostic.md`.

### Follow-up: FI score counted-allele check against the project FI input

To distinguish the duplicate `other_allele` field from the counted allele itself, I took a deterministic approximately 1-in-1,024 sample of PGS005228 using SHA-256 of each score row's rsID, chromosome and position, then compared the sample with the project's **harmonized** Atkins FI GWAS. Sampling was defined before the source comparison. The sample contained 2,323 score rows with unique rsIDs. Of these, 2,322 had a unique exact-rsID, chromosome and position match in the harmonized FI file, and the PGS `effect_allele` matched the FI input's `A2` for all 2,322. The remaining rsID, rs6912853, was absent from the filtered input and lies in the chromosome 6 MHC interval excluded by that input's documented QC.

This supports the plausibility of the **counted allele for this matched sample**; it does not validate the duplicated other-allele field, establish that every unsampled record is correct, independently verify the score weights or predictive performance, or resolve the score for variants absent after QC. The comparison is not independent—the score and FI GWAS are based on the same original source. The scoring file has not been modified or repaired. The full deterministic sample and outcomes are in `fi_pgs_allele_sample_audit.csv`, with the sampling rule, counts and source hashes in `fi_pgs_allele_sample_summary.json`.

### Full counted-allele comparison with the harmonized FI source

I extended the earlier deterministic sample to **all 2,379,429 rows** of the checksum-verified PGS005228 file and scanned all 6,291,635 records in the project's filtered, harmonized FI GWAS. Exact rsID plus GRCh37 chromosome and position matched 2,377,040 score rows; for every one of those matched rows, the PGS `effect_allele` was identical to the FI input's `A2`. There were no duplicate score rsIDs, coordinate disagreements, or counted-allele mismatches in the matched set. The remaining 2,389 score rsIDs were absent from the filtered FI input: 2,359 are in its prespecified chromosome 6 MHC exclusion interval and 30 are outside that interval, with their exclusion reasons not established here. Absence does not mean zero effect or a bad score variant.

The full-file counts independently agree with all 2,323 earlier deterministic-sample outcomes. This establishes strong **within-source counted-allele string consistency** for matched variants. It does not validate the repeated `other_allele` field, confirm every score weight or original reference/alternate representation, or establish independent prediction or calibration. The score was not modified; no unmatched allele was filled in automatically. Counts, checksums and algorithm scope are in `fi_pgs_full_counted_allele_audit.json`.

### Full FFS score counted-allele comparison

I also compared all **2,761,868** records in the checksum-verified PGS005229 FFS score with all **8,883,488** rows in the article-linked FFS summary-statistics file. Of the score records, **2,373,832** had usable score coordinates and a unique exact rsID, GRCh37 coordinate and biallelic allele-string match; every one had the score's counted allele equal to the FFS source `A1`. No counted-allele mismatches or coordinate disagreements were observed in that matched subset.

The remaining **388,036** score rows have no chromosome/position and no other allele in the original score file. As a weaker follow-up, each had a unique rsID match in the FFS source and the score counted allele matched the source `A1`. Without score-side coordinates this is **rsID-and-allele evidence only**, not an independently verified position match. Neither coordinates nor alleles were filled into the original score. The interpretation of `A1` as the source effect allele remains a documented-format inference. These checks do not validate weights, LD handling, genotyped samples, calibration or prediction.

The complete counts and checksums are in `ffs_pgs_full_counted_allele_audit.json`; the weaker rsID-only audit is `ffs_pgs_missing_coordinate_rsid_audit.json`. In conjunction with the FI audit above, the two checks show strong within-source consistency for the published scores while preserving their different missing-field limitations.

### Full-score weight-orientation diagnostics

I compared the published FI and FFS score weights with the corresponding **marginal GWAS beta aligned to each score's counted allele**, restricted to the exact rsID/GRCh37-position/allele matches in the full audits. This is an in-sample source-consistency diagnostic: the GWAS underlies the respective score. Both scores use SBayesR, so their jointly estimated, LD-aware weights need not equal the marginal GWAS betas or always have the same sign. The FFS source's A1 effect convention remains an inference.

| Score | Matched rows | Nonzero pairs with same sign | Nonzero pairs with opposite sign | Descriptive Pearson r: score weight vs. aligned marginal beta |
|---|---:|---:|---:|---:|
| FI PGS005228 | 2,377,040 | 2,139,292 | 199,961 | 0.898 |
| FFS PGS005229 | 2,373,832 | 2,107,307 | 218,043 | 0.730 |

The sign columns exclude pairs where either displayed score weight or source beta is zero. The files contain six-decimal score weights; 23,005 FI and 48,482 FFS matched weights are displayed as zero, which need not mean the underlying unrounded effect is exactly zero. Source beta rounding can also yield zeros. Correlations and sign percentages across linked variants are descriptive and do not provide independent predictive accuracy, a clinical effect, or an independent statistical comparison between FI and FFS. Differences in correlations can arise from source scale, LD, weights, filtering and SNP support. None of the 2,389 FI rows absent from the filtered source or 388,036 FFS rows without usable score-side coordinates is included in these calculations.

Complete counts, checksum-matched input provenance and descriptive source-p strata are recorded in `fi_pgs_weight_orientation_audit.json` and `ffs_pgs_weight_orientation_audit.json`. The FI calculation also reproduced the prior deterministic-sample sign counts. No weights or score files were modified, and no individual genotype or outcome was evaluated.
