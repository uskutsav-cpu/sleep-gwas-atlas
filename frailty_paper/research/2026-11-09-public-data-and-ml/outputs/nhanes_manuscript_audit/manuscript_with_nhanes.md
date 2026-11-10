<!--
WORKING DRAFT — provisional Abstract, Introduction, Methods framework, Results
for frozen FI atlas reuse, secondary latent-factor analyses and read-only
aging-context extraction, plus a provisional Discussion.
Do not submit or treat as a complete manuscript. The systematic review has not
been screened, broader frailty-specific harmonization/analysis is incomplete,
and a complete manuscript-ready frailty_v1 package does not exist. The Methods text distinguishes the
completed protocol/search/plan locks from analyses that have not yet run.
Quantitative Results appear only when source tables are frozen under the
project plan; the current FI subsection reuses an immutable atlas-v1.0 extract.
-->

# Global genetic correlations of sleep and circadian phenotypes with frailty and aging

## Abstract (provisional)

**Background:** Sleep traits and frailty definitions are heterogeneous, and
prior genetic studies have reported sleep–frailty overlap without establishing
independent replication or mechanism.

**Methods:** We prepared a protocol-led review workflow and a staged genetic
analysis. PubMed searches were run on 22 September 2026; the protocol was
frozen before screening. Frozen atlas LDSC results were reused for 12 sleep or
circadian traits paired with the primary Frailty Index (FI), retaining the
original 396-pair correction family. A separate sensitivity analysis tested
the 12 traits against seven published latent frailty factors, with correction
across its fixed 84-pair family. A separate retrospective NHANES observational sensitivity analysis examined sleep-trouble disclosure and an adapted FI with direct sleep and fatigue items omitted. A separate FI × 12-sleep local-sharing
sensitivity family was locked before local-result access; it remains incomplete
and paused, and contributes no local estimates to this manuscript snapshot.

**Results:** Nine of 12 sleep–FI estimates met the inherited all-396
Benjamini–Hochberg criterion (q≤0.05). In the distinct latent-factor family,
48 of 84 estimates met its 84-pair BH criterion, and 37 also met the
same-family Bonferroni sensitivity. In NHANES, sleep-trouble disclosure was associated with the sleep/fatigue-excluded adapted FI (adjusted difference 0.0698, 95% CI 0.0528–0.0869; n=3,336). This observational result does not evaluate genetic overlap. The refreshed PubMed retrieval contains
61,009 source-record occurrences; after deduplication, 56,117 remain
unscreened. No independent
pairwise sleep–frailty genetic replication or participant-level score validation is available. Exploratory follow-ups do not establish a complete family-level local-sharing result or a validated shared causal locus.

**Interpretation:** The frozen results show genome-wide genetic correlation
between several sleep traits and the FI, while latent-factor findings remain
sensitivity evidence. Unknown participant overlap, European-ancestry source
dominance, unavailable matched sleep/fatigue-excluded GWAS statistics for this analysis and incomplete literature
screening limit interpretation. Genetic correlation does not establish
causality. This is a provisional results snapshot, not a completed review or
submission-ready manuscript.

## Introduction

Sleep is multidimensional. Self-reported insomnia, sleep duration, daytime
sleepiness, chronotype and sleep timing describe related but non-equivalent
features of sleep and circadian behavior. Genetic studies have examined these
phenotypes separately, including large studies of insomnia and sleep duration,
and studies using actigraphy-derived sleep measures [@jansen2019; @dashti2019;
@jones2019; @wang2019]. Combining them into a single sleep exposure can
therefore obscure phenotype-specific genetic relationships.

Frailty is also measured in more than one way. Deficit-accumulation frailty
indices, the physical frailty phenotype and electronic health record scores
capture overlapping but distinct aspects of vulnerability. Recent genetic
studies have extended this work from the Frailty Index and Fried phenotype to
the Hospital Frailty Risk Score and to latent factors that distinguish general
frailty from specific groups of deficits [@atkins2021; @ye2023; @mak2025;
@foote2025]. A 2026 CLSA study analyzed three-category Fried phenotype in
23,105 participants and reported one genome-wide significant lead association;
the authors reported no validation sample, effect estimates or confidence
intervals, and no SNP-heritability estimate [@borhan2026ordinal]. It therefore
adds phenotype-side context, but does not provide a matched genome-wide
sleep–frailty replication estimate. Estimates across these definitions should
be interpreted in the context of their source populations and phenotype
construction, rather than treated as interchangeable measurements of one
outcome.

Observational studies and prior systematic reviews have reported associations
between sleep problems and frailty, but have differed in which sleep features,
frailty definitions, age groups and study designs they included [@wai2020;
@insomniafrailty2025; @sleepqualityfrailty2025; @sleepdurationfrailty2025].
The reviews focused on older-adult populations and narrower combinations of
sleep and frailty measures. Their findings support the relevance of the topic
while leaving room for an updated evidence map spanning sleep and circadian
phenotypes, multiple frailty constructs and related aging outcomes. The
observational literature alone cannot establish shared genetic architecture or
causal direction.

Genetic association studies have already examined sleep–frailty relationships using both genome-wide correlation and Mendelian randomization (MR). Song et al. reported LDSC and local sharing for insomnia and the Atkins Frailty Index (FI), although their unconstrained genetic-correlation estimate exceeded the theoretical range and changed under a constrained-intercept analysis [@song2024sleepfrailtygenetics]. Deng et al. reported an LDSC correlation for chronotype with FI and evaluated multiple sleep traits by MR [@deng2024sleepfrailtymr]. Li et al. reported FI genetic correlations for sleep duration and daytime napping as well as MR estimates [@li2024sedentarysleepfrailty]. Other 2024 and subsequent studies examined sleep disturbances, FI and Fried frailty with MR, including multivariate sleep and aging outcomes [@lu2024sleepdisturbancefrailtymr; @che2025frailtysleepmr; @zhang2025sleepfrailtyagingmr; @gao2025sleepfrailtymr; @huang2026osafrailtymr; @zhang2026sleepagingfrailtymr]. MR estimates address a different estimand from LDSC genetic correlation and do not replicate rg. The 2025 Che et al. FI and sleep-duration meta-analyses both included UK Biobank participants, while exact participant intersections were not reported; source overlap for the other studies remains incompletely reconciled. The Huang et al. paper is present in the frozen PubMed snapshot as three query-source occurrences deduplicated to one record, and remains unscreened along with the rest of the queue. The complete focused comparison and unresolved source/overlap details are recorded in `original_analysis/prior_genetic_sleep_frailty_audit.md`.

The present study therefore does not claim the first genetic sleep–frailty
link, nor a new local signal or molecular mechanism. Its intended contribution
is a phenotype-specific evidence map across sleep and circadian traits,
frailty dimensions and aging outcomes, with replication and mechanistic
analyses restricted to data that pass source, overlap and quality-control
gates. Genome-wide genetic correlation can summarize shared effects across
traits, subject to the assumptions and limitations of the method
[@buliksullivan2015]; it does not identify direction, a causal variant, or a
molecular mechanism.

We are assembling a systematic evidence map and a staged genetic analysis of
sleep and circadian phenotypes in relation to frailty and selected aging
outcomes. The analysis distinguishes deficit-accumulation, physical and
health-record frailty measures; examines physical components and latent
frailty dimensions where source data and quality-control gates permit; and
prioritizes independent replication before local-sharing and functional
follow-up. Mechanistic analyses will be treated as exploratory and will proceed
only for loci that pass the prespecified upstream evidence gates. No causal
interpretation will be made from genetic correlation, fine-mapping,
colocalization or molecular-QTL overlap alone.

## Methods

<!-- PROVISIONAL METHODS FRAMEWORK. The planned procedures below come from the
     frozen v1 plan. They are not evidence that the corresponding analyses have
     been run. Reconcile with actual commands, versions, and frozen outputs
     before submission. -->

### Study design and analysis plan

This project combines a systematic evidence map with a staged analysis of
summary-statistic genetic associations. The analysis plan (frailty
genetics v1) was frozen before access to new frailty-specific association
results. It designates the 12-trait sleep and circadian panel in the original
atlas as fixed, with the Frailty Index (GWAS Catalog accession GCST90020053)
as the primary frailty endpoint. The Fried Frailty Score (GCST90295968) is a
separate construct-level replication endpoint, not an exact phenotype
replication; it will not be meta-analyzed with the Frailty Index. Other
frailty definitions, components, latent factors and aging phenotypes remain
secondary or contextual unless their source checks pass and a separate
multiplicity family is locked before result access. The executed latent-factor
family was separately frozen before its h² and rg analyses as a 12-sleep-trait
by seven-factor sensitivity family, with a fixed 84-pair correction
denominator; it is not a confirmatory replication family.

### Systematic review and evidence map

The protocol defines adults as the target population and includes observational
and human genetic-association studies of sleep or circadian phenotypes in
relation to a named frailty construct. A secondary evidence map covers physical
frailty components and related aging outcomes; these outcomes will not be
interpreted as equivalent to frailty. Three prespecified PubMed strategies
used a publication-date cutoff of 22 September 2026. After an interrupted
EFetch, acquisition resumed on 24 September; current ESearch IDs were
reconciled and the final neighborhood segment was refreshed. Exact search
strings, counts, retrieval logs, raw records, checksums and the preserved
pre-refresh log are retained in the review folder. The protocol was frozen
before screening, but after the PubMed search had been acquired; it was not
prospectively registered. Embase, Scopus, Web of Science and, if selected,
PsycINFO exports have not yet been supplied, and screening has not begun.

For the current PubMed-only record build, three prespecified query results
were combined with query and source-record identifiers retained. The tracked Python script
`frailty_paper/scripts/08_build_review_records.py` parsed 311 XML files (61,009 source-record occurrences, including the separately retained two-record index update) and deduplicated
by exact PMID, then normalized DOI where titles agreed or a unique title
candidate existed, then normalized title only when identifiers did not
conflict; otherwise records were retained separately. Duplicate memberships
remain in the audit trail, and the retained-record priority was primary
sleep–frailty, genetic sleep–frailty, then frailty-neighborhood query. This
build yielded 56,117 retained PubMed records and 4,892 duplicate occurrences.
These counts exclude Embase, Scopus, Web of Science and PsycINFO; licensed
exports have not been supplied, so cross-database merging and deduplication
remain pending. Two reviewers will screen titles/abstracts and full texts
independently; disagreements will be resolved by discussion or a third
reviewer and recorded. Full-text exclusions will receive a primary reason.
Observational risk of bias will be appraised with the JBI checklist appropriate
to each design, using item-level Yes, No, Unclear or Not applicable decisions
and evidence locators. Missing reporting will be judged Unclear rather than low
risk. Screening, extraction and appraisal have not begun.

### Exploratory NHANES phenotype-overlap sensitivity

Separately from all genetic analyses, we used public NHANES 2011–2014 data to examine whether the observational association of reported sleep trouble with a health-deficit index persists after removing directly overlapping sleep and fatigue items. The analysis was retrospective and not prospectively registered. The exposure was whether a participant had ever told a doctor or other health professional about trouble sleeping (SLQ050); it is not a clinical insomnia diagnosis. The outcome was an adapted 49-item NHANES frailty index based on Liu et al. (2025; https://doi.org/10.1186/s12877-025-05752-9), alongside a 47-item version omitting the two depression-questionnaire items on sleep or excessive sleep and tiredness/little energy. Each version was the mean of observed deficit scores and required at least 80% of its own items; the 47-item version required at least 38. We retained participants meeting all four FI versions' requirements for the paired comparisons. This adapted FI is not the Atkins FI or the original Fried physical phenotype.

We created the full survey design using positive examination weights (WTMEC2YR divided by two), cycle-specific masked strata and primary sampling units, then restricted to adults aged 60 years or older. Survey-weighted linear models adjusted for age, sex and cycle; the stated confidence intervals and p-values use 32 design degrees of freedom. A secondary complete-covariate analysis additionally adjusted for reported race/ethnicity, education, income-to-poverty ratio and smoking; this is not a causal adjustment set. Direct component omission cannot remove indirect overlap through health conditions, medicines, disability or care-seeking. Full item definitions, missingness, source provenance and sensitivity analyses are supplied in the NHANES supplementary audit.

### GWAS sources, phenotype hierarchy and cohort overlap

The source registry records each candidate GWAS accession, publication,
phenotype, file provenance, checksum, genome build, ancestry, sample size,
effect encoding and known cohort participation. Before harmonization, source
and release versions, variant conventions, access terms and study-specific
sample overlap will be checked. Unknown overlap will remain unknown; it will
not be treated as zero. UK Biobank, FinnGen, 23andMe, CHARGE and other named
cohorts will be recorded when supported by source evidence. Results from
related phenotypes will be labeled comparable-phenotype evidence, not exact
replication. A same-cohort result will not be called independent replication.
The completed secondary factor family used the seven European-ancestry,
GRCh37 Catalog GWAS GCST90624046–GCST90624052. The next accession,
GCST90624053, is a pneumonia GWAS rather than a frailty factor and was not
included; source scanning found 466,894 rows with missing or non-finite beta
or standard error. No eligible Fried, HFRS or physical-component sleep-pair
estimate is available in this package.

### Variant harmonization and quality control

Original source files will be preserved. New summary statistics will be
harmonized to the locked GRCh37/hg19 framework unless an individual source
requires a documented and validated liftover. The planned pipeline will
standardize variant identifiers and chromosomes; verify positions, effect and
other alleles, effect scales, standard errors, P values, allele frequencies
and per-variant sample sizes where available; detect duplicate or invalid
records; and conservatively handle strand-ambiguous variants. Primary locus
analyses require complete dense autosomal summary statistics, aligned alleles
and documented build. Lead-only and HapMap3-only files will not substitute for
dense locus inputs. Per-trait QC reports will preserve raw and retained row
counts, missingness, duplicates, ambiguous variants, warnings and source
metadata. The fixed 12-trait sleep panel passed source scanning,
harmonization, munging and the locked h² gates. The primary Frailty Index and
seven latent-factor sources also have recorded h² results; the 20-row master
table is an aggregation of those estimates, not a new LDSC run. Other
frailty-resource harmonization and QC remain incomplete or blocked.

### SNP heritability and global genetic correlation

For each eligible phenotype, SNP heritability and LDSC intercept diagnostics
will be estimated using ancestry-matched LD resources. The frozen input gate
requires a finite heritability Z statistic of at least 4 and a finite LDSC
intercept no greater than 1.20; failed or missing estimates will retain their
failure state rather than be represented as zero. The primary global analysis
is the fixed set of 12 sleep/circadian traits paired with the primary Frailty
Index. Genetic correlation, standard error, test statistic, P value, direction
and diagnostics will be reported. For integration with the original atlas,
multiple-testing interpretation will inherit its locked Benjamini–Hochberg
family of all 396 sleep-by-non-sleep pairs; q values will not be recalculated
on only the 12 displayed frailty comparisons. Any Fried-construct analysis
will remain separate and use its separately specified conservative family.
As a correction sensitivity, report Bonferroni-adjusted P values using that
same 396-pair denominator; this does not replace the inherited BH primary
decision.
The separate latent-factor analysis used two-sided LDSC for all 12 × 7
prespecified pairs only after all seven factor h² gates passed. Benjamini–
Hochberg q values were calculated across the full fixed family of 84, retaining
the sensitivity-only interpretation because UK Biobank participation is
present and exact participant overlap is unresolved. As an alternate
multiple-testing sensitivity, Bonferroni adjustment over that same 84-pair
family retained 37 of the 48 BH-significant estimates; the primary decisions
remain the locked BH q-values. The 72 selected
sleep–aging-context estimates in the Results are read-only extracts from the
frozen atlas; their original all-396 q values were retained and no new LDSC
analysis was run. These contextual rows are not replication.
The full Fried-score summary-statistics file has since been acquired and
structurally checked. It was not analyzed because the exact-release genome
build and effect-allele conventions remain unresolved, and the UK Biobank
sample overlap with most sleep sources is expected while exact intersections
are unknown. No construct-level replication estimate is reported.
Genetic correlation describes shared genome-wide effects and does not establish
causal direction [@buliksullivan2015].

### Replication, local sharing and shared-variant analyses

Replication will be classified as exact phenotype replication, comparable
phenotype replication or unavailable. Independence requires review of
ancestry, phenotype direction and cohort overlap. Distinct frailty constructs
will not be combined by effect-size meta-analysis. LAVA local genetic
correlation analyses are restricted to pairs passing the upstream
heritability and input gates, using the existing GRCh37 European-ancestry
reference with matching variant order and alleles [@werme2022]. A separate
sensitivity-only family was locked before frailty-specific local-result
access for the primary FI paired with each of the 12 fixed sleep traits across
all 2,495 regions in the repository locus file (29,940 pair–locus slots).
Inputs are the full dense harmonized summaries; the pinned LAVA runtime is
version 0.1.5 at commit `e729a245f7b6923967a96804fbf5246eadf2d6c6`, with the
GRCh37 European-ancestry UK Biobank LD v1.1 reference. The local-univariate
threshold remains the original atlas threshold, 0.05/(45×2,495) =
4.45335114673792×10^-7; it was not recalculated for this smaller family. A
local bivariate test is eligible only when both traits pass that threshold in
the same region. For each pair, the sample-overlap parameter is estimated from
the LDSC cross-trait covariance intercept divided by the square root of the
two univariate LDSC intercepts. Exact participant intersections remain
unknown, and this estimate may also reflect residual population structure or
other bias. Benjamini–Hochberg correction is applied across all 29,940 fixed
slots, assigning P=1 to ineligible slots and failed eligible tests; all
eligibility and failure states are retained. Family completeness requires no
more than 5% untested local-univariate slots, no more than 1% of loci with a
processing or per-trait univariate failure, and no more than 5% failed tests
among eligible bivariate slots. This family is sensitivity-only, not
independent replication or a confirmatory local analysis. A missing local
signal will not be interpreted as evidence of no sharing when QC or power
gates fail. Its execution and family-level QC must finish before local
estimates can be interpreted. The lock and prepared-input manifest are
recorded in `original_config/lava_frailty_sensitivity_v1.yaml` and the external
analysis workspace; the final run manifest and family-level tables are
pending completion. The read-only full-family audit at 23:07:06–23:07:58 UTC validated
18,234 of 29,940 receipts; inventory reached 18,237 and runner state reported 18,237 during
the scan. All 12 trait-level failure lower bounds exceed the frozen 1%
allowance. The coordinator requests six workers; the persisted resource
safeguard holds effective concurrency at two. Receipt, claim and duplicate
identity issues were zero, and both worker PIDs were alive on distinct jobs.
The run remains active to preserve the complete scheduled receipt set, but
the locked gate precludes local-sharing inference for every trait.

PLACO analyses require full dense harmonized inputs and prespecified eligible
pairs. Where the original 396-pair family applies, the locked criteria are
PLACO P no greater than 5×10^-8/396, marginal P no greater than 1×10^-4 and
z-squared no greater than 80; the separate conjFDR criterion is q no greater
than 0.05. Other testing families require their denominator to be fixed before
results are accessed. A PLACO association is statistical pleiotropy, not proof
of a shared causal variant or mechanism [@ray2020].

Fine-mapping will be undertaken only after the upstream families and a
nonduplicated locus union are complete. The locked contract specifies
SuSiE-RSS version 0.14.2, 95% credible sets, at most 10 signals, at least 50
variants per locus, and no locus truncation, with ancestry-matched signed LD,
validated allele order, convergence and diagnostics [@zou2022]. Trait–trait
colocalization will use coloc version 5.2.3 with primary priors p1=p2=10^-4 and
p12=10^-5 plus the locked sensitivity set. Strong shared-signal support
requires PP.H4 of at least 0.80 and a PP.H4:PP.H3 ratio of at least 5. Fine-
mapping and colocalization identify candidate variant and shared-signal
evidence under their assumptions; neither establishes causality or mediation
[@giambartolomei2014]. These analyses have not begun.

### Molecular context and reproducibility

QTL and cell-type resources will be queried only for loci that pass the
upstream evidence gates. Each analysis must record the resource, study,
ancestry, build, tissue or cell annotation, assay, overlap evidence and
limitations. Molecular overlap will be treated as contextual support rather
than independent replication unless cohort independence is demonstrated.
Analysis commands, software versions, configurations, input hashes and output
hashes will be recorded; final tables and figures will be generated only from
the completed, validated and frozen result package. Specific software versions
and executed commands are only partially inventoried in Supplementary Table
16. The separately documented FI h² verification used LDSC 3.0.1 with
checksum-pinned EUR LD-reference files; it reproduced FI h² at reported
precision but did not recreate the legacy atlas inputs or all analysis runs.
The full execution environment and command inventory therefore remain
incomplete; no submission-level reproducibility claim is made.

<!-- Remaining work: complete licensed-database searching and dual screening,
     obtain eligible alternative frailty sources, pass downstream evidence
     gates, reconcile software/commands, and revisit claims after those gates. -->

## Results

### Exploratory observational phenotype-overlap sensitivity

Among 3,472 older adults with positive examination weights, 3,336 had the sleep disclosure item and met the completeness rule for all compared index versions; 969 reported ever disclosing sleep trouble. The survey-weighted age-, sex- and cycle-adjusted difference in the adapted 47-item FI was 0.0698 (95% CI 0.0528–0.0869; p=1.54×10⁻⁹) for participants reporting sleep-trouble disclosure versus those not reporting it. The analogous difference for the full 49-item index was 0.0736 (95% CI 0.0566–0.0907). Index composition and denominators differed, so the difference between these coefficients is not a causal decomposition. In the common complete-covariate subsample (n=3,035), the more extensively adjusted difference was 0.0677 (95% CI 0.0512–0.0843). The independently rerun computation reproduced these figures; using the survey-regression residual-df convention instead of the reported design-df convention gave 0.0698 (95% CI 0.0527–0.0869; p=4.42×10⁻⁹).

These results reduce concern that direct sleep and tiredness items fully account for this **observational** relationship under this adapted NHANES measure. They do not quantify or correct phenotype overlap in the Atkins FI GWAS, establish genetic correlation or causality, or validate the genetic scores. Other health deficits can still reflect sleep-related conditions, fatigue, medication use and healthcare contact.

### Exploratory CLSA physical-component follow-up

After the CLSA public PheWeb2 service became available, on 6 November 2026 we screened the 33 distinct previously selected prototype variants against European-ancestry, combined-sex marginal associations for average grip strength, maximum grip strength and four-metre walking time. Existing reciprocal GRCh37-to-GRCh38 position mappings were matched to the target by exact position, allele pair and rsID. All 99 variant–phenotype pairs were returned. Reported variant-level sample sizes were 22,809 for each grip measure and 24,231 for walk time, rather than the slightly larger catalog totals. The public methods describe rank inverse normal transformation, so effects are not in kilograms or seconds. The frontend labels ALT as the effect allele; full source genotype coding was not independently reconstructed. Five displayed p values were below 0.05, but none survived Holm correction across the complete 99-test family (all adjusted p values were 1). Rounded source p values were used for correction rather than ratios of rounded effects and standard errors.

These physical measurements are not FI, FFS or a joint sleep–frailty endpoint; the two grip outcomes overlap. The screen was exploratory, and the rs113851554 response was inspected during endpoint discovery before the full screen. Participant and family independence from every discovery/selection source remains unverified. The results do not test score calibration, discrimination or aggregate-score association, and absence of corrected marker-level significance does not establish absence of smaller effects. The public catalog had no code or description containing “frail” or “insomn”; this does not establish that the separate ordinal-frailty GWAS is available or that all potentially related labels have been exhausted. Complete results and source receipts are supplied in the [CLSA supplement](../clsa_public_followup_20261106/README.md).

### Frozen Frailty Index estimates and workflow verification

The existing atlas-v1.0 analysis contains a frozen set of 12 global genetic
correlations between the sleep/circadian panel and the Frailty Index. In the
original locked multiple-testing family of 396 sleep-by-non-sleep pairs, 9 of
the 12 FI comparisons had Benjamini–Hochberg q≤0.05. The largest positive
estimates were for insomnia (r_g=0.6405, SE=0.0235), sleep apnea
(r_g=0.4687, SE=0.0321) and short sleep (r_g=0.4543, SE=0.0272). Sleep
duration (r_g=−0.2220, SE=0.0275) and sleep efficiency (r_g=−0.1520,
SE=0.0365) had negative estimates. Chronotype, accelerometer-derived sleep
duration and sleep timing did not meet q≤0.05. The complete 12-row result
extract and locked-family q-values are retained in
`original_analysis/frozen_atlas_frailty_global_rg.tsv`; q-values
were not recomputed on the displayed FI subset.
All nine comparisons meeting the inherited BH q≤0.05 criterion also met the
same-family Bonferroni sensitivity criterion (P×396≤0.05); the three other
comparisons met neither criterion. The row-level calculation is retained in
`original_analysis/frozen_fi_correction_sensitivity.tsv` with a checksum sidecar.

A separate `frailty_v1` workflow reran all 12 two-trait LDSC analyses. The raw
P values agreed with the frozen atlas at floating-point precision: 11 were
numerically identical and one differed only in decimal serialization
(relative difference 1.72×10^-16). The manuscript retains the frozen atlas
estimates and all-396 q values. This same-source workflow check is not
independent cohort replication.

The FI source is European ancestry, N=175,226. Its observed-scale SNP
heritability was h²=0.1093 (SE=0.0050; Z=21.86; LDSC intercept=1.020),
matching the fresh, checksum-pinned verification run at the reported
precision. UK Biobank contributes to the FI source and to 11 of the 12 sleep
sources, while exact participant intersections remain unknown. These global
correlations do not establish independent replication, local sharing, shared
loci, molecular mechanisms or causal direction.

### Secondary sleep × latent-frailty factors

The prespecified secondary family comprised 12 sleep traits × seven latent
frailty factors (84 pairs). All seven factor GWAS passed the inherited h²
input gate (finite Z≥4 and LDSC intercept≤1.20); Factor 3 had the largest
intercept (1.195). Across the 84 LDSC genetic correlations, 48 had
Benjamini–Hochberg q≤0.05 using the complete fixed family. The general factor
showed positive genetic correlation with insomnia (r_g=0.6658, SE=0.0256,
q=1.73×10^-147), short sleep (r_g=0.5278, SE=0.0319,
q=5.76×10^-60), and sleep apnea (r_g=0.4272, SE=0.0375,
q=8.53×10^-29), and negative correlation with sleep duration
(r_g=−0.2587, SE=0.0284, q=6.66×10^-19). Among the six specific factors,
the largest absolute estimate was sleep apnea with Factor 4 (r_g=0.4276,
SE=0.0419, q=2.19×10^-23). The insomnia–general-factor estimate has direct
construct overlap: insomnia is one of the 30 deficit GWAS used to construct
the general factor [@foote2025]. The project and factor studies used different
GWAS releases, but UK Biobank overlap is expected and the exact participant
intersection is unknown; this estimate is part-whole dependent and is not
independent validation. The factor model also includes tiredness/lethargy, a
related but nonidentical indicator to our daytime-sleepiness phenotype. All
latent-factor estimates remain sensitivity-only because factor GWAS contain
UK Biobank data and exact participant overlap with sleep GWAS remains unknown.
They do not constitute independent replication or identify causal effects.
The complete 12 × 7 estimate matrix with q-value markers is
shown in the [latent-factor heatmap](original_analysis/sleep_latent_frailty_rg_heatmap.pdf).
The [Figure 3 forest plot](original_analysis/figure3_latent_factor_forest.pdf)
presents the same sensitivity-family estimates with Wald 95% confidence
intervals; neither display tests differences between correlated factors.
Same-family Bonferroni sensitivity results are reported in
[Supplementary Table 15](supplementary_tables/table_s15_correction_sensitivity.tsv).
A separate descriptive [insomnia/sleep-apnea endpoint plot](original_analysis/insomnia_sleep_apnea_frailty_dimensions_forest.pdf)
shows the available FI and factor estimates with their distinct q-value
families; it is not a test of differences between endpoints.

### Read-only aging-context estimates from the frozen atlas

We extracted 72 previously estimated global genetic correlations covering
12 sleep traits and six selected aging traits from the frozen atlas, without
rerunning LDSC. The original all-396-pair q-values were retained; 26 of 72
comparisons had q≤0.05. Counts were 9/12 for parental lifespan, 8/12 for
healthspan, 7/12 for hand-grip strength, 1/12 each for longevity and
Parkinson disease, and 0/12 for Alzheimer disease. For example, insomnia was
positively correlated with healthspan (r_g=0.3432, q=3.70×10^-15) and
negatively correlated with parental lifespan (r_g=−0.2740, q=7.98×10^-13);
these directions follow each source phenotype's effect coding and should not
be read as causal effects. A separate source-level crosswalk classifies 44
pairs as shared-cohort-expected and 28 as possible/unknown; exact participant
intersection is unknown for all 72. This is contextual reuse, not independent
replication or a new frailty analysis.

## Discussion (provisional)

### Principal findings

In the frozen atlas, insomnia had the largest positive sleep–Frailty Index
genetic correlation among the 12 sleep phenotypes (r_g=0.6405); sleep apnea
and short sleep also had positive correlations, whereas continuous sleep
duration and actigraphy-derived efficiency had negative estimates. The
secondary latent-factor family showed a strong positive correlation between
insomnia and the general frailty factor and additional associations across
specific factors. The aging-context extraction also retained associations
with healthspan, parental lifespan and grip strength. Taken together, the
primary-FI estimates, secondary latent-factor sensitivity family and
read-only aging-context correlations are descriptive patterns of global
genetic correlation. They do not establish robustness across frailty
definitions or provide independent confirmation. They also do not show that
any sleep phenotype causes frailty or that one latent dimension explains more
of an association than another.

### Interpretation and evidence limits

The latent-factor estimates are useful for describing which specific factors
share genome-wide signal with sleep traits, but the current analyses do not
formally compare correlated genetic-correlation estimates. In particular,
counts of q-significant pairs across factors are not a test that one factor is
more strongly related to sleep than another. The cross-endpoint forest plot is
descriptive, combines results with distinct correction families, and should
not be interpreted as a test of endpoint differences. Independent assessment
requires nonoverlapping or appropriately overlap-adjusted data and a method
that supports direct comparison.

Independence is the main limit on interpretation. UK Biobank participation is
expected for the Frailty Index and most sleep sources, and the latent-factor
GWAS also include UK Biobank data. Exact participant intersections have not
been established. The aging-context estimates are read-only extracts from the
same frozen atlas and therefore provide context, not replication. We have not
established an independent sleep–frailty replication pair. At the 2026-09-26
23:07:58 UTC audit, the locked FI × 12-sleep LAVA sensitivity family remained
active and incomplete (18,234 of 29,940 receipts validated). Trait-level
failure lower bounds exceeded the frozen 1% allowance for all 12 pairs. The
complete scheduled run and family collation remain pending, and no local
estimates are interpreted. This does not establish an absence of local sharing.
There is no frailty-specific
shared-variant, fine-mapping, colocalization or molecular-QTL result; these
analyses remain gated on upstream evidence and eligible inputs.

The present evidence also covers only a subset of frailty constructs. Fried
frailty, HFRS and the five physical components have not contributed eligible
sleep-pair estimates in this package. The physical-component source provenance
and effect models remain insufficiently resolved, the acquired Fried file
remains gated from harmonization pending exact-release build/effect metadata,
and the exact HFRS endpoint remains unverified. Thus,
the latent-factor patterns cannot yet be used to claim which physical deficits
drive sleep–frailty overlap. The systematic-search acquisition is not a
completed review: licensed database exports, dual screening and risk-of-bias
assessment remain outstanding. These gaps prevent a complete synthesis of
genetic and observational evidence.

### Observational sensitivity and genetic inference

The NHANES sensitivity study is a useful separate test of item overlap in observed data, with an independently reproducible survey analysis and a transparent 47-item adaptation. The result is cross-sectional, retrospectively defined and vulnerable to health-selection and healthcare-disclosure effects. Its 80%-observed rule permits variable denominators. Mak et al. reported a GWAS of a modified FI excluding eight mental-wellbeing items including insomnia and fatigue (https://doi.org/10.1111/acel.13868); the full summary statistics for that modified phenotype have not been identified in their listed release. That broader exclusion cannot isolate the contribution of the two sleep and fatigue items. Obtaining those statistics and establishing sample comparability could enable a broader genetic sensitivity analysis; a same-participant GWAS excluding only sleep and fatigue would address the more specific question. NHANES does not provide either GWAS or independent score validation in the present analysis.

### Exploratory regional functional follow-up

Separately from the frozen FI analyses, an exploratory postselected FOXP2-region comparison used MVP insomnia and FinnGen hospital frailty summary statistics with external ancestry-oriented LD. In a one-megabase window, one modelled signal pair had PP.H4=0.979 under p1=p2=10⁻⁴ and p12=10⁻⁵; PP.H4 changed to 0.322 at p12=10⁻⁷. The narrower window omits a distal frailty signal, and wider analyses were sensitive to the LD reference. MVP sample size was approximated, and ancestry-matched in-sample LD and complete cohort independence were not established. This model-dependent statistical follow-up uses a distinct hospital-record phenotype and does not establish a shared causal variant, an effect on FOXP2, a relevant cell type, or score validation. It includes no molecular QTL in this signal-pair comparison. Detailed diagnostics and limitations are in the [regional audit](functional_followup/foxp2_gnomad_coloc_review.md).

### Published observational context in an independent setting

In a previously published study of 4,744 community-dwelling adults in Taiwan's Healthy Aging Longitudinal Study, Wu et al. used repeated-measures models and reported an adjusted association between insomnia and Fried-criteria frailty (odds ratio 1.41, 95% CI 1.16–1.72; https://doi.org/10.1016/j.jagp.2023.01.028). This is the authors' published observational result, not an analysis of our genetic scores. It concerns a physical frailty phenotype rather than the adapted NHANES or Atkins FI; its sleep and fatigue definitions and possible overlap require further comparison. The association cannot establish causation or genetic replication, and the published analytic sample does not establish how many participants have usable genotypes for this project.

### Conclusion

The completed analyses support phenotype-specific genome-wide genetic
correlations between sleep traits and the Frailty Index, with related signals
for a general frailty factor and several specific factors. The estimates are
not independent replication and do not establish causal or mechanistic
relationships. The current paper should remain a provisional analysis report
until independent replication, additional eligible frailty definitions,
systematic-review screening and the prespecified downstream evidence gates
are addressed.

### Figure legends

The working legends and current status for Figures 1–5 are maintained in
[`figure_legends.md`](figure_legends.md). Figures 1–3 remain provisional or partial;
Figures 4–5 are transparent eligibility-status panels, not local or molecular
findings. All five figures require final refresh after the remaining evidence
gates are resolved.

### NHANES supplementary materials

The accompanying `supplement/` directory contains the 49-item source and coding audit, survey model results, participant flow, item missingness by cycle and exposure, figure 1 and its source data, and the independent raw-GWAS score-source audit. The full original reconstruction package remains unchanged. See `audit_and_readiness.md` for source limitations, denominator and survey checks, and the genetic-score provenance assessment.

### Supplementary-table legends

The working legends and availability status for Supplementary Tables 1–16 are
maintained in [`table_legends.md`](table_legends.md). Only tables with source-
linked, currently available results have been generated; unavailable tables
are explicitly listed with their evidence gates.

### Supplementary material

The [provisional supplement](supplement.md) contains seven source-linked tables
and documents which requested tables remain unavailable or unjustified. It is
not the complete frozen supplement.

## Exploratory prediction of laboratory DNA effect magnitude

A separate molecular benchmark used 633,076 development observations, including 605,885 absolute log2 reporter-effect measurements and 27,191 auxiliary saturation coefficients with a separate output unit. Genomic neighborhoods and shared central sequence contexts were kept within partitions; neighborhoods overlapping the outside tests were excluded. Small neural regression heads were trained by backpropagation, with linear/tree baselines and a generic-cell comparison. The pretrained DNA encoder remained fixed. A documented development-only amendment added direct absolute-error training before outside model performance was calculated.

The primary outside test used 9,602 original fitted-effect measurements across primary neural progenitors and HEK293 cells, representing 4,820 verified physical DNA pairs from [McAfee et al. (2023)](https://doi.org/10.1016/j.xgen.2023.100404). Sequence-only mean absolute error was 0.0764 versus 0.0763 log2 units for the training-cell median in neural progenitors, and 0.1013 versus 0.0976 in HEK293 cells. Effect-magnitude Spearman correlations were −0.015 and 0.087, respectively. The additional fixed-model test used 547 baseline neuronal effects from [Retallick-Townsley et al. (2026)](https://doi.org/10.1038/s41467-026-77114-x), after duplicate reconciliation, physical-sequence verification and genomic-overlap exclusion. Error was 0.0829 versus 0.0773 for the training-cell median, with Spearman correlation 0.049. The paired genomic-window bootstrap confidence interval for error improvement over this median was entirely below zero. Thus useful cross-study magnitude prediction was not established. Original raw barcode QC and fitted-effect estimation were not independently reconstructed.

AVI, a permissive composite containing AlphaGenome-derived and other annotations, reduced macro error relative to sequence-only inputs by 2.71% in the primary study and 1.06% in the additional neuron study. The additional study failed the predeclared tolerance for reduced within-region ranking correlation; that statistic used only four eligible groups and does not establish biological harm. Both study-level gates were required, so AVI was not retained. The original primary result, which passed its individual gate, remains unchanged. These comparisons cannot attribute benefit uniquely to AlphaGenome.

Among 4,782 DNA pairs tested in both primary-study cell types, sequence-only cell-difference error was 0.1029, with correlation 0.018. A clearly labeled posthoc comparison found error 0.0909 when predicting no between-cell difference, versus 0.1029 and 0.1000 for sequence-only and AVI-added models. This comparison does not support learned variant-specific cell interactions and did not change the model or retention rule.

Nominal 90% intervals covered 97–99% of outside observations but were broad. The frozen input policy accepted none of the primary measurements because their cell groups lacked multi-study development support, and all 547 additional baseline neuronal measurements. Subsequent predeclared outside-performance protection withheld all cell groups. Coverage and abstention therefore do not establish useful precision or a universal coverage guarantee. Unsupported immune-stimulation conditions receive no authorized prediction; baseline projections compared with stimulated outcomes are diagnostic only.

The complete supplement preserves the original 1,351 gene-evidence rows and covers 498 distinct variants across five modeled cell contexts: 488 existing molecular candidates combined with 33 score markers, with 23 shared variants. There are 470 supported single-substitution inputs and 28 unsupported indels. Authorized prediction and rank fields are blank; separately labeled raw projections are unvalidated diagnostics. Nearby molecular credible-set comparisons remain unanswered. This work supplies neither person-level sleep/frailty prediction nor independent human genetic-score validation, genetic correlation, disease causation or causal target-gene evidence. The NHANES result remains separate observational evidence. [Outside-study error, coverage and answer availability](../sleep_frailty_ai_v7/figures/all_studies_error_and_answer_availability.png) are reported together; full source audits and adverse results are retained in the molecular supplement.


## Coverage-amended molecular prediction benchmark (V9)

A later, explicitly developmental amendment reused the already revealed McAfee and Retallick measurements, preserving the original V7 results above. The V9 magnitude dataset contained 613,161 measurements. One line-aware architecture, with either fixed or partially adapted pretrained representations, was compared using a single common epoch across seven whole-study omissions. The selected fixed-reader, five-epoch candidate improved development mean absolute error over a training-only cell/culture median by 0.26%; this selected development result is not independent evidence of accuracy. Three-seed final and source-omitted ensembles were trained with and without AVI using matched rows, epochs and shared initial weights. Independent numerical replay and fitted-gradient checks verified computation, not biological usefulness.

Before joining new outcomes, models, error predictors, empirical 90% interval calibration, acceptance rules and all new-study predictions were frozen by cryptographic hash. The new experiment comprised 610 eligible original absolute fitted log2 coefficients from [Matoba et al. (2026)](https://doi.org/10.1038/s41593-026-02454-2), after source, construct orientation, reciprocal genome-coordinate and overlap checks. The assay-active source subset includes detected and undetected allelic effects; there was no significance filter. Its repeated cultures derive from one primary neural-progenitor donor, whose independence from the earlier primary-cell donor is not established. Original barcode QC and coefficient fitting were not independently reconstructed, signed source-allele conventions remain unverified, and full pretrained genomic-neighborhood independence is not established.

Sequence-only absolute error was 0.190428, compared with 0.189128 for the training-only cell/culture median (-0.69% relative improvement; negative values mean worse prediction). Effect-magnitude Spearman correlation was -0.0170. The frozen full-model usefulness gate passed: **False**. Nominal 90% intervals covered 61.48% of new measurements; 26 cases passed the input policy before outcomes were joined, and 0 remained authorized after outside-performance protection. These results do not establish a universal confidence guarantee.

AVI changed error by -0.24% relative to the matched sequence-only model. All retention gates passed: **False**. Local confirmation used 12 nearby pairs in 9 one-megabase windows; sufficient local confirmation was **False**. Because AVI combines AlphaGenome-derived and other predictors, this experiment cannot attribute a benefit uniquely to AlphaGenome. The earlier V7 retention decision remains unchanged.

Direct paired-difference results on held NN-calibration regions compare predictions with zero, cell-median and culture-median differences. Those outcomes also informed source-omission interval scaling, so their coverage is descriptive. The new one-cell experiment cannot independently establish variant-specific differences between cell types. No project genomic context receives an authorized answer from this assay-active-only validation. This molecular work does not validate insomnia/frailty scores in people, establish genetic correlation or causation, or strengthen target-gene evidence by itself. The NHANES sensitivity result and all original genetic analyses remain unchanged. Detailed methods, complete adverse comparisons and readiness gaps are in the V9 molecular supplement.

## Cell matched molecular prediction and reserved experimental validation

A separate development expansion added 41,451 paired reporter measurements from 13,817 OL27 DNA designs in three cell models, preserving the earlier analyses and all original data partitions. Three frozen author-trained DNA models supplied features for a small regression head. One common epoch selected across seven omitted sources improved development error over training-only cell/culture medians by 0.77%; this is selected development evidence. All OL33-connected and independent blood-test regions remained outside fitting.

The reserved blood experiment comprised 236 designs representing 49 variants, with K562 as the primary known-line comparison; unfamiliar blood lines remained diagnostic. A second test comprised 39,723 OL33 measurements across three known lines from the same Siraj paper and lab. Full source DNA, genomic build and overlap checks preceded fitting, and models, confidence rules and both added-input arms were frozen before target join.

In the independent blood K562 test, sequence-only mean absolute error was 0.140531 versus 0.136613 for the training-only culture median, a relative improvement of -2.87%. Nominal 90% intervals covered 99.58%; 0 cases met the input policy before outcomes were joined. Full accuracy gates passed: **False**; useful confidence gates passed: **False**.

In the held OL33 library test, sequence-only mean absolute error was 0.171032 versus 0.175737 for the training-only culture median, a relative improvement of 2.68%. Nominal 90% intervals covered 94.86%; 0 cases met the input policy before outcomes were joined. Full accuracy gates passed: **False**; useful confidence gates passed: **False**.

The matched permissible composite AVI input met the prespecified retention gates: **False**. This comparison cannot attribute an effect uniquely to AlphaGenome. Complete negative and positive comparisons, interval size, cell/nearby differences and audit limitations are in the [cell matched molecular supplement](supplement/molecular_prediction_cell_matched.md). This laboratory benchmark supplies no human sleep/frailty score validation or causal gene evidence. All original genetic and NHANES observational results remain unchanged.

### Nearby changes and differences between cell lines

In the held library, the sequence model correctly ordered absolute measured effects in 52.73% of 1,119 nearby comparisons across 181 one-megabase windows. The bootstrap interval for ordering was 49.38% to 55.82%, including chance ordering. Its paired effect-difference error improvement over predicting no nearby difference was 0.51%, with an interval crossing zero. The independent blood publication supplied only two eligible nearby comparisons, which cannot establish local replication.

For 39,723 comparisons of the same physical designs between the three cell lines, ordering accuracy was 43.30%. The average correlation within fixed pairs of cell lines was -0.0153. Predicting no between-cell difference gave lower absolute error (0.175347) than the model (0.192371). Improvement over the less accurate training-median difference therefore does not establish variant-specific cell learning. These results do not support cell-specific prioritization of sleep or frailty variants.

Adding AVI increased independent blood prediction error by 0.54% relative to sequence alone; the paired window bootstrap interval for the error increase was 0.03% to 1.12%. The small held-library gain did not satisfy the predefined retention rule. The composite input remains excluded from authorized biological predictions.

The nominal 90% intervals achieved high coverage with a median halfwidth of about 0.45 log2 units, while both arms withheld all cases. This is conservative withholding without useful answer availability. The accuracy requirements were a minimum 5% error improvement over both training-only typical-effect baselines, positive lower confidence bounds and a within-cell-line correlation of at least 0.1 with a positive lower bound. The added-input requirements also included at least 1% improvement in both primary tests and preservation of nearby ranking.
