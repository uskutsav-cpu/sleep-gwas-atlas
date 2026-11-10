# CLSA public physical-component follow-up — 6 November 2026

The previously unavailable Canadian browsers are online. We queried the public CLSA PheWeb2 API for the project's 33 distinct candidate variants against three measured physical traits. All 99 requested variant–phenotype combinations were present. Five reported p values were below 0.05; none survived Holm correction across all 99 tests (all adjusted values were 1). This exploratory screen provides no multiplicity-adjusted evidence at these individual markers for these outcomes. It does not establish absence of smaller associations, independent frailty replication, a shared causal mechanism, or predictive validity of either score.

## Sources and phenotype inventory

- [CLSA PheWeb2](https://clsa-pheweb.cerc-genomic-medicine.ca/) and its public [phenotype catalog API](https://clsa-pheweb.cerc-genomic-medicine.ca/api/v1/phenotypes/), accessed 6 November 2026. The response has 4,260 rows representing 745 distinct phenotype codes across ancestry/sex strata. No code or phenotype description contains `frail` or `insomn`. This lexical catalog check is not proof that a related construct is absent under every possible label. It does not establish that the separate Borhan ordinal frailty GWAS has been released here.
- [CLSA about page](https://clsa-pheweb.cerc-genomic-medicine.ca/about). Its public application source states that genotype coordinates were converted from GRCh37 to GRCh38, continuous traits underwent rank inverse normal transformation, and Regenie models included genetic PCs, genotyping batch, sex, age and age squared. Therefore the displayed effects are not kilograms or seconds. The frontend explicitly passes `alt` as the table's effect allele. We used that documented display convention; the upstream processed-file conversion has not independently been reproduced. The raw frontend and catalog are archived with this follow-up.
- [CARTaGENE phenotype catalog](https://cartagene-pheweb.cerc-genomic-medicine.ca/api/phenotypes.json) was also obtained: 701 entries, with sleep duration, apnea and walking-activity labels but no `frail` or `insomn` label. These are different outcomes; we did not substitute them for the project's clinical insomnia or FI/FFS outcomes. [CARTaGENE's own documentation](https://cartagene-pheweb.cerc-genomic-medicine.ca/about) says its analysis is European-ancestry, uses Regenie, and applies rank inverse normal transformation to continuous traits. Its data and results should not be conflated with CLSA.

## Target outcomes

| Code | Measurement | Catalog N | Variant-level N reported for these candidates |
|---|---|---:|---:|
| GS_EXAM_AVG_COM | Average grip strength across trials | 22,817 | 22,809 |
| GS_EXAM_MAX_COM | Maximum grip strength across trials | 22,817 | 22,809 |
| WLK_TIME_COM | Time to complete a 4-metre walk | 24,242 | 24,231 |

We selected European-ancestry, combined-sex summaries. These counts are not interchangeable with an age-60+ sample, nor are they a joint sleep–frailty participant count. Continuous grip and walking are components related to physical frailty; they are not either frailty index, the full Fried phenotype, or dichotomous muscle weakness. Higher transformed grip means greater strength; higher transformed walk time means slower walking. Average and maximum grip overlap and must not be counted as independent cohorts or independent replications. We did not fit or reconstruct a combined phenotype from marginal statistics.

## Results and precision

| Variant | Outcome | ALT allele | Displayed ALT effect | Displayed p | Holm p, family of 99 |
|---|---|---|---:|---:|---:|
| rs113851554 | Average grip | T | −0.031 | 0.014 | 1 |
| rs113851554 | Maximum grip | T | −0.028 | 0.024 | 1 |
| rs323509 | Maximum grip | C | 0.014 | 0.025 | 1 |
| rs28528997 | Average grip | G | 0.015 | 0.037 | 1 |
| rs323509 | Average grip | C | 0.012 | 0.046 | 1 |

The ALT allele entries above must be read against the exact response and association table; the automated manifest retains the returned reference/alternate alleles. The source reports rounded effects, standard errors, p values, frequencies and imputation quality, limiting numerical precision. We used reported p values for correction; rounded effect/SE ratios need not reproduce them. Displayed imputation-quality values range from 0.93 to 1.0; this is source summary metadata, not a new genotype-level QC audit. None of the 99 tests passes a Bonferroni threshold of 0.05/99 either. The smallest displayed p is 0.014, far from that threshold, so the lack of significance is not a borderline rounding decision.

All results, including negative results, are in `candidate_component_associations.csv`. The 33 candidates came from earlier analyses; this is not a prospectively preregistered replication. During endpoint discovery we inspected the rs113851554 response before executing the complete screen. This selection history must be retained. Participant and family overlap with every discovery/selection source, including prior gait follow-up sources, has not been exhaustively established; do not label the associations independent validation merely because the cohort name differs.

## Alleles and score limits

The existing project crosswalk provides unique reciprocal GRCh37→GRCh38 position mappings. The public API accepted one exact coordinate and allele-pair representation with the expected rsID for each of the 33 variants. These checks support identity of the queried target contrasts, not equivalence of every source genotype model. The separate `score_allele_aligned_component_associations.csv` reports marginal effects oriented to each prototype's own counted allele, after checking its exact mapped pair. The 22-marker equal-count and weighted definitions have different counted alleles at eight sites and are kept separate. A sign flip assumes the same biallelic contrast; multiallelic source encoding and FFS source transformation caveats remain. We did not calculate a patient score, refit weights, impute missing phenotypes, or compare predictive performance.

For example, the rs113851554 T allele has a nominal association with lower grip in these displayed summaries. It is not significant after the planned screen correction and does not establish a sleep–frailty mechanism or causal gene. The absence of corrected marker-level results cannot be interpreted as a test of the aggregate score or as a powered equivalence test.

## Reproduction

Run `python3 work/clsa_public_candidate_screen_20261106.py` from the project root. It uses Python's standard library and the system `curl`, preserves TLS verification, caches responses, and requires the project's archived GRCh38 variant crosswalk and three distinct score manifests. Four concurrent requests are used; the cache prevents repeated downloads. The script writes 99 target rows, aligned rows for each score definition, phenotype metadata and an audit receipt. It rejects identity conflicts and does not treat unavailable entries as zero effects or null p values. The browser source and catalog snapshots are included under `source_metadata/`.

The script requires the larger project inputs, so this directory is a source-and-result follow-up rather than a self-contained genotype analysis environment. `receipt.json` records response URLs and hashes; `source_metadata/checksums.json` records metadata snapshots. No restricted or participant-level data were accessed.

## Manuscript use

Suggested supplementary wording: “After public service became available, we screened 33 previously selected variants against European-ancestry CLSA summary associations for average grip strength, maximum grip strength, and four-metre walking time. All 99 variant–phenotype pairs were retrieved with exact coordinate, allele-pair and rsID matches. Five displayed p values were below 0.05, but no association survived Holm correction across 99 tests. These exploratory marginal association lookups used rounded public summaries and do not validate the prototype scores or a composite frailty outcome.”

NHANES remains a separate observational sensitivity analysis. The missing sleep/exhaustion-excluded frailty GWAS and independent individual-level score evaluation remain the highest-priority genetic gaps.
