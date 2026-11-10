# NICOLA and HALST: access and score-analysis plan

Prepared for an eligible research team, 1 November 2026. No cohort participant data have been accessed, no request has been submitted and no genetic score has been evaluated. Existing NHANES findings are a separate observational analysis.

## Decision required before analysis

Request approval and confirm whether linked sleep, frailty and genotype data can be analyzed in the same authorized environment. A public paper or genetic-only archive cannot establish that linkage. Freeze the exact 15- and 22-variant model files, allele conventions and weights before examining target associations; the 22-variant weighted models and 22-marker equal-count model are distinct and must not be silently interchanged.

## NICOLA request

Seek project-specific access through the [NICOLA research team](https://nicola.qub.ac.uk/FileStore/) and establish whether its secure setting supports an approved link to the [EGA genetic release EGAD00010002762](https://ega-archive.org/studies/EGAS00001007915). The genetic archive is not a substitute for the separate clinical data/linkage approval. NICOLA's policy describes possible cost recovery; no no-fee quote or permitted cross-environment process has been established.

Request the specific Wave 1 and follow-up sleep-question items, their dictionaries and visit dates; all 30 items and construction algorithm underlying the [published FI](https://academic.oup.com/eurpub/article/34/4/710/7632123); relevant weight, activity, grip and walking measures where available; demographic covariates; genotype dosages or calls, build, per-variant quality, batch and ancestry PCs; relatedness and permission to link all of these. Its 30-item FI includes daytime sleepiness. A derived FI omitting that item has at most 29 deficits, and further exclusions lower that maximum; whether it remains a valid FI, its completeness rule and component mapping need specific justification rather than assuming equivalence to the original instrument.

The [published methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC11293819/) state that the original 30-item score was calculated for respondents with complete data. For a paired item-removal sensitivity, first preserve that original complete-case population and compute both versions on the same people; report their denominators explicitly. A broader sample that becomes eligible after removing an item would be a separate analysis and cannot be treated as a matched contrast. The published categorical cutoffs have not been validated for a shortened score. The supplement and precise variable definitions still require verification before implementation.

The EGA landing page currently gives two inconsistent counts: 2,969 in its description and 2,978 in its dataset table. The [genetic resource preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC13174777/) reports 2,978 after initial QC. Obtain the current release-specific number from the custodian. None is the number with the requested outcomes and all variants. Ask for a confidential feasibility table: genotype–sleep–FI intersections by wave, age, consent and broad ancestry; missingness and outcome distributions; exact 15/22-variant coverage, allele contrast and imputation quality; overlap with discovery cohorts or family members; and total fees and secure-computing conditions. Do not infer zero overlap from different cohort names.

## HALST request

An authorized investigator should seek access under the [HALST cohort procedures](https://pmc.ncbi.nlm.nih.gov/articles/PMC5837206/) and confirm eligibility, IRB, cost and restrictions before submission. Request the individual-level sleep items and visit dates, the original Fried-type components and scoring code, covariates, direct or imputed genotype data, metadata and quality, ancestry PCs and relatedness; confirm permissions to link clinical and genomic data.

The [published study](https://ir.nhri.edu.tw/handle/3990099045/14916?locale=en-US) used 4,744 people with insomnia and frailty measures; this is not a verified genotyped intersection. Obtain the corresponding phenotype–genotype intersection and events, not merely separate counts. A sensitivity dropping the exhaustion criterion would produce a modified physical construct; report it separately from full Fried-type frailty and do not imply that it is a validated FI or an exhaustion-excluded GWAS. HALST may be informative for Taiwanese ancestry, but European-derived score performance and allele/LD portability cannot be assumed.

## Common analysis after authorization

1. Verify the analysis agreement permits the named study, data linkage, local environment and any AI-assisted handling; follow custodian rules. Create an audit of cohort participation and relatives across every original discovery and prior selection/follow-up dataset. Do not label a dataset independent without that audit.
2. Harmonize per variant using validated assembly, position, both alleles, strand, ploidy, frequency, imputation quality and source/target multiallelic contrast. Resolve source FFS effect direction and the existing count-versus-weighted allele differences. Prespecify missing-marker and proxy rules; an altered proxy score is a new model.
3. Select a primary sleep outcome and a primary frailty outcome before target association testing; report exact questions, recall period, clinical/symptom definition, eligibility, sample flow, overlap and missingness. Analyze full and item-excluded outcomes on the same eligible sample as a sensitivity, recognizing construct and denominator changes.
4. Evaluate each frozen score against sleep and frailty in specified separate models; compare with age, sex and ancestry-adjusted baseline models, accounting for relatedness, genotyping batch and study-specific design. Prespecify the multiple-testing family. Report effect estimates, uncertainty and incremental performance; only assess calibration/discrimination when the sample, endpoint and events support them. Avoid model tuning on the final test sample.
5. Report NIL/low-coverage results and selection, transportability and power limits. Genetic association or score prediction is neither evidence of causation nor proof of a causal gene. Observational NHANES and published HALST associations are not genetic replications.

## Suggested message to a mentor

> We identified NICOLA and HALST as possible additional validation cohorts for our experimental sleep–frailty genetic scores. Both have documented sleep, frailty and genetic resources, but their exact linked sample, marker coverage, access terms and cost need confirmation. Could we discuss whether your institution can help request a feasibility assessment or collaborate with either custodian? I have prepared a cohort-specific list of the data and checks needed.
