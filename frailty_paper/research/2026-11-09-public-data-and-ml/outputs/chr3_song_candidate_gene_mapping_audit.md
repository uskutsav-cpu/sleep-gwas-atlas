# The chromosome 3 shared-risk lead: interpreting competing gene labels

**rs34290943 is a candidate marker, not yet an identified causal gene.** [Song et al.](https://www.frontiersin.org/journals/aging-neuroscience/articles/10.3389/fnagi.2024.1358996/full) mapped it to **DALRD3** in their cross-trait analysis. Our Ensembl GRCh37 positional lookup labels its nearest overlapping gene body **LAMB2**, and the allele-specific transcript lookup reports a synonymous annotation in two LAMB2 transcripts along with other transcript annotations. These are distinct annotation statements; the sources inspected do not establish why Song's mapping selected DALRD3 or which gene mediates either phenotype. A synonymous annotation alone does not establish a functional effect.

Song reported a second candidate, rs10865954, mapped to MIR191. Our lookup placed that marker within KLHDC8B; its C alternate also has transcript annotations involving KLHDC8B and C3orf84. The author mapping, gene-body overlap and VEP consequence annotations use different criteria or annotation resources, and their outputs should not be presented as interchangeable. Do not infer the original mapping rule without its code or sufficient metadata. [Our prior annotation audit](candidate_allele_consequence_review.md).

## Association results to keep separate

Song's Table 1 reports two-trait MTAG p-values of 5.81 × 10⁻¹⁰ for frailty and 5.63 × 10⁻¹⁰ for insomnia at rs34290943, and CPASSOC p = 7.55 × 10⁻⁹. These are **cross-trait method outputs**, not the single-trait p-values in the original source GWAS, not proof of a common causal variant, and not risk estimates for either allele. Song's reported effect allele is C. The project’s later shared-count orientation is T, the other allele; a p-value does not contain the effect direction.

For the T allele, the project’s aligned follow-up gives the following single-trait p-values. These are distinct models and scales; several datasets share or may share participants. The comparisons and the 22-marker candidate set were selected after prior analyses.

| Result source | Trait | Source p |
|---|---|---:|
| Atkins | FI | 2.60 × 10⁻⁵ |
| Ye | Fried score | 1.10 × 10⁻⁴ |
| Mak | FI | 1.10 × 10⁻⁸ |
| Mak | modified physical frailty | 1.98 × 10⁻⁴ |
| FinnGen | hospital frailty score | 0.00401 |
| UK Biobank | insomnia symptom measure | 0.000719 |
| MVP | clinical insomnia | 0.000485 |
| FinnGen | clinical insomnia | 0.00110 |

All eight aligned point estimates are positive, but the nominal, unadjusted p-values are not eight independent replications, and the original FI includes sleep-related content. The case-control and quantitative coefficients cannot be summed; their [aligned effects and uncertainty](shared22_eight_source_evidence_matrix.csv) should be interpreted separately. The source build of the Song table is assumed to be GRCh37, and exact source alleles and multiallelic representations remain under review.

For Song's other candidate rs10865954, the separately obtained MVP-insomnia and Finnish HFRS marginal p-values are **0.333 and 0.217** respectively; neither is nominally significant. This is not proof of no association or a contradiction of the published cross-trait result: the analysis methods, cohorts and phenotypes differ, and power and allele representations matter. Its presence in the earlier 15-marker illustration therefore does not make it a member of the later 22-marker joint nominal screen. [Source rows](song124_mvp_hfrs_joint_screen.csv).

The two markers are approximately 50 kb apart. An rsID-level 1000 Genomes EUR lookup reported r² ≈ 0.061, but this is not an allele-specific, conditional or causal analysis and rs10865954 has a multiallelic reference record. It does not establish that the two signals are independent or share a causal variant. [Reference lookup](selected_variant_pair_ld_reference.csv).

## Interpretation and next step

Treat the locus and all gene labels as hypotheses. Verify each GWAS build and allele pair, characterize ancestry-appropriate regional LD, and use suitable full summary statistics to fine-map each trait and test shared versus distinct signals. Only then use relevant molecular-QTL and regulatory evidence to prioritize genes; a failure to reject heterogeneity in a gene-mapping method is not by itself proof of one causal variant or gene. The project's local genetic-correlation runs did not pass quality gates, so Song's reported 3p21.31 local result has not been independently verified by those runs. No gene function or individual risk is established by this audit.
