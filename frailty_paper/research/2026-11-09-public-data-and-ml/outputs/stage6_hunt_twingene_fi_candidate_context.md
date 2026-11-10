# HUNT insomnia and TwinGene frailty evidence for five score candidates

29 October 2026. Five of the 15 project score variants are in the [Lane et al. HUNT insomnia replication table](stage6_lane2019_candidate_replication_review.md) and have TwinGene frailty-index results in the project's existing [Atkins et al. cohort extract](candidate_twingene_fi_cohort_lookup.csv). This is a descriptive cross-trait context check, not a test of either full score.

| Variant | HUNT insomnia p | TwinGene FI p | Source UK Biobank FFS p |
|---|---:|---:|---:|
| rs113851554 | 0.000182 | 0.109 | 0.19 |
| rs1942262 | 0.00414 | 0.528 | 0.003 |
| rs4886140 | 0.0119 | 0.710 | 0.45 |
| rs35881094 | 0.177 | 0.274 | 0.20 |
| rs10865954 | 0.462 | 0.641 | 0.00000078 |

Three of these five have unadjusted HUNT insomnia p < 0.05; none has unadjusted TwinGene FI p < 0.05. rs113851554, near MEIS1, has the strongest HUNT insomnia evidence in this small set but no nominal TwinGene FI or UKB FFS evidence in the displayed analyses. Conversely, rs10865954 has a strong UKB FFS association but not a nominal HUNT insomnia association. The null results are imprecise and cannot establish absence of effects. The UKB FFS result contributes to the project's score development; the FI meta-analysis used for candidate selection includes TwinGene. HUNT individual-level overlap with every project cohort has not been exhaustively checked. Thus these p-values are not independent replication of shared risk, not a fair comparison of predictive scores, and not proof of a common causal gene or allele.

FI and FFS are distinct constructs. The FI includes health deficits and the discovery definition includes an insomnia item; the physical FFS uses a different construction. Each GWAS uses different measurement and power. HUNT measures self-reported insomnia symptoms, not frailty. The source TwinGene allele-direction adapter has a documented conflict with the original README; p-values are invariant to its sign, but do not use these rows to infer effect concordance. The [direction audit](twingene_direction_genomewide_audit.md), [coordinate audit](stage6_lane_aou_six_liftover_review.md), and prototype's multiallelic warnings explain further limitations. These variants and outcomes were selected post hoc; no multiple-testing-adjusted replication claim is made.

The [five-variant table](stage6_hunt_twingene_fi_candidate_context.csv) also includes UKB cohort and meta-FI p-values, HUNT odds ratios and Partners p-values. It is reproduced with `python3 work/cross_cohort/hunt_twingene_context.py`; the [receipt](stage6_hunt_twingene_fi_candidate_context_receipt.json) records input and output checksums. Source cohorts: [Lane et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6415688/), [Atkins et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8441299/), and [Ye et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/).
