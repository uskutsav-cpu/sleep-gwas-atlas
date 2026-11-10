# Fifteen insomnia–frailty candidates: combined evidence

This is a research review table, not an individual risk score. The candidates were selected from the top 30 published cross-trait results after availability filtering, not from a prespecified independent screen. Each p-value is unadjusted. A result for one phenotype does not establish a shared causal variant; lack of significance does not rule one out. Nearest gene is a positional annotation only.

| Variant | Counted allele | Nearest gene body | FI meta p | FFS p | UKB insomnia p | TwinGene FI p | FI cohort signs |
|---|---|---|---:|---:|---:|---:|---|
| rs113851554 | T | MEIS1 | 0.199 | 0.19 | 2.35e-21 | 0.109 | same |
| rs8067974 | G | AC002117.1 | 7.77e-07 | 0.027 | 1.52e-11 | 0.643 | same |
| rs1909650 | G | PCDH8P1 | 6.32e-07 | 0.097 | 9.02e-10 | 0.842 | opposite |
| rs11152363 | A | TCF4 | 0.0871 | 0.032 | 2.18e-07 | 0.651 | opposite |
| rs1032524 | T | IPO9-AS1, NAV1 | 0.000353 | 0.0001 | 2.56e-06 | 0.675 | same |
| rs323509 | A | RP11-6N13.1 | 1.8e-06 | 8e-10 | 6.83e-09 | 0.0775 | opposite |
| rs78206187 | G | AC009474.1 | 0.000694 | 0.75 | 6.3e-08 | 0.0321 | same |
| rs1942262 | A | TCF4 | 0.00594 | 0.003 | 2.57e-06 | 0.528 | opposite |
| rs4886140 | A | RPP40P2 | 7.89e-05 | 0.45 | 1.05e-07 | 0.71 | opposite |
| rs12536335 | G | FOXP2 | 2.63e-06 | 0.0017 | 7.65e-07 | 0.007 | opposite |
| rs35881094 | G | LINC01122 | 0.602 | 0.2 | 4.37e-05 | 0.274 | opposite |
| rs17550821 | A | LINC01122 | 0.886 | 0.072 | 2.13e-08 | 0.492 | opposite |
| rs7734045 | C | AC009014.3 | 0.00545 | 0.015 | 1.59e-06 | 0.183 | same |
| rs11790060 | C | RP11-165J3.6 | 0.00106 | 6e-05 | 2.2e-05 | 0.973 | same |
| rs10865954 | T | KLHDC8B | 0.00181 | 7.8e-07 | 2.65e-06 | 0.641 | opposite |

The FI meta-analysis includes the UKB and TwinGene FI components displayed in the companion CSV; those columns are not independent replications. FI and FFS both contain UKB participants, and other sample overlap is not fully known. TwinGene sign is based on the project’s audited sign adapter, not an author-confirmed correction; “same” and “opposite” compare raw estimated directions and do not test a difference statistically. Many TwinGene estimates have wide intervals. FFS A1 effect coding remains an informed source-format inference; nine variants have more than two reference alleles and need exact-contrast verification.

For context, rs12536335 lies within the FOXP2 gene body and showed nominal evidence in five previously reviewed sources, but its FI component directions are opposite; it is a follow-up candidate, not a confirmed shared causal variant. rs113851554 lies within MEIS1 and has strong insomnia p-values, but FI and FFS evidence is less clear in these inputs. These gene locations do not show mediation or causal function.

The CSV includes aligned effects and standard errors where available, all three insomnia sources, gene-annotation query URLs, and listed counted alleles. Source-study provenance and unresolved exact allele contrasts are documented in the companion notes. The two frailty beta scales and the insomnia beta scales differ; do not add them or interpret p-values as risk probabilities. The existing sparse prototype has not been fitted or evaluated in a cohort. The prior TwinGene LDSC precision gate and the FI–sleep LAVA QC gates remain failed.

See `fi_ffs_component_crosswalk.md` for score definitions; `candidate_twingene_sensitivity.md` for the TwinGene adapter; and `candidate_nearby_gene_notes.md` for the gene-distance method.

Allele-specific VEP consequence results, including the NAV1 3′ UTR annotation for the rs1032524 C alternate and FOXP2 intronic annotation for rs12536335 G, are documented in `candidate_allele_consequence_review.md`. These annotations do not establish shared disease mechanisms or score validity.
