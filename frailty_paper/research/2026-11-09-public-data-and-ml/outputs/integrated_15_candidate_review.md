# Fifteen candidate variants: integrated research review

This table brings together the project’s selected FI–insomnia candidate set with FFS, three insomnia sources, the full FinnGen hospital frailty GWAS and nearby-gene annotations. It is a **post-selected research shortlist**, not a validated risk score or 15 confirmed shared causal variants. All p-values are unadjusted, and small p-values do not measure personal risk.

The counted allele is the orientation used for effect sizes in the companion CSV. Nearby gene-body labels describe position, not proven target genes. FI contains UKB and TwinGene; FFS and UKB insomnia also include UK Biobank. FinnGen insomnia and HFRS use overlapping releases. The FFS effect-allele convention and TwinGene adapter remain inferred. Nine variants have multiallelic records in broader reference data despite matching the two strings in the specific source checks.

| Variant | Counted allele | Nearest gene-body label | FI p | FFS p | UKB insomnia p | FinnGen insomnia p | MVP insomnia p | FinnGen HFRS p |
|---|---|---|---:|---:|---:|---:|---:|---:|
| rs113851554 | T | MEIS1 | 0.199 | 0.19 | 2.35e-21 | 1.16e-10 | 4.22e-15 | 0.983 |
| rs8067974 | G | AC002117.1 | 7.77e-07 | 0.027 | 1.52e-11 | 0.223 | 0.261 | 0.777 |
| rs1909650 | G | PCDH8P1 | 6.32e-07 | 0.097 | 9.02e-10 | 0.102 | 1.58e-05 | 0.000386 |
| rs11152363 | A | TCF4 | 0.0871 | 0.032 | 2.18e-07 | 0.018 | 0.144 | 0.617 |
| rs1032524 | T | IPO9-AS1;NAV1 | 0.000353 | 0.0001 | 2.56e-06 | 0.408 | 0.161 | 0.77 |
| rs323509 | A | RP11-6N13.1 | 1.8e-06 | 8e-10 | 6.83e-09 | 0.000158 | 0.0509 | 0.00425 |
| rs78206187 | G | AC009474.1 | 0.000694 | 0.75 | 6.3e-08 | 0.000405 | 0.000319 | 0.275 |
| rs1942262 | A | TCF4 | 0.00594 | 0.003 | 2.57e-06 | 0.443 | 0.574 | 0.368 |
| rs4886140 | A | RPP40P2 | 7.89e-05 | 0.45 | 1.05e-07 | 0.166 | 0.0154 | 0.00226 |
| rs12536335 | G | FOXP2 | 2.63e-06 | 0.0017 | 7.65e-07 | 0.00625 | 4.57e-08 | 7.85e-09 |
| rs35881094 | G | LINC01122 | 0.602 | 0.2 | 4.37e-05 | 0.254 | 0.133 | 0.602 |
| rs17550821 | A | LINC01122 | 0.886 | 0.072 | 2.13e-08 | 0.247 | 0.413 | 0.632 |
| rs7734045 | C | AC009014.3 | 0.00545 | 0.015 | 1.59e-06 | 0.0352 | 0.0184 | 0.000325 |
| rs11790060 | C | RP11-165J3.6 | 0.00106 | 6e-05 | 2.2e-05 | 0.162 | 0.0627 | 0.17 |
| rs10865954 | T | KLHDC8B | 0.00181 | 7.8e-07 | 2.65e-06 | 0.0741 | 0.333 | 0.217 |

## How to use the shortlist

- rs12536335 and rs7734045 have matching estimated directions and nominal p < 0.05 in the original five reviewed sources and nominal HFRS associations. That is a descriptive screen after selection; FinnGen HFRS and insomnia overlap. rs12536335 also has opposite FI component directions under an inferred TwinGene adapter. Neither is a confirmed causal marker.
- rs323509 is strongly associated with FFS in the source study and has a nominal HFRS association; its MVP insomnia p=0.05089 and FI component directions disagree. It sits in a chromosome 5 region with allele-specific reference LD to a published FFS lead and an insomnia marker. LD is not colocalization.
- rs8067974 is an example where the positional label can miss another plausible target: a published frailty-domain SMR analysis links this top SNP to brain-cortex HEXIM1 expression, but the result does not establish shared insomnia–frailty causality.
- rs113851554 has strong insomnia evidence but weak FI and FFS marginal evidence in the original inputs; it is not a proven frailty variant. Nonsignificance elsewhere likewise does not prove absence of an effect.

The companion `integrated_15_candidate_review.csv` retains signed trait-specific effects, HFRS SEs, coordinates, directional and multiallelic flags, VEP annotations and nearby FFS-lead lookup. Effects are on different trait scales and cannot be summed together. The separate 15-marker prototype was evaluated only for arithmetic on public reference genotypes, without outcomes; it has not been calibrated or validated. See `reference_15_score_diagnostic.md`, `hfrs_insomnia_ldsc_intercept_audit.md` and the individual candidate follow-ups for fuller limitations.
