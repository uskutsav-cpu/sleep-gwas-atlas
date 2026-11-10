# Which allele points toward higher FI, FFS and insomnia?

`candidate_trait_direction_alleles.csv` translates the three existing prototype betas into the allele associated with a higher value of each trait **within the reported two-allele contrast**. If a beta for the originally counted allele is positive, that allele is listed; if negative, its reported comparator is listed. This does not change any of the original score weights or create an individual risk probability. For multiallelic sites, reversing the counted allele as if all genotypes were biallelic is not safe until study-specific variant representation and genotype coding have been verified.

Of 15 post-selected markers, **12 have the same estimated higher-direction allele across the FI meta-analysis, FFS and original UKB insomnia inputs**. Three do not:

| Variant | Higher FI direction | Higher FFS direction | Higher insomnia direction | Nearby gene-body label |
|---|---|---|---|---|
| rs113851554 | T | G | T | MEIS1 |
| rs35881094 | G | T | G | LINC01122 |
| rs17550821 | G | A | G | LINC01122 |

These are point-estimate directions, not statistically established opposite effects. In particular, weak or imprecise estimates can change sign. All three variants fail at least one nominal p < 0.05 test in the three-source table. Both rs113851554 and rs35881094 have multiallelic reference records; source coding requires further verification.

Eight of the 12 same-direction markers have unadjusted p < 0.05 in **all three** original inputs. This is a descriptive, post-selection filter, not a corrected discovery or independent replication. In the five-source follow-up, rs12536335 and rs7734045 had the same sign and nominal p < 0.05 across FI, FFS and the three reviewed insomnia sources. That additional pattern is still subject to overlapping cohorts, multiple testing and selection. For rs12536335, the FI UK Biobank and TwinGene component directions disagree under the project's inferred TwinGene adapter. The nearby FOXP2 and AC009014.3 annotations for those two markers are positional, not evidence of causality.

Do not combine the three raw weights into a single score: the FI, FFS and insomnia effects describe different measurements and scales. A usable score for a specific outcome needs frozen marker selection, verified allele coding, LD-aware weights and evaluation in independent genotyped participants. The FFS allele convention and TwinGene adapter remain inferred, nine of the 15 rsIDs are multiallelic in the reference, and the project’s prespecified local-correlation QC gates remain failed. The original inputs, p-values, alleles and caveats are preserved in `candidate_risk_gene_review_15.csv` and `prototype_five_source_descriptive_triage.csv`.
