# Putting FI cohort differences in context

The UK Biobank–TwinGene difference estimates for rs28528997 and rs12536335 are among the larger differences observed in a deterministic sample of the released FI summary statistics. This is **context for follow-up**, not an empirical replication p-value or proof of biological heterogeneity.

For each of 97,901 sampled, matched variants, I calculated a standardized cohort difference `(beta_TwinGene − beta_UKB) / sqrt(SE_TwinGene² + SE_UKB²)` after aligning alleles and applying the inferred TwinGene sign adapter. About **6.1%** had an absolute difference statistic greater than 1.96. Among the **1,155** sampled variants with |FI meta beta / SE| ≥ 3, the corresponding proportion was **10.6%**. These are descriptive proportions of correlated variants, not independent tests or estimates of the rate of true between-cohort differences. The samples, phenotype definitions, allele representations and standard-error assumptions can all affect them.

| Candidate | Absolute standardized FI cohort difference | Proportion of 1,155 stronger-FI sampled variants with difference at least this large |
|---|---:|---:|
| rs28528997 | 4.35 | 1 / 1,155 (0.087%) |
| rs12536335 | 3.89 | 2 / 1,155 (0.17%) |
| rs2815764 | 2.81 | 21 / 1,155 (1.82%) |
| rs1909650 | 1.43 | 300 / 1,155 (26.0%) |
| rs817509 | 1.22 | 391 / 1,155 (33.9%) |
| rs2155292 | 0.40 | 912 / 1,155 (79.0%) |

The |meta beta / SE| ≥ 3 group was chosen to give context among stronger FI associations; it does **not** reproduce the candidates' actual cross-trait selection. The six candidates were chosen using other association results and have been repeatedly inspected. Sampling SNPs does not make their LD-independent loci independent; sampled variants can also be correlated with the candidates. These proportions therefore must not be interpreted as corrected p-values, discovery probabilities or formal evidence for candidate outliers. The group statistic uses rounded meta betas and standard errors. UKB and TwinGene FIs were separately normalized and differ in composition, age and ascertainment, and the zero-covariance assumption is unverified. The inferred sign adapter is supported by within-file reconstruction but conflicts with the source README and remains unconfirmed by the authors. Multiallelic representation adds uncertainty for five of six candidates.

The appropriate next checks are to confirm the original cohort allele convention and representation, compare well-matched phenotypes in independent cohorts, and investigate whether residual technical or measurement differences explain the estimates. This analysis cannot establish the effect of removing sleep items or validate a risk score. The earlier `six_leads_fi_cohort_comparison.md` gives the betas, confidence intervals and direct difference tests.

## Sampling and verification

We reused the deterministic rsID-hash sample from `twingene_direction_genomewide_audit.md`, matched the 6,291,635-row filtered FI meta file and 16,557,800-row cohort archive on rsID, chromosome, position and allele pair, and excluded the same eligible duplicates before additional numerical checks. The matched count exactly reproduced the previous audit's 97,901; no additional rows were excluded by the meta-SE and reconstruction checks. `fi_cohort_heterogeneity_reference_context.json` provides the full distribution summaries and all candidate values. The sample is drawn from the filtered available set, not the entire genome; no LD pruning or independent outcome data were used.
