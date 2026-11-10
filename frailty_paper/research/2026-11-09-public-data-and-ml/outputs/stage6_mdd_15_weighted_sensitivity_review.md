# Fifteen-marker weights against case-enriched depression symptoms

Using the [case-enriched insomnia-symptom and fatigue GWAS screen](stage6_mdd_clinical_candidate_screen_review.md), we conducted an exploratory marginal-summary check of the full 15-marker prototype. This uses three different published prototype weight sets—FI, inferred FFS and insomnia—not a single jointly trained model. All 15 position-and-allele pairs were present in both GWAS. The 22-marker score was not evaluated because the cached full-reference genotypes for that list were not available in the current workspace.

| Outcome | FI weights Z | FFS weights Z | Insomnia weights Z |
|---|---:|---:|---:|
| Insomnia symptom | 0.48 | 0.03 | 1.09 |
| Fatigue symptom | −0.58 | −0.87 | −1.05 |

These approximations are all statistically uncertain; they do not establish zero effect, show that weights differ, rank FI versus FFS prediction, or validate either score. The outcomes are depression-related symptoms in case-enriched samples, not primary insomnia, FI, FFS or a joint sleep–frailty outcome. The target samples overlap across outcomes and may overlap prior project sources; non-overlap has not been established. Selection and repeated project-wide analyses make the nominal p-values exploratory and inappropriate as prespecified replication claims.

For each score, GWAS beta/SE is aligned to the prototype counted allele by the source A1/A2 pair. Let `a` be the per-allele score weights multiplied by the cached 503-person 1000 Genomes European genotype standard deviations, and `R` the corresponding reference dosage correlation (cross-chromosome correlations set to zero). The reported statistic is `aᵀz / sqrt(aᵀRa)`. This approximates a score association under strong assumptions about sample overlap across SNPs, comparable phenotype scaling, and LD and genotype variance between reference and target; it is not a regression on individual scores. The source GWAS includes varying cohorts and per-variant sample sizes. `SE × reference genotype SD` spans 0.0142–0.0164 for insomnia and 0.0222–0.0254 for fatigue; this is a rough precision diagnostic, not proof of the approximation. The FFS direction is inferred and multiallelic variants require source-specific checks. An ancestry-matched target LD matrix and actual participant data would provide stronger evidence.

Reproducibility: `work/stage6_mdd_15_weighted_sensitivity.py`, [summary statistics](stage6_mdd_15_weighted_sensitivity.json), and the prior [reference provenance](stage6_hfrs_sparse_score_summary_receipt_EUR503.json). The saved reference genotype file was generated for this project with alleles counted according to the 15-marker prototype. The GWAS source data remain subject to PGC scientific-use and redistribution terms, as documented in the source review.

A later [22-marker approximation](stage6_mdd_22_weighted_sensitivity_review.md) uses saved reference frequency and count-variance summaries rather than the missing individual-level reference cache. It is statistically uncertain and uses a different denominator convention, so it cannot rank the two prototypes.
