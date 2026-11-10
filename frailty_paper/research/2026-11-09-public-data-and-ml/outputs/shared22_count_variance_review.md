# How regional marker correlations affect the 22-marker count

In 503 European-labeled 1000 Genomes reference individuals, the 22-marker shared-direction count has sample variance **9.223 allele-count units squared**. A descriptive covariance decomposition attributes 1.203 units², or **13.0% of that sample variance**, to the five marker pairs on the same chromosome within 1 Mb. This is a property of this constructed count in this reference sample—not variance in insomnia or frailty explained, predictive accuracy, or the effect of removing markers.

| Contribution to count variance | Units² | Share of total |
|---|---:|---:|
| Sum of the 22 individual dosage variances | 8.025 | 87.0% |
| Twice the covariance of the five pairs within 1 Mb | 1.203 | 13.0% |
| Twice the covariance summed over the other 226 pairs | −0.005 | −0.1% |
| Total | 9.223 | 100.0% |

Rounding affects the displayed percentages. The five regional-pair contributions are positive in this sample. The largest are rs7734045–rs7731706 (0.391 units²) and the FOXP2-labeled rs12536335–rs12671330 pair (0.362 units²). These components are **algebraic allocations**, not independent causal effects. The sum over distant pairs is near zero because positive and negative covariances can cancel; it does not mean every such pair has zero correlation. Sample structure, relatedness and sampling noise can affect pooled covariances, especially across distant loci.

For count `S = sum(X_i)`, the calculation uses `Var(S) = sum Var(X_i) + 2 sum_{i<j} Cov(X_i,X_j)`. It uses sample covariance of each person's counted-allele dosages. The separate [regional LD review](shared22_reference_ld_interpretation.md) reports phased-haplotype r and r²; these are related but different quantities. The previous regional LD analysis considered only same-chromosome pairs within 1 Mb; this variance decomposition includes all 231 dosage-covariance pairs.

Each candidate was matched on its assumed GRCh37 position, rsID and two-allele contrast to a 1000 Genomes phase 3 reference record. All 503 individuals had called diploid genotypes for all 22. We reproduced the prior count SD and verified that the variance components sum to the observed count variance. [Machine-readable results](shared22_count_variance_diagnostic.json) include all pair contributions and the source-score hash; the calculation script is `work/shared22_count_variance.py`. No individual genotypes or outcomes were saved.

This is a small, pooled reference panel, and the marker selection, source allele conventions and multiallelic representations remain uncertain. Neither these contributions nor the 1 Mb grouping establish the number of independent causal signals. A properly developed LD-aware score and independent validation in each target population are required to evaluate predictive value; the existing count remains an unvalidated research prototype.
