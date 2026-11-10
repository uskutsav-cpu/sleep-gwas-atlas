# What drives the 22-variant experimental scores?

The three 22-variant weighted sums are affected differently by their constituent variants. In 503 European-labeled 1000 Genomes reference individuals, rs28528997 accounts for about **28.0% of the sample variance of the FFS-weighted sum**, using a covariance allocation described below. The three largest allocations together account for about **50.5%**. By comparison, the top three allocations sum to 30.9% for the FI-weighted sum and 26.3% for the insomnia-weighted sum. These percentages describe variation **in the constructed sums**, not variance of frailty, insomnia, disease risk, heritability, or predictive performance.

| Sum and largest allocations | Allocated fraction of score variance | Standalone variance term / score variance |
|---|---:|---:|
| FI: rs12671330 | 12.1% | 9.6% |
| FI: rs12536335 | 10.4% | 6.6% |
| FI: rs2155292 | 8.5% | 9.1% |
| FFS: rs28528997 | 28.0% | 26.1% |
| FFS: rs11636381 | 11.9% | 12.9% |
| FFS: rs2815764 | 10.5% | 8.4% |
| Insomnia: rs1909650 | 9.9% | 7.6% |
| Insomnia: rs12536335 | 8.8% | 5.8% |
| Insomnia: rs817509 | 7.7% | 6.5% |

The rs28528997 FFS coefficient is particularly influential in this reference sample, but this does not mean that the variant has a causal or clinically important effect. Its rsID is multiallelic in the broader reference and the FI cohort comparison found opposing UK Biobank and TwinGene point estimates after an inferred TwinGene sign adapter. That discrepancy concerns the **FI** coefficient; it does not by itself invalidate the separate FFS coefficient. Both the original coding and the FFS allele convention need independent confirmation. See [cohort comparison](six_leads_fi_cohort_comparison.md) and [candidate review](expanded_22_candidate_risk_gene_review.md).

## Method and interpretation

For each separate sum, let X be the matrix of counted-allele dosages and w the corresponding trait's marginal GWAS weights. With C the sample covariance matrix of dosages, the variance is wᵀCw. A variant's allocated fraction is wᵢ(Cw)ᵢ / (wᵀCw): this assigns each covariance term equally to its two variants and fractions add to 100%. The standalone term is wᵢ²Cᵢᵢ / (wᵀCw); its sum need not be 100%. Allocated terms can be negative and are **not** causal effects or independent contributions. Removing or replacing a variant generally changes both variance and covariance, so an allocated fraction is not necessarily the effect of removing it.

The pooled off-diagonal covariance terms sum to 13.4%, 12.7%, and 13.5% of score variance for FI, FFS and insomnia respectively. These reflect correlations in this finite pooled reference sample, including possible LD, ancestry structure, relatedness and sampling variation; they are not a clean estimate of LD alone. The panel contains several European subgroups, and the calculations are descriptive with no uncertainty interval. Sample frequencies and covariances may differ in other cohorts or ancestries. Use these findings to plan sensitivity analyses, not to pick winners by the largest allocation.

All 22 exact two-allele contrasts were matched in the panel and had complete calls; source-specific multiallelic representation and effect orientation remain unresolved. Weights were selected after prior results and are marginal, not jointly optimized for LD. There are no participant outcomes in this analysis. A real evaluation requires frozen weights, appropriate genotype QC and independent outcome data. The three sums use different scales and should not be added together. The reference-allele frequencies and all variant-level decomposition terms are in [the CSV](reference_22_score_variance_decomposition.csv); calculation metadata and subgroup diagnostics are in [the JSON](reference_22_score_variance_decomposition.json). No individual genotypes or scores were retained. Script: `work/score22_variance_decomposition.py`.
