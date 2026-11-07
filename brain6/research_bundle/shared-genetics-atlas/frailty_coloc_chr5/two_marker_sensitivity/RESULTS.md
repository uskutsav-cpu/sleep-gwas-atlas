# Chromosome 5: two-marker LD residual sensitivity

The existing FI, FFS and registry-insomnia associations are compatible with different **relative patterns** across the region, but this simplified calculation cannot establish distinct or shared causal signals. Using previously selected markers and p-derived signed Z values:

| Target given the other marker | FI residual Z, EUR / FIN | FFS residual Z, EUR / FIN | Registry insomnia residual Z, EUR / FIN |
|---|---:|---:|---:|
| rs7703746 given rs1363103 | 0.80 / 1.01 | 2.91 / 3.10 | 2.25 / 2.30 |
| rs1363103 given rs7703746 | 3.26 / 3.31 | 2.79 / 2.98 | 0.68 / 0.83 |
| rs2431108 given rs1363103 | −1.89 / −2.28 | −2.29 / −2.79 | −2.37 / −2.54 |
| rs1363103 given rs2431108 | 3.63 / 3.85 | 4.70 / 4.96 | 1.51 / 1.77 |

All Z values count rs7703746 G, rs1363103 T or rs2431108 T. The signed LD is positive for the first pair and negative for the second. These are **unadjusted descriptive residuals**, not validated conditional-GWAS p-values, independent loci, or formal tests of differences between phenotypes. They vary with the conditioning marker and reference panel. The 99-person FIN panel is nested within the 503-person EUR panel and is not an independent cohort replication.

The residual `(Z_target − r × Z_condition)/sqrt(1−r²)` relies on a simplified two-marker Gaussian model with matching study LD, association-statistic covariance and sample composition. The saved correlations are estimated from phased reference haplotypes, not the study cohorts' genotype LD; the reference is small and ancestry/LD mismatch, rounded FI beta/SE, FFS model-specific statistics and sample overlap can materially affect the result. FFS A1 orientation is supported by the [37-lead source audit](../../replication/ye_2023/allele_adapter/RESULTS.md) but not directly checked for every regional row. FFS, FI and insomnia differ in scale and ascertainment. No multivariate conditioning on all regional markers or multiple-signal fine-mapping was performed.

The [complete table](results.tsv) also reports an alternative beta/SE-derived Z for all 24 trait–reference–ordered-pair combinations and descriptive normal tails. The difference between source p-derived and Wald Zs is a model/precision sensitivity, not grounds to choose the more favorable result. All markers and comparisons are post-selected, and no multiplicity-corrected discovery is claimed. The FFS coloc gate, prior-sensitive FI–insomnia model and failed LAVA quality gates remain unchanged.

The [protocol](PROTOCOL.md), [analysis](calculate.py), [source-hash receipt](receipt.json), and [independent arithmetic check](verification.log) preserve this calculation. EUR LD values were cross-checked against separately saved reference calculations; the underlying reference is shared, so this is an arithmetic consistency check rather than independent reference validation.
