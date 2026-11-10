# How concentrated are the 15-marker research scores?

In the European 1000 Genomes reference panel, each of the three prototype sums is influenced by a handful of the 15 selected variants. The largest allocated contribution to variation in the FFS sum is **rs323509 (24.7%)**; for the insomnia sum it is **rs113851554 (20.0%)**. These are shares of the variation **of the constructed sums in this reference panel**, not shares of frailty or insomnia risk, variance explained in a phenotype, or causal effects.

| Research sum | Variants with largest allocated shares of score variance | Their allocated shares |
|---|---|---|
| FI | rs323509; rs1909650; rs12536335 | 12.7%; 11.8%; 11.7% |
| FFS | rs323509; rs10865954; rs11790060 | 24.7%; 18.9%; 13.1% |
| Insomnia | rs113851554; rs11152363; rs78206187 | 20.0%; 8.4%; 8.2% |

For this descriptive calculation, each person's sum is the existing trait-specific coefficient multiplied by the counted-allele dosage, summed across 15 variants. We used 503 European-labeled 1000 Genomes Phase 3 samples with complete biallelic calls at the selected contrasts. The fractions are computed using the sample dosage covariance matrix: for score S and marker contribution Xᵢ = βᵢGᵢ, the allocated share is Cov(Xᵢ,S)/Var(S). These shares sum to 100% because each pairwise covariance is split equally between its two markers. This allocation is a bookkeeping convention, not a unique causal or predictive contribution; shares can in general be negative or greater than 100%. It incorporates linkage and any other observed sample covariance, including population structure and finite-sample noise. The score weights are post-selected, the FFS effect convention is inferred, and nine rsIDs have multiallelic records in broader reference data.

The sum of individual markers' standalone variances, ignoring covariance, accounts for 98.2% of the pooled FI score variance, 87.0% of FFS, and 87.5% of insomnia; the remaining net portions (1.8%, 13.0%, 12.5%) are cross-marker covariance in this sample. These are net sums of positive and negative covariances, not the total amount of linkage or evidence of independent loci. Between-chromosome sample covariance and population structure can contribute. The trait-specific sums use different units and should not be added together.

The largest FFS share is consistent with the earlier leave-one-out analysis showing that a few markers can materially change the relationship between prototype sums. For rs323509, its allocated FFS share ranged from approximately 19.8% to 29.5% across five EUR-labeled subgroups; for rs113851554, its insomnia share ranged from 18.3% to 23.1%. These are small subgroup descriptions, not ancestry-transferability comparisons or confidence intervals. Within-subgroup sample sizes are only 91–107, participants may be related, and no outcome information was available. No clinical prediction, calibration, prospective risk, causal mechanism or score performance was evaluated.

## Reproducibility

All 15 exact rsID, GRCh37 position and listed biallelic contrasts matched the same publicly indexed 1000 Genomes Phase 3 reference VCF as the earlier analysis; all 503 genotypes were called. The three recomputed score correlations matched `reference_15_score_diagnostic.json` to 10⁻¹², and the decompositions were checked to sum to total sample variance. The full per-variant values are in `reference_15_score_variance_decomposition.csv`; subgroup results and input/output hashes are in `reference_15_score_variance_decomposition.json`. No person-level genotypes or scores were saved. These diagnostics do not resolve source allele conventions or validate the scoring weights.
