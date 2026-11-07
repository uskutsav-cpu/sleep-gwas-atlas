# FI–insomnia regional coloc: rounding stress test

The FI beta and standard-error fields used in the saved matched files have at most four decimal places. To assess how this precision could affect the three post-selected single-signal colocalizations, we independently perturbed each FI beta and SE by a uniform random amount between −0.00005 and +0.00005, recalculated the standard coloc approximate Bayes factors and repeated 500 times per region using a fixed seed. FinnGen inputs, variant sets, effect priors and association priors were unchanged. The assumed perturbation range is consistent with rounding to four decimal places, but the original rounding procedure has not been established; independent uniform errors are an illustrative stress test, not the posterior distribution of the underlying statistics. Perturbations were not constrained to reproduce the source p-values. The resulting ranges are neither confidence nor credible intervals.

| Region, GRCh37 | Baseline PP4, p12=10⁻⁵ | Perturbed PP4 range | Baseline PP4, p12=10⁻⁶ | Perturbed PP4 range | Conditional lead across 500 perturbations |
|---|---:|---:|---:|---:|---|
| chr3 | 0.00199 | 0.00195–0.00204 | 0.000199 | 0.000195–0.000204 | rs1499984 500 times |
| chr5 | 0.512 | 0.470–0.557 | 0.0951 | 0.0815–0.112 | rs1363103 331, rs13164188 169 |
| chr7 | 0.587 | 0.538–0.653 | 0.125 | 0.104–0.158 | rs10953765 269, rs1456031 231 |

The chr5 and chr7 top conditional variant can change under these small input perturbations, reinforcing that the ranks are provisional. Chr3's stable lead is conditional on H4 despite very low PP4 and does not establish a shared causal variant. Prior sensitivity remains substantial in all three analyses; these perturbations do not resolve single-causal-variant assumptions, unmeasured variants, phenotype overlap, ancestry differences, or independence. Full-precision source data are needed for a definitive precision check.

`run.R` records the seed, input files and model, and `replicates.tsv` contains the baseline and all 500 realizations for each pair of priors. The ordinary regional results are in [chr3 and chr7](../RESULTS.md) and [chr5](../../frailty_coloc_chr5/RESULTS.md).
