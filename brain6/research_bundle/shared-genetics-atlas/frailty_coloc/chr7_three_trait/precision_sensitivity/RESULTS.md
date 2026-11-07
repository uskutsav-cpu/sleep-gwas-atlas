# Chromosome 7: FI precision stress test

The FI–insomnia cluster persisted in **all 200** illustrative rounding perturbations, but the candidate marker changed: rs12671330 remained the leader in **119 of 200** runs, while ten other markers led at least once. Thus this specific stress test leaves the model's regional cluster relatively stable but shows that its top-marker ranking is fragile.

| Result | Baseline | Range across 200 perturbations |
|---|---:|---:|
| HyPrColoc FI–insomnia cluster posterior | 0.8822 | 0.8781–0.8878 |
| FI–insomnia ABF H4, p12=10⁻⁵ | 0.9528 | 0.9510–0.9552 |
| FI–insomnia ABF H4, p12=10⁻⁶ | 0.6688 | 0.6602–0.6805 |
| FI–ADHD ABF H4, p12=10⁻⁵ | 0.0441 | 0.0334–0.0513 |
| FI–ADHD ABF H4, p12=10⁻⁶ | 0.00460 | 0.00344–0.00538 |

The remaining leading markers were rs12705982 (34 runs), rs4730637 (11), rs1005958 (7), rs1378771 (6), rs60234044 (6), rs4236599 (5), rs12154391 (4), rs2396766 (4), rs936145 (3), and rs1005959 (1). These variants are correlated and were selected post hoc; these frequencies are not probabilities that they are causal.

Each run independently perturbed every FI beta and SE by a random amount between −0.00005 and +0.00005, while keeping the other two traits and the 1,055-marker support unchanged. The assumed interval is compatible with rounding to four decimal places, but the actual source rounding procedure and correlations among errors are unknown; perturbations were not constrained to reproduce original p-values. Ranges are **not confidence intervals, credible intervals, or estimated error rates**. They do not establish robustness to different cohorts, missing variants, phenotype definitions, overlap, LD or other model assumptions. The default HyPrColoc model still assumes independent study statistics and at most one causal variant per trait. FI includes sleep-related items, and the failed LAVA QC gates and independent-replication gap remain unresolved.

The [protocol](PROTOCOL.md), [200 replicates plus baseline](replicates.tsv), [code](run.R), [run log](run.log) and [verification](verification.json) preserve the seed and inputs. The baseline reproduced the previously stored pairwise posteriors, and row counts and cluster outputs were checked; this is not statistical calibration. The original source-derived input remains on the SSD and is subject to its noncommercial research and redistribution restrictions.
