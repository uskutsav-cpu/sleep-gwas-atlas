# Conditional p-value enrichment: descriptive check

On 1,138,262 common HapMap3 markers, 26,557 (2.33%) have an approximate registry-insomnia p≤0.01. When markers are conditioned on a smaller frailty p-value, this proportion rises:

| Conditioning trait and cutoff | Observed selected markers | Observed insomnia p≤0.01 | Simulated correlated unit-variance reference | Simulated correlated inflated-variance reference |
|---|---:|---:|---:|---:|
| FI p≤0.01 | 35,854 | 1,491 / 35,854 = 4.16% | 1.26% | 2.50% |
| FI p≤0.001 | 9,029 | 554 / 9,029 = 6.14% | 1.23% | 2.95% |
| FFS p≤0.01 | 44,317 | 1,720 / 44,317 = 3.88% | 1.03% | 2.37% |
| FFS p≤0.001 | 12,135 | 684 / 12,135 = 5.64% | 1.20% | 2.50% |

The simulated columns use one million zero-mean Gaussian draws per pair with the fitted correlation parameter, either unit variances or the genome-wide variance estimates used by the PLACO+ diagnostic. They are mathematical reference scenarios, **not estimates of the actual genetic null**. They omit linkage disequilibrium, realistic trait architecture, participant relatedness, uncertainty in the fitted parameters, and phenotype ascertainment. The observed markers are LD-correlated, so the denominators are not independent observations and the table does not provide valid conventional confidence intervals. Some simulation subsets have few events; the complete counts and all three p cutoffs are in `conditional_fractions.tsv`.

These findings motivate further investigation of cross-trait enrichment, but do not establish its biological cause, the exact role of sample overlap or FI sleep items, or the presence of any particular shared causal variant. The FFS input is the aligned source Z vector and does not depend on the corrupted PLACO output. This diagnostic is **not** the published conditional/conjunctional FDR method, and no FDR-controlled SNPs were estimated. The repository's official [pleioFDR implementation](https://github.com/precimed/pleiofdr) requires an appropriate reference and computational setup; its documented full example recommends at least 16 GB RAM, while this machine has 8 GB. Sample-overlap correction and reference LD/pruning would need review before a valid analysis.

The [protocol](PROTOCOL.md), `run.R`, `conditional_fractions.tsv` and `receipt.json` record thresholds, directions, seed, inputs, hashes and verification. The reverse conditioning direction and all thresholds are retained in the table. Validity remains dependent on the upstream GWAS and normal-tail approximations from rounded Z scores; it is not a substitute for independent two-trait replication or calibrated local analysis.

A [physical-block sensitivity](block_sensitivity/RESULTS.md) gives similar descriptive uncertainty estimates using 5 Mb and 10 Mb blocks. These windows do not guarantee independence or validate a genetic-null test.

A [method and input-readiness audit](conjfdr_readiness/RESULTS.md) records the official reference and runtime requirements and quantifies ties and threshold sensitivity from the existing rounded Z inputs; a published-style conjFDR run remains unperformed.
