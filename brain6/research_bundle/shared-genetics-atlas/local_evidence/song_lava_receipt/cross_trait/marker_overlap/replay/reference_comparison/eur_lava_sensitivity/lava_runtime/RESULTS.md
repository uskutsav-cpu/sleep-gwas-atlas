# Direct LAVA univariate check using public-reference LD

We reconstructed the selected 2,738-marker 1000 Genomes European dosage LD matrix in memory and supplied it to the pinned LAVA runtime for locus 464. We deliberately retained the original UKB reference allele frequencies for LAVA's binary insomnia conversion, so this is a **hybrid reference sensitivity**, not a fully independent 1000 Genomes analysis. The FI variance reproduced the earlier independent SVD-based calculation to within 10⁻¹², checking the numerical implementation on the same support.

| Phenotype | LAVA variance diagonal | Univariate p, unadjusted |
|---|---:|---:|
| FI | +1.881×10⁻⁶ | 0.001905 |
| Insomnia | +2.312×10⁻⁶ | 3.87×10⁻⁹ |

Neither phenotype triggered LAVA's negative-variance drop in this one calculation. These are **post-selected, uncorrected univariate values**, not evidence that a frozen univariate or bivariate testing family passes, nor a local genetic correlation. No bivariate test was run. The 503-person reference is small, is not established as matched to either GWAS, and has shown material sensitivity to resampling and population composition; the original local family quality gates remain failed. A positive result under this reference cannot establish that the original negative FI variance was erroneous or that either result is calibrated.

## Reproducibility

The [pre-run protocol](PROTOCOL.md), [Python reconstruction](run.py), [LAVA runtime script](run.R), [runtime output](lava_result.json), and [receipt](receipt.json) record the analysis. The dosage LD matrix passed in memory directly to R; neither it nor individual genotypes were saved. The R script substitutes a runtime-only matrix reader and restores the original afterward; it does not modify the source reference. The separate [verification](verify.R) independently calculates the reported F and chi-square tails from saved projected sums and checks the FI result against the prior SVD analysis. These checks do not establish reference validity, sample independence, sample-overlap correction or false-positive calibration. The source VCF/index and the full historical inputs remain subject to the [parent limitations](../RESULTS.md).
