# PLACO+ tail behavior under its stipulated random-effect nulls

Under the three zero-mean bivariate-normal sub-null distributions specified by the [PLACO+ paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12547305/), the rejection probabilities calculated for the project's fixed nuisance parameters are below the scan's numerical threshold. This complements the earlier fixed-effect stress test: it does **not** establish that real sleep–frailty GWAS satisfy the stipulated distributions.

At the two-pair screening threshold α=2.1963×10⁻⁸, the following ratios are the calculated **per-marker rejection probability divided by α**. The scan's |Z|≤√80 cap is included.

| Fitted pair | Neither trait associated (H00) | Only frailty associated (H01) | Only sleep associated (H02) |
|---|---:|---:|---:|
| Insomnia–FI | 0.036 | 0.721 | 0.315 |
| Insomnia–FFS | 0.021 | 0.821 | 0.200 |

Here H01 and H02 use zero-mean, inflated-variance random-effect Z distributions and the fitted correlation as specified in the protocol. The corresponding uncapped integrals reconstruct the official package's p₁+p₂−p₀ calculation at both the screening threshold and nominal 0.05. Independent conditioning on either coordinate agrees numerically. This shows the implementation behaves conservatively in these particular stipulated sub-null distributions; it does not prove that each real SNP's conditional null p-value is calibrated or that the fitted parameters and GWAS match the model.

The earlier [fixed-mean stress calculation](../RESULTS.md) held one non-null Z mean at specific values and found higher rejection under some such scenarios. Those are different data-generating assumptions, not contradictory numerical results. Actual effect-size distributions, overlapping samples, LD, nuisance-estimation uncertainty and phenotype construction remain unresolved. Neither calculation estimates the real-data FDR or validates the 82 numerical crossings. The conservative max-p screen, independent-data requirements and failed LAVA gates remain unchanged.

An [observed marginal-tail description](../empirical_margins/RESULTS.md) finds more large Z values than a single fitted normal reference. Because the observed distribution contains true signals, linked variants and selection, that comparison is not a test of the component null distributions or of PLACO+ error control.

The [protocol](PROTOCOL.md), [all model probabilities](model_probabilities.tsv), [calculation](calculate.R), [runtime](runtime.txt) and [verification](verification.json) document the model, tolerances and checks. No real SNP was reclassified.
