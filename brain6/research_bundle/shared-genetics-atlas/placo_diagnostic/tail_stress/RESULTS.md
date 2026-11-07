# PLACO+ fixed-effect tail stress

For a particular one-trait-only Gaussian stress model, the exploratory PLACO+ rule can reject more frequently than its nominal genome-wide screening threshold. This is stronger evidence that the 82 numerical crossings must stay **exploratory**, but it does not estimate the real data's false-discovery rate or show that the method is generally invalid under the authors' intended model.

We fixed the nuisance parameters estimated for the FI and FFS scans, used unit-variance Gaussian Z scores with their fitted correlation, and set the true mean of one Z to μ while setting the other to zero. The event includes the scan's exclusion of either |Z| > √80. The official PLACO+ rule crosses the two-pair numerical cutoff α=2.1963×10⁻⁸ at absolute Z products of approximately 20.135 for FI and 20.206 for FFS. The product distribution is integrated directly; no rare-event Monte Carlo extrapolation is used.

| True non-null Z mean μ | FI probability of crossing | FFS probability of crossing | FI / α | FFS / α |
|---:|---:|---:|---:|---:|
| 0 | 7.80×10⁻¹⁰ | 4.61×10⁻¹⁰ | 0.036 | 0.021 |
| 1 | 2.07×10⁻⁸ | 1.40×10⁻⁸ | 0.94 | 0.64 |
| 2 | 5.79×10⁻⁷ | 4.41×10⁻⁷ | 26.4 | 20.1 |
| 3 | 9.12×10⁻⁶ | 7.63×10⁻⁶ | 415 | 348 |
| 4 | 8.59×10⁻⁵ | 7.68×10⁻⁵ | 3,911 | 3,499 |

These are **per-marker probabilities under a stipulated composite-null scenario**—one phenotype can have a real fixed association while the other has none. They are not family-wise error rates and do not identify which real markers, if any, satisfy this scenario. The fitted correlation is used as the hypothetical sampling-error correlation, although its interpretation in these heterogeneous GWAS requires further study. Genome-wide LD, sample overlap, uncertainty in the nuisance estimates, actual effect-size distributions, source precision, phenotype construction and selection are not simulated. The results do not establish a bound over all null parameters or an adjusted p-value for any candidate.

For numerical verification we independently conditioned on each Gaussian coordinate in turn, repeated the integral with tighter tolerance, and checked sign symmetry. At nominal α=0.05 and μ=2, the quadrature yields rejection rates 0.174 for FI and 0.167 for FFS, consistent with the earlier 1,000-draw simulation within a prespecified five-standard-error tolerance. Code, all μ=0–6 results, and source hashes are in this directory. See the [simulation review](../NULL_SIMULATION_REVIEW.md), [official method paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12547305/), and [completed candidate audit](../completed_audit/RESULTS.md).

A [cap-bound diagnostic](cap_bound/RESULTS.md) shows why a minimal worst-case guard based only on the existing Z cap cannot recover these candidates: it is always at least as conservative as max-p.

A subsequent [model-consistent random-effect null calculation](model_null/RESULTS.md) finds rejection probabilities below the same threshold under PLACO+'s stipulated zero-mean inflated-variance distributions. This distinguishes that assumption from the fixed-mean stress cases; neither calculation verifies the real GWAS null model.
