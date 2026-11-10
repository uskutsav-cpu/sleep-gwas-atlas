# If only directly listed Metabochip markers were available

This is a genotype-only feasibility sensitivity, **not an observed BRHS result**. The [manifest screen](brhs_metabochip_score_marker_screen.md) identified one of the 15 markers and two of the 22 on the manufacturer array, with unverified chip revision and alleles. Imputed and other genotype releases may contain additional markers, so this hypothetical subset must not be treated as actual UCLEB data coverage.

Using the same 503 European-labelled 1000 Genomes Phase 3 reference participants as the prior prototype comparison, we compared each full score with the partial sum containing only those directly listed markers. Pearson correlations were:

| Prototype | Partial vs full correlation |
|---|---:|
| 15-marker FI-weighted | 0.234 |
| 15-marker FFS-weighted | 0.283 |
| 15-marker insomnia-weighted | 0.178 |
| 22-marker shared-direction allele count | 0.222 |
| 22-marker FI-weighted | 0.220 |
| 22-marker FFS-weighted | 0.162 |
| 22-marker insomnia-weighted | 0.250 |

The partial sums track their corresponding full sums weakly in this reference sample, so calculating a one- or two-marker sum and calling it the original score would be misleading. These correlations depend on allele frequency, linkage disequilibrium and score weights in the reference population. They are not measures of accuracy in BRHS, predictive performance, fraction of disease risk, or the validity of imputation. The FFS direction remains inferred; the score definitions were selected after prior exploratory work. 1000 Genomes samples include several European subpopulations and possible relatives. Exact variant and allele pairs were matched in the reference VCF, but target data require their own quality and allele checks.

The [script](../work/brhs_chip_subset_reference.py) retrieves the public indexed reference VCF for the 33 distinct markers, requires complete calls and a unique biallelic identifier/position/allele match, applies each prototype's counted-allele coding, and saves only aggregate results. It independently reproduces the three previously recorded full 15-versus-22 score correlations to numerical precision. [Machine-readable results and input hashes](brhs_chip_subset_reference.json). This reproducibility check does not prove the GWAS effect coding or validate a prediction model.
