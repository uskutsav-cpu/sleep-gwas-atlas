# Reference-group sensitivity of the 22-marker prototype

The same fixed 22 FI weights produce different distributions of raw genotype sums in the five labeled 1000 Genomes superpopulation reference groups. The within-group sample standard deviation of the FI-weighted sum ranges from **0.0316 to 0.0447**; the FFS-weighted sum ranges from **0.0144 to 0.0192**. This illustrates that centering, scaling and allele frequencies can depend on the chosen reference sample. It does **not** measure disease risk, predictive performance, clinical validity or relative risk between groups.

| 1000 Genomes reference label | People | FI-sum SD | FFS-sum SD | Insomnia-sum SD |
|---|---:|---:|---:|---:|
| AFR | 661 | 0.0316 | 0.0144 | 0.0591 |
| AMR | 347 | 0.0394 | 0.0161 | 0.0686 |
| EAS | 504 | 0.0320 | 0.0162 | 0.0600 |
| EUR | 503 | 0.0447 | 0.0192 | 0.0778 |
| SAS | 489 | 0.0413 | 0.0186 | 0.0715 |

The source labels are broad, internally diverse reference-panel groupings, not biological categories or a proxy for any individual's identity. These small, potentially related and admixed samples are not representative of all populations; no uncertainty intervals or ancestry-PC adjustment were calculated. Their differences in frequency and LD do not show how the candidate effects transfer to another population. The source GWAS and weights were derived predominantly in European-ancestry datasets, and outcome-based performance must be tested in each intended population rather than inferred from this table.

For example, the C allele frequency for rs6058750 is 0.278 in the EUR-labeled reference sample and 0.935 in the AFR-labeled sample. That is a difference in a specific reference allele count, **not** a disease-risk difference or evidence that this variant has a causal effect. For multiallelic rsIDs, the matched two-allele record does not prove that the original GWAS used the same representation. Labels and frequencies for all candidates are in [the frequency table](22_score_reference_group_allele_frequencies.csv).

## Method and limitations

We matched the same 22 specified biallelic records by rsID, assumed GRCh37 coordinate and alleles in 1000 Genomes phase 3, with complete diploid calls for 2,504 individuals. The three separate scores use unchanged FI, FFS or UKB insomnia marginal weights. EUR score correlations reproduced the earlier 503-person analysis within 1e−12. We did not store individual genotypes or sums. The [receipt](22_score_reference_group_diagnostic.json) records aggregate group counts, score SDs, allele frequencies, cross-score correlations and source-weight hash; script: `work/score22_global_reference_diagnostic.py`.

The weights are selected after looking at related results, are not LD-adjusted, and have unresolved study-specific coding and FI cohort-consistency questions. No FI, FFS, insomnia or disease outcomes were used. The sums are not calibrated probabilities and must not be applied to clinical decisions. A meaningful portability evaluation needs independent outcome data, source/target overlap checks, appropriate genotype and ancestry QC, calibration and predictive performance with uncertainty, as described in [the validation design](independent_outcome_validation_cohort_feasibility.md).
