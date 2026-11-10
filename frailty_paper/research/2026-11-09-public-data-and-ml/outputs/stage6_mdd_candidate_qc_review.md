# Depression-symptom GWAS candidate quality checks

The freely available clinical depression symptom meta-analysis provides a secondary screen of the existing score candidates, not an independent population-based insomnia or frailty validation. Its insomnia and fatigue phenotypes are symptoms assessed in depression cohorts, and participant overlap with discovery studies remains unresolved. Source: [Adams et al.](https://doi.org/10.1017/S0033291724001880); [summary statistics](https://doi.org/10.6084/m9.figshare.22745573). Observe the source's research-use and redistribution conditions.

All 33 distinct candidate variants matched on GRCh37 position and unordered alleles in both symptom files. After orienting the reported control frequency to the same allele, the maximum absolute difference versus the 503-person 1000 Genomes European reference was 0.055 for insomnia and 0.060 for fatigue; none exceeded 0.10. This is a useful coding check, not proof of correct strand, target ancestry, or absence of cohort differences. The reference is finite, and controls can differ from it.

Across the candidates, reported INFO ranges were 0.851–1.000 for insomnia and 0.917–1.020 for fatigue; per-marker effective sample sizes were 18,455–18,857 and 7,706–7,832 respectively. One reported INFO exceeds 1; the interpretation of the meta-analysis field has not been established, so these values should not be treated automatically as conventional bounded imputation r-squared. Effective sample sizes and the clinical selection also limit precision and generalizability.

## Frequency sensitivity for the 15-marker prototype

We repeated the exploratory summary-statistic calculation using a case-count-weighted approximation of the reported case/control allele frequency and Hardy–Weinberg genotype variance, keeping 1000 Genomes European within-chromosome correlations and setting cross-chromosome correlations to zero. The corresponding earlier calculation used empirical reference genotype standard deviations. Meta-analysis frequency summaries need not be exact pooled frequencies, and LD, per-marker sample composition, and model assumptions remain uncertain.

| Weights | Insomnia symptom: original → alternative Z | Fatigue symptom: original → alternative Z |
|---|---:|---:|
| FI | 0.48 → 0.43 | −0.58 → −0.53 |
| FFS | 0.03 → −0.02 | −0.87 → −0.82 |
| Insomnia | 1.09 → 1.03 | −1.05 → −1.00 |

The qualitative interpretation is unchanged: these post hoc approximate Z values do not establish associations or exclude effects. They are not participant-level prediction, a calibrated score test, or validation of frailty. A subsequent [22-marker target-frequency sensitivity](stage6_mdd_22_weighted_sensitivity_review.md) retrieved the necessary reference genotypes and is described separately; it does not solve the target-LD or participant-data limitations.

Reproducible local code: `work/stage6_mdd_candidate_qc.py` and `work/stage6_mdd_15_target_frequency_sensitivity.py`; aggregate outputs: `outputs/stage6_mdd_candidate_qc.json` and `outputs/stage6_mdd_15_target_frequency_sensitivity.json`. Candidate-level source data remain local under the source terms.

## Reference-panel resampling diagnostic

To see how the previous 15-marker approximation reacts to the finite reference sample, we resampled the 503 reference individuals 2,000 times with replacement, recalculating genotype standard deviations and within-chromosome LD while holding every GWAS Z statistic fixed. With insomnia symptom, the 2.5th–97.5th percentiles of the resulting approximate Z values were 0.42–0.54 for FI weights, −0.014–0.073 for FFS weights, and 1.03–1.14 for insomnia weights. With fatigue, they were −0.63 to −0.53, −0.93 to −0.82, and −1.10 to −1.01, respectively. No resample had zero variant variance.

These ranges describe only sensitivity to this particular reference-panel resampling procedure. They are **not confidence intervals for associations or effects**: the GWAS estimates, LD transportability, model assumptions, sample overlap and phenotype differences remain fixed and unaccounted for. The resampling treats individuals as exchangeable, which may be imperfect for a reference panel containing related people. The FFS-weighted insomnia approximation changes sign around zero; its sign should not be interpreted. Code: `work/stage6_mdd_15_reference_bootstrap.py`; aggregate results: `outputs/stage6_mdd_15_reference_bootstrap.json`.
