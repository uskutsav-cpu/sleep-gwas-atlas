# Effect-difference and heterogeneity validation

Audit date: 2026-10-08. No historical heterogeneity flag is deleted.

The historical calculation is Δ=rg_replication−rg_discovery, SE(Δ)=√(SE_discovery²+SE_replication²), Z=Δ/SE(Δ), Q=Z² and P=erfc(|Z|/√2). Arithmetic is correct for all 41 estimated comparisons **conditional on zero covariance between the two rg estimators**. Since the sleep GWAS is reused, this independence assumption is unverified. The proper variance is SE_discovery²+SE_replication²−2Cov(rg_discovery,rg_replication). A cross-trait LDSC intercept is not the covariance of two rg estimators. Joint genomic block jackknife or another valid covariance estimator is needed for calibrated difference inference. Standard LDSC deals with within-pair overlap under its model; it does not establish independence of repeated estimates ([Bulik-Sullivan et al., 2015](https://www.nature.com/articles/ng.3406)). Difference inference with correlated summary estimates can require a joint jackknife ([Howe et al., 2022](https://www.nature.com/articles/s41588-022-01062-7)). These methodological sources were checked on 2026-10-08.

The original seven flags are nominal P<0.05 in the selected 23-success subset. They are neither seven family-wise significant findings nor calibrated evidence of different biology. Among all 41 tested comparisons,16 meet nominal P<0.05;10 meet exploratory BH across41;4 meet exploratory Bonferroni across41. These new multiplicity summaries are retrospective diagnostics and do not alter the discovery or replication correction families.

Every one of the seven flagged pairs remains direction-concordant and Bonferroni-significant on the replication side. In all seven the replication estimate is smaller than the discovery estimate. Source differences, ascertainment, selection of discovery extremes, winner's curse and unknown estimator covariance are possible explanations; the evidence does not isolate a causal explanation.

The `tables/HETEROGENEITY_MASTER.tsv` retains all 41 estimated comparisons, including the34 outside the preserved seven-success flagged subset, plus source definition comparisons, both estimates, both intervals, effect differences, zero-covariance difference intervals, and hypothetical covariance-grid P values. `heterogeneity_7.tsv` is the preserved historical flagged subset. It does not hide the other34 estimated comparisons.

## Executed covariance sensitivity

For each of the41 comparisons, the audit substitutes Cov=ρ·SE_discovery·SE_replication with ρ∈{−0.9,−0.5,0,0.5,0.9}. This is a hypothetical sensitivity grid, not an estimated participant-overlap correction. No ρ is preferred based on the outcomes.

| Assumed estimator ρ | Nominal P<0.05, all 41 | Exploratory Bonferroni, all 41 |
|---:|---:|---:|
| −0.9 | 9 | 1 |
| −0.5 | 10 | 2 |
| 0 | 16 | 4 |
| 0.5 | 22 | 11 |
| 0.9 | 32 | 24 |

The count depends materially on unknown covariance. Thus `DIRECTION_REPLICATED_EFFECT_MAGNITUDE_QUALIFIED` is appropriate for the seven flags; absence of a nominal difference for other pairs does not demonstrate equivalence. No pooling or meta-analysis of discovery and replication estimates was performed. Native joint-covariance inference is `BLOCKED_MISSING_JOINT_JACKKNIFE_INPUTS`.

## Preserved seven comparisons

| Sleep trait | External phenotype | Discovery rg (SE) | Replication rg (SE) | Δ | Nominal P, covariance0 | Match |
|---|---|---:|---:|---:|---:|---|
| insomnia | Abdominal pain | 0.5591 (0.0359) | 0.4098 (0.0332) | -0.1493 | 0.00226361 | EXACT |
| insomnia | K21 Gastro-oesophageal reflux disease | 0.5216 (0.0378) | 0.364 (0.0375) | -0.1576 | 0.00307768 | EXACT |
| insomnia | Diaphragmatic hernia | 0.4325 (0.0369) | 0.1915 (0.039) | -0.2410 | 7.16469e-06 | EXACT |
| insomnia | Noninfectious gastroenteritis | 0.5211 (0.0497) | 0.3062 (0.0602) | -0.2149 | 0.00590812 | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| shortsleep | K21 Gastro-oesophageal reflux disease | 0.3819 (0.0419) | 0.1961 (0.0357) | -0.1858 | 0.000737252 | EXACT |
| shortsleep | Noninfectious gastroenteritis | 0.4238 (0.0503) | 0.2062 (0.0534) | -0.2176 | 0.00301499 | COMPARABLE_WITH_DOCUMENTED_DIFFERENCES |
| shortsleep | J44 Other chronic obstructive pulmonary disease | 0.4258 (0.0549) | 0.2584 (0.0331) | -0.1674 | 0.00902021 | EXACT |

Every row has per-source comparison and ascertainment qualifications in the machine-readable master. Differences in Pan-UKB questionnaire/ICD/phecode definitions and Finnish register or MVP veteran ascertainment are documented characteristics, not demonstrated causes of the effect-size differences. Exact phenotype matching does not eliminate ascertainment differences.
