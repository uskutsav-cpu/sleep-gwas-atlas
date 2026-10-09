# Phase 0/1 report: Sleep/Circadian Genetic Atlas

**Input status:** Derived from the listed local LDSC output tables. This report does not independently validate source-GWAS provenance or reproduce LDSC.

## 1. Trait registry

- Registered traits: **45**
- Traits by domain:
  - `aging`: 6
  - `cancer`: 6
  - `cardio`: 4
  - `immune`: 6
  - `metabolic`: 5
  - `neuro`: 2
  - `psychiatric`: 4
  - `sleep`: 12
- Source-verification status by declared value:
  - `SOURCE_PENDING`: 3
  - `SOURCE_VERIFIED`: 42

## 2. SNP heritability and QC gate

- Traits with parsed LDSC h2: **27**
- PASS: **27**; DROP: **0**

Highest h2 Z-scores (descriptive only):

| trait | scale | h2 | se | z | intercept | verdict | qc_reason | ldsc_n_eff_h2 | ldsc_n_eff_h2_gt_12000 | mixer_univariate_status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| bmi | observed | 0.1908 | 0.0053 | 36.0 | 1.187 | PASS | pass | 130000.0 | True | NOT_RUN |
| scz | liability | 0.2282 | 0.0074 | 30.84 | 1.079 | PASS | pass | 28820.0 | True | NOT_RUN |
| napping | observed | 0.0812 | 0.0029 | 28.0 | 1.05 | PASS | pass | 36750.0 | True | NOT_RUN |
| sbp | observed | 0.1325 | 0.0048 | 27.6 | 1.153 | PASS | pass | 100400.0 | True | NOT_RUN |
| snoring | liability | 0.1054 | 0.0041 | 25.71 | 1.031 | PASS | pass | 40260.0 | True | NOT_RUN |
| shortsleep | observed | 0.0647 | 0.0027 | 23.96 | 1.027 | PASS | pass | 20400.0 | True | NOT_RUN |
| sleepdur | observed | 0.0664 | 0.0029 | 22.9 | 1.038 | PASS | pass | 29620.0 | True | NOT_RUN |
| adhd | liability | 0.2405 | 0.0107 | 22.48 | 1.042 | PASS | pass | 30840.0 | True | NOT_RUN |
| insomnia | liability | 0.1008 | 0.0045 | 22.4 | 1.015 | PASS | pass | 31630.0 | True | NOT_RUN |
| frailty | observed | 0.1093 | 0.005 | 21.86 | 1.02 | PASS | pass | 19150.0 | True | NOT_RUN |

Gate definition: h2 Z >= 4 and LDSC intercept <= 1.20. A low h2 Z and an inflated intercept are reported as distinct failure modes.

## 3. Genome-wide genetic correlation

- Parsed sleep-disease pairs: **136**
- FDR < 0.05: **51**

Strongest positive estimates (descriptive only):

| sleep_trait | disease_trait | rg | se | p | fdr |
| --- | --- | --- | --- | --- | --- |
| insomnia | frailty | 0.6405 | 0.0235 | 4.322e-163 | 5.877e-161 |
| insomnia | adhd | 0.4033 | 0.0318 | 5.836e-37 | 2.6460000000000003e-35 |
| snoring | bmi | 0.3779 | 0.0187 | 1.029e-90 | 6.998999999999999e-89 |
| insomnia | healthspan | 0.3432 | 0.042 | 3.178e-16 | 5.402e-15 |
| snoring | healthspan | 0.2927 | 0.0368 | 1.828e-15 | 2.26e-14 |

Strongest negative estimates (descriptive only):

| sleep_trait | disease_trait | rg | se | p | fdr |
| --- | --- | --- | --- | --- | --- |
| insomnia | parental_lifespan | -0.274 | 0.0367 | 8.869e-14 | 8.615e-13 |
| sleep_efficiency | ovarian_cancer | -0.2289 | 0.0712 | 0.0013 | 0.004533 |
| sleepdur | frailty | -0.222 | 0.0275 | 6.817e-16 | 1.03e-14 |
| napping | parental_lifespan | -0.1793 | 0.0274 | 5.969e-11 | 4.273e-10 |
| snoring | parental_lifespan | -0.1647 | 0.0282 | 4.912e-09 | 3.181e-08 |

## 4. Downstream handoff

1. Use only FDR-significant, well-powered pairs with source provenance retained for Phase 2 local genetic correlation.
2. Confirm binary-trait prevalence citations before reporting liability-scale h2. Treat LDSC N_eff x h2 only as a screening statistic; actual MiXeR eligibility requires univariate MiXeR diagnostics.
3. Preserve each QC ledger, munging log, LDSC log, and table hash with the figure so Phase 1 results remain auditable.
