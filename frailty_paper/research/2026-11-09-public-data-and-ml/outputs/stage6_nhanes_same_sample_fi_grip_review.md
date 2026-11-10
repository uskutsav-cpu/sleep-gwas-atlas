# Same-sample comparison: sleep trouble, deficit index and grip

After the separate associations were known, we compared the same participants using the previously adapted 24-item index and combined grip strength. Both models used the same sleep question, age and sex covariates, and NHANES examination weights, strata and PSUs.

| Cycle | Same participants | Adapted index difference (95% CI), 0–1 units | Grip difference (95% CI), kg |
|---|---:|---:|---:|
| 2011–12 | 1,383 | +0.072 (0.054 to 0.090) | +0.32 (−1.01 to 1.64) |
| 2013–14 | 1,474 | +0.085 (0.062 to 0.108) | −0.84 (−2.49 to 0.81) |

The reported exposure is having ever told a health professional about sleep trouble, not a clinical insomnia diagnosis. The index association is positive in both cycles; the grip estimates include both positive and negative differences and are not precise enough to conclude no association. The comparison supports examining different frailty constructs, but **it is not a comparison of a frailty index with the Fried frailty score**: grip is only one physical component. The index is an unvalidated adaptation that excludes direct sleep and fatigue items but includes symptoms and functional limitations that may relate to sleep. Its coefficients and grip coefficients use different scales, and differing statistical significance is not evidence that their true effects differ. No formal cross-outcome difference test was performed.

These cross-sectional associations may reflect confounding, healthcare access and selection; age/sex adjustment does not establish causation. Complete-case analysis selects on both outcomes and could introduce bias. Earlier results informed this post hoc analysis; it is not an independent discovery, a score-validation study or a test of genetic mechanisms. Participants with no grip measure are excluded, and NHANES survey weights do not remove all missing-data bias.

[Analysis plan](stage6_nhanes_same_sample_fi_grip_plan.md), [results](stage6_nhanes_same_sample_fi_grip_results.csv), [R code](../work/nhanes_same_sample_fi_grip.R), and [index definitions and limitations](stage6_nhanes_adapted_fi_review.md). For a proper FI versus FFS comparison, obtain or construct both measures at the same wave, ensure sleep/exhaustion overlap is addressed, and evaluate genetic and phenotype associations on the same participants with a prespecified analysis.
