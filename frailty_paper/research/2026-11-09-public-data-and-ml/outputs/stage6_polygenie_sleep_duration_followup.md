# Published frailty and insomnia scores: exploratory sleep-duration follow-up

The public [PolyGenie GCAT interface](https://polygenie.igtp.cat/) contains three continuous questionnaire outcomes labelled sleep hours. The saved questionnaire tables from the earlier score analysis contain the following displayed estimates for the existing genome-wide Fried-frailty and insomnia polygenic scores, using the top-versus-bottom decile setting. These are **different scores** from this project's two candidate scores. The analysis is exploratory and selected after inspecting the interface.

| Published score | Questionnaire label | Displayed beta | Nominal p |
|---|---|---:|---:|
| Fried frailty | Sleep Hours | −0.0059 | 0.923 |
| Fried frailty | Sleep Hours: Working Days | −0.0056 | 0.933 |
| Fried frailty | Sleep Hours: Non-working Days | +0.0072 | 0.923 |
| Insomnia | Sleep Hours | −0.0978 | 0.0916 |
| Insomnia | Sleep Hours: Working Days | −0.0718 | 0.271 |
| Insomnia | Sleep Hours: Non-working Days | −0.1208 | 0.0997 |

None of these six nominal p-values is below 0.05. This does **not** establish no association, equivalent scores, or replication failure: sample sizes per phenotype and decile, missingness, confidence intervals and exact questionnaire wording and transformation were not available in the saved table and have not been independently verified. Although the labels say hours, do not interpret the numerical beta as an hours difference without verifying the model and units. Sleep duration is not insomnia diagnosis, and these outcomes do not measure frailty. Different signs or p-values are not a statistical comparison between scores; the three duration measures are related. The six results were chosen after inspection amid a much larger PheWAS, so they should not be treated as prespecified confirmatory tests.

The cohort is a distinct named target but participant and relative non-overlap with all discovery cohorts remains unverified. The scores use discovery GWAS also relevant to this project, and the resulting score–phenotype associations cannot validate the 15- or 22-marker prototypes or establish prediction, calibration or causality. The [main PolyGenie review](stage6_polygenie_sleep_frailty_review.md) describes the source scores, cohort, outcome settings and other limitations.

Reproducibility: select rows whose `Code` is `cont_sleephours`, `cont_laborablesleephours`, or `cont_nolaborablesleephours` in the `response.prs-table-content.data` arrays of the saved questionnaire responses. Source responses, previously retrieved 29 October 2026 UTC, are `work/polygenie_Fried_Frailty_Questionnaire.json` (SHA-256 `cc498f14ee04ac9f5e547a2985caae5943213b8f9d9e2821892d9a11e56b4266`) and `work/polygenie_Insomnia_Questionnaire.json` (SHA-256 `9abcba56abe4fde2aaba5dc8cd2ea9e4a8a4d85242c913d1aa8ad003c937dcaf`). The response retains the displayed rounded values, not individual records.
