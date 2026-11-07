# Sleep-only and separate-question sensitivity

**Historical coding:** The table below predates the [HUQ051 visit-code correction](../visit_code_audit/RESULTS.md). Corrected full FI–sleep and FI-minus-sleep-and-fatigue correlations are 0.31078 and 0.23804, respectively; the complete rerun is linked there.

On the same **2,433 people** in the approximate NHANES frailty pilot, removing the sleep item alone reduced its correlation with that very item; removing the low-energy item as well reduced it further. A distinct question about ever telling a health professional about trouble sleeping showed a smaller numerical change. These are descriptive, survey-weighted within-person associations using working scores, not causal effects or genetic correlations.

| FI version | Correlation with DPQ030 sleep symptom | Correlation with SLQ050 disclosure of sleep trouble |
| --- | ---: | ---: |
| Full working FI | 0.311 | 0.250 |
| Without sleep item | 0.256 | 0.238 |
| Without sleep and low-energy items | 0.238 | 0.232 |

The paired full-minus-sleep-removed difference is **0.054** for DPQ030 (95% survey interval 0.050–0.058) and **0.013** for SLQ050 (0.010–0.016). Removing low energy after sleep produces additional differences of **0.019** (0.014–0.023) and **0.005** (0.002–0.009), respectively. The original, direct question is part of the full FI, so its higher correlation is partly a mathematical consequence of shared measurement. The estimates also reflect rescaling, correlations among deficits, selection and the particular response coding; they do **not** quantify causal effects of sleep or fatigue or isolate all measurement overlap. These post-hoc intervals are descriptive and not corrected for the broader project.

The two sleep questions have different meanings: [DPQ030](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DPQ_G.htm) asks about trouble sleeping **or sleeping too much** in the last two weeks, whereas [SLQ050](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/SLQ_G.htm) asks if the participant **ever told a health professional** about trouble sleeping. Neither establishes an insomnia diagnosis. SLQ050 is not a component of the FI but is collected from the same people and may depend on healthcare access; the FI also includes healthcare-use items. It is not an independent replication. The FI remains a broad deficit index after exclusions, and the original score construction and missing-data rules are not independently validated.

Both questions were available for all participants in the selected primary domain. Results use the full mobile-exam survey design, combined-cycle weights and 32 design degrees of freedom. The [protocol](PROTOCOL.md), [counts](counts.tsv), [all correlations](correlations.tsv), [paired contrasts](contrasts.tsv), [code](followup.R) and [run log](verified_run.log) retain the analysis. All 32 source hashes were rechecked; direct weighted formulas reproduced correlation estimates, finite-difference gradients checked survey contrast standard errors, and all six reported differences matched the corresponding saved correlations. These are arithmetic and file checks, not validation of the scoring model. The [parent pilot](../RESULTS.md) documents sampling, score construction, the mobility skip issue and the distinction from UKB genetic FI–FFS correlation. No person-level derived data were exported.
