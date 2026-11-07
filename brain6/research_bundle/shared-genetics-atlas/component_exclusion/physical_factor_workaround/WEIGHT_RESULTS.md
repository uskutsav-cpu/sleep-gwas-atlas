# Definition sensitivity of the explicit genetic index

The equal-weight index remains the reported index. A post-hoc calculation across all nonnegative weights, including zero weights and individual indicators, shows how much its point estimate depends on the chosen mixture of physical traits.

| Insomnia source | Possible point correlations across nonnegative weights | Equal-weight result |
|---|---:|---:|
| Registry | −0.0975 to +0.0051 | −0.0736 |
| Questionnaire | −0.2910 to −0.0037 | −0.2309 |

These are mathematical ranges conditional on rounded observed correlations, not confidence intervals, validated weight choices or uncertainty bounds. The extremes use the insomnia outcomes to optimize weights and must not be selected as a new preferred index. The upper endpoint for both sources is activity alone; the most-negative endpoint excludes activity and optimizes the walking/grip mixture. The questionnaire point estimate remains negative across these weights, but uniform statistical evidence does not follow: the activity-only estimate is close to zero with a wide interval.

A fixed path gives activity weight t and walking/grip weights (1−t)/2:

| Activity weight | Registry correlation | Questionnaire correlation |
|---|---:|---:|
| 0 | −0.0950 | −0.2854 |
| 0.25 | −0.0820 | −0.2533 |
| 1/3: equal weights | −0.0736 | −0.2309 |
| 0.50 | −0.0522 | −0.1709 |
| 0.75 | −0.0190 | −0.0750 |
| 1: activity alone | +0.0051 | −0.0037 |

Increasing activity's weight moves both associations toward zero in this particular path. This does not show that activity is unimportant or unrelated to sleep: the activity phenotype, sampling error and model/source differences remain. It shows that the phrase physical function does not uniquely define one genetic construct. The [protocol](WEIGHT_PROTOCOL.md), [calculation](weight_sensitivity.py) and [numerical outputs](joint_fit/weight_sensitivity.json) retain all results. Support enumeration was independently checked against a 0.01-step simplex grid; the grid extrema differed by less than 0.000005. No inferential p values were assigned to optimized weights.

The [scientific figure](explicit_index.png) compares the three individual measures with the fixed equal-weight index for both insomnia sources. A [vector version](explicit_index.svg) and [R plotting source](draw_index.R) are available. Displayed intervals are unadjusted normal 95% intervals from 200-block jackknife SEs. The [index results](INDEX_RESULTS.md) record the limited two-source correction, source assumptions and interpretation; this figure and weight sensitivity supply neither a matched exhaustion-deleted frailty score nor independent validation.
