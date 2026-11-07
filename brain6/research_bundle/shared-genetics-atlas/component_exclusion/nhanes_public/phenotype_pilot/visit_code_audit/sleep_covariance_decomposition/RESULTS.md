# How the included sleep item contributes to working FI–sleep covariance

In the corrected, exploratory NHANES sample, the weighted covariance between the working FI and the same sleep question is **0.009318**. Because the question is a component of the FI, the score can be decomposed exactly, person by person, into its sleep term and all remaining terms:

**FI = sleep/n + (1 − 1/n) × FI-without-sleep**, where n is the number of observed FI items. In the 2,433-person analysis domain, n ranges from 41 to 49, with median 49. On that score and response scale:

| Term in the exact covariance identity | Weighted covariance with sleep | 95% survey interval |
|---|---:|---:|
| Included sleep term, sleep/n | 0.001746 | 0.001548–0.001944 |
| Remaining weighted FI terms | 0.007572 | 0.005838–0.009307 |
| Full working FI | 0.009318 | 0.007481–0.011155 |

The first term is about **18.7% of the observed covariance** in this specific mathematical decomposition. This is a descriptive arithmetic share, not a causal contribution, variance explained, correction for shared measurement, or the change in correlation after deleting an item. The two terms are not assumed independent, and their intervals should not be added. The remainder may still correlate with sleep because of related deficits, behavior, illness, measurement and selection; the varying denominator also matters. The previously reported weighted FI–sleep correlation is **0.3108** and was reproduced, but covariance terms cannot be read as additive fractions of that correlation.

The sleep question includes trouble sleeping **or sleeping too much** in the preceding two weeks and is not an insomnia diagnosis. This is a survey-weighted observational analysis with provisional working scores, not an item-excluded GWAS or a genetic correlation. The primary sample omits participants for whom the quarter-mile question was skipped, including those reporting difficulty walking without equipment; its selection can affect these estimates. The broader FI and adapted Fried score remain unvalidated reconstructions, and the calculation does not isolate physical frailty or establish causality.

The [protocol](PROTOCOL.md), [scoring copy](score_copy.R), [analysis](analyze.R), [covariances](covariances.tsv), [sample summary](summary.tsv), [run log](run.log), and [receipt](receipt.json) document the arithmetic. All 32 original public-file hashes matched prior receipts; person-level score identity and direct weighted covariances were checked in memory, and the previously corrected correlation was reproduced. The only changes to the scoring copy were its output directory and adding the number of observed FI components to memory. No participant-level derived records were saved.
