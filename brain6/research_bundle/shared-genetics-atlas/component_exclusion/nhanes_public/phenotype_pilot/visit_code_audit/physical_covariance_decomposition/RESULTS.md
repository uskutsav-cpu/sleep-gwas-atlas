# Shared physical questions in the working FI–Fried covariance

In the same corrected 2,433-person NHANES pilot, the weighted covariance between the working FI and adapted Fried score is **0.05682**. Two questions—difficulty walking a quarter-mile and difficulty lifting or carrying—are included in both scores. On this sample, an exact per-person FI identity separates the two FI terms from the remaining FI terms:

**FI = (walking + lifting)/n + (n − 2)/n × FI-without-these-two-items**, where n is the number of observed FI deficits. Both shared responses are observed in the analyzed domain.

| Component of the FI in that identity | Weighted covariance with the full adapted Fried score | 95% survey interval |
|---|---:|---:|
| Two shared-question terms | 0.00944 | 0.00820–0.01068 |
| Remaining FI terms | 0.04738 | 0.04048–0.05428 |
| Full FI | 0.05682 | 0.04907–0.06456 |

The two shared-question terms account for approximately **16.6% of this observed covariance under this particular algebraic decomposition**. They are correlated with the full adapted score, which itself contains both questions; the remainder may also correlate with those and other components. This is not a causal contribution, the isolated effect of shared measurement, or the change in correlation caused by removing the questions. Covariance components are not additive components of the correlation, their uncertainty intervals are not independent, and their intervals should not be added. The existing paired score-deletion sensitivity changes the score content and denominators and answers a different question.

The original weighted FI–adapted-Fried correlation of **0.5642** was reproduced. These are provisional, survey-weighted **phenotypic** results from noninstitutionalized NHANES participants, not the UKB-based **genetic** FI–FFS correlation. The adapted Fried score uses self-reported functional proxies, and the complete-score domain excludes participants whose quarter-mile response was skipped, including those reporting difficulty walking without equipment. Score construction, selection and common influences remain limitations. No sleep-item-excluded GWAS, causal effect, independent replication or validated physical-only construct is implied.

The [protocol](PROTOCOL.md), [scoring copy](score_copy.R), [analysis](analyze.R), [covariance table](covariances.tsv), [sample summary](summary.tsv), [run log](run.log) and [receipt](receipt.json) preserve the analysis. All 32 public input files matched prior hashes. The score identity and weighted covariances were checked per person and by direct arithmetic in memory, and the previously corrected correlation was reproduced. No participant-level derived data were exported.
