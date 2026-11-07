# Age and sex sensitivity of FI–Fried phenotype correlations

**Historical coding:** The table below predates the [HUQ051 visit-code correction](../visit_code_audit/RESULTS.md). The corrected age/sex-adjusted full-score estimate is 0.540; its detailed rerun is linked there.

In the same 2,433-person NHANES complete-score pilot, the weighted correlation between the approximate FI and adapted Fried score was **0.563**. The weighted partial correlation after linear adjustment for age and recorded sex was **0.539** (95% survey interval 0.487–0.591). The paired adjusted-minus-unadjusted difference was −0.024 (95% interval −0.035 to −0.014).

| Same-person comparison | Unadjusted r (95% interval) | Age/sex-adjusted partial r (95% interval) |
| --- | ---: | ---: |
| Full approximate FI and adapted Fried score | 0.563 (0.513–0.613) | 0.539 (0.487–0.591) |
| Reduced FI43 and two-component FFS2 sensitivity | 0.275 (0.225–0.325) | 0.248 (0.197–0.299) |

FI43 removes sleep, fatigue, two shared physical questions, BMI and grip from the working FI; FFS2 retains only room-walking difficulty and low activity. The latter is **not a validated Fried score**. Their adjusted-minus-unadjusted difference is −0.027 (95% interval −0.036 to −0.018). The score changes also alter content, denominators and reliability; the comparison cannot attribute differences to any one component.

The partial correlations describe associations between survey-weighted linear residuals after adjustment for top-coded age and the dataset's binary recorded sex variable. They do not establish causation, mediation or adequate control for health, socioeconomic or other confounding. This is a post-hoc sensitivity; intervals are not corrected for the wider research project. Complete-score inclusion excludes people reporting difficulty walking without special equipment because the required question is skipped, and those excluded differed on observed characteristics ([selection audit](../selection_audit/RESULTS.md)). Thus these results do not describe everyone aged 60 or older. Both scores are provisional reconstructions, and neither this observational analysis nor the item deletions provide an item-excluded GWAS or an independent genetic replication.

The [pre-result protocol](PROTOCOL.md), [code](analyze.R), [all estimates](results.tsv), and [run log](run.log) preserve the analysis. The full survey design uses four-year examination weights, cycle-specific strata and nested PSUs; intervals use 32 design degrees of freedom and linearized moment covariance. The weighted residual correlations were independently recalculated using weighted least squares, two finite-difference step sizes agreed, and all 32 source files matched their saved hashes. These checks do not independently validate the score definitions. No participant-level derived records were saved.
