# Correction: healthcare visits in the exploratory NHANES pilot

The original working FI inadvertently treated **HUQ051 code 8** as missing. The [2013–14 CDC codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/HUQ_H.htm) defines that code as 16 or more visits. Under the pilot's working threshold, it should be scored as a deficit of 1. The original saved outputs and follow-up reports remain historical results from the earlier coding; use the corrected estimates below for comparisons that were rerun.

| Primary comparison | Original weighted r | Corrected weighted r (95% survey interval) |
| --- | ---: | ---: |
| Full approximate FI vs adapted Fried score | 0.56313 | 0.56416 (0.51405–0.61427) |
| Full approximate FI vs sleep symptom | 0.31057 | 0.31078 (0.26185–0.35971) |
| FI without sleep and fatigue vs sleep symptom | 0.23758 | 0.23804 (0.18822–0.28787) |
| FI without sleep and fatigue vs adapted score without quarter-mile proxy | 0.43255 | 0.43424 (0.37858–0.48991) |

Among positive-examination-weight adults aged at least 60, 115 had code 8; **60** were in the primary 2,433-person complete-score sample, and their FI values changed. Eligibility was identical for all 19,151 design records in the old and corrected primary runs. The paired corrected-minus-original FI–Fried correlation difference was 0.00103 (95% survey interval −0.00038 to 0.00244). The corrected paired difference between the two FI–sleep correlations was 0.07274 (95% interval 0.06609–0.07939). These are phenotype sensitivities, not genetic correlations, clinical validation or proof of no meaningful bias from other assumptions.

The valid HUQ051 code 3 represents **4–5 visits**, spanning the working five-visit threshold. The existing primary implementation assigns it 0.5, while the `visits_low`/`visits_high` scenarios only vary the earlier cycle's 4–9 category. They therefore do not resolve the 2013–14 boundary ambiguity. The question wording also changed between cycles. The approximate FI and adapted Fried scores, mobility skip, self-reported measures and selection all remain limitations.

A [same-sample boundary sensitivity](boundary_sensitivity/RESULTS.md) subsequently compared assignments of 0.5 and 1 for this 2013–14 category; the full FI–adapted-Fried correlation changed from 0.56416 to 0.56282. It does not recover the unknown individual counts or validate the source scoring.

Three prior follow-ups were also rerun with the corrected coding. The age/sex-adjusted full FI–adapted-Fried partial correlation was **0.540** (95% interval 0.488–0.592), versus 0.539 under the original coding. The reduced FI43–FFS2 correlation was **0.276** unadjusted and **0.250** adjusted. Removing the shared physical questions from both working scores yielded r=**0.239**; FI-minus-sleep-only versus the sleep item yielded r=**0.257**. Full corrected [age/sex](corrected_covariate_sensitivity/results.tsv), [shared-item](corrected_shared_items/correlations.tsv), and [sleep-question](corrected_sleep_specific/correlations.tsv) tables, scripts and logs are in their respective folders. These follow-ups remain post-hoc sensitivities, and FFS2 is not a validated Fried score.

The [pre-result plan](PROTOCOL.md), [original](score_original.R) and [corrected](score_corrected.R) scripts, full [original outputs](original/correlations.tsv), [corrected outputs](corrected/correlations.tsv), [paired comparison](comparison.tsv), [comparison code](compare.R) and [log](comparison.log) preserve both versions. The scripts differ in their output path and this single code mapping; source hashes were checked before the run. The same participant IDs and primary eligibility flags were checked in memory, the change occurred only at code 8, and no individual records were exported. These checks do not establish original authors' exact scoring or rule out other coding errors.

A [working FI–sleep covariance decomposition](sleep_covariance_decomposition/RESULTS.md) quantifies the direct included-item term using the corrected score and same survey domain. It is a mathematical decomposition, not a causal or genetic analysis.

A [shared-physical-item covariance decomposition](physical_covariance_decomposition/RESULTS.md) quantifies the two common FI question terms against the full adapted Fried score. It is an algebraic observational breakdown, not a causal or genetic result.

A [same-sample mobility-skip sensitivity](mobility_same_sample/RESULTS.md) compares assigning skipped walking items zero versus one in the same 2,986 eligible respondents. The working FI–adapted-Fried correlation changes materially, highlighting missing-data and shared-item sensitivity; the extremes are not formal bounds or a genetic analysis.
