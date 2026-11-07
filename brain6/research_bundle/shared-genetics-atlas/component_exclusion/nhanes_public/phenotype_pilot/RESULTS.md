# Exploratory NHANES frailty and sleep phenotype comparison

**Correction:** The figures below use an earlier visit coding that mistakenly treated HUQ051 code 8 as missing. A [corrected rerun](visit_code_audit/RESULTS.md) gives full FI–adapted-Fried r=0.56416 and retains the same 2,433-person primary sample. Corrected sleep, shared-item and age/sex follow-ups are also linked there; the tables below preserve historical results.

In a **provisional, survey-weighted observational pilot**, an approximate FI and an adapted Fried score were positively correlated in the same NHANES respondents. The analysis uses public 2011–2014 data, not the UKB GWAS sample. It is a working reconstruction with unresolved scoring choices, not a validated reproduction of either published instrument.

The primary complete-score domain contains **2,433 respondents aged 60 or older**. The full mobile-exam survey design contains 19,151 participants before domain subsetting; combined-cycle weights, cycle-specific strata and nested sampling units were used. Confidence intervals below use survey linearization and 32 design degrees of freedom. Counts are unweighted; correlations and intervals are weighted.

| Same-person comparison | Weighted Pearson r | 95% interval |
| --- | ---: | ---: |
| Full approximate FI and adapted Fried score | 0.563 | 0.513 to 0.613 |
| Full FI and sleep symptom | 0.311 | 0.262 to 0.359 |
| FI without sleep and low-energy items and sleep symptom | 0.238 | 0.188 to 0.287 |
| Modified FI and adapted Fried score without its quarter-mile proxy | 0.433 | 0.377 to 0.488 |

The paired difference between the two FI–sleep correlations is **0.073** (95% interval 0.066 to 0.080) in the same participants. The full FI literally contains the sleep response, so a higher correlation is partly expected from shared measurement; removing both sleep and fatigue also changes the FI denominator and construct. The result does not isolate the contribution of either item, establish causality or show a purely physical mechanism. The modified FI still contains nonphysical deficits and shares other functional questions with the adapted Fried score. These intervals are exploratory and do not correct for the wider project or establish equivalence or clinical thresholds.

A [post-hoc sleep-specific follow-up](sleep_specific/RESULTS.md) separately removes the sleep item and then low energy, and compares a different question about disclosure of sleep trouble to a health professional. It remains the same cohort and provisional scores, not an independent replication or causal decomposition.

A [shared-physical-item sensitivity](shared_items/RESULTS.md) finds lower correlations after removing common physical questions on the same sample. Those changes also alter score content and scale, so they do not identify the isolated effect of item overlap.

A [same-sample age and sex sensitivity](covariate_sensitivity/RESULTS.md) finds an adjusted partial correlation of 0.539 for the full working scores. This is observational and does not remove the sample-selection or instrument limitations.

The primary analysis leaves the quarter-mile item missing for people skipped because of walking-equipment needs. A [separate sensitivity](correlations.tsv) assumes those skipped responses indicate a deficit in the affected walking questions; it includes **2,986** people and yields full FI–adapted-Fried r=0.715. Both the scoring **and the included people** change, so this is not a causal test of the skip assumption or a directly comparable effect difference. Assigning the ambiguous earlier-cycle 4–9 healthcare visits to either adjacent FI category gives full FI–adapted-Fried correlations 0.564 in the same primary-sized sample. Weight-loss intentionality can also be unknown; the primary coding leaves that component missing unless BMI itself gives a positive result.

A [selection audit](selection_audit/RESULTS.md) shows that the complete-score sample excludes all respondents reporting difficulty walking without special equipment, because the questionnaire skips the required quarter-mile question for them. Included and excluded respondents also differ in age and observed sleep responses. This limits the interpretation of the primary correlation.

The previously re-executed project **genetic** FI–FFS correlation is 0.7534 (SE 0.0199) from UKB-containing GWAS, documented in the [global audit](../../../global_audit/fi_ffs_rerun/RESULTS.md). It is a different estimand, cohort and pair of instruments; it should not be compared numerically as if it were the NHANES phenotypic correlation. This pilot cannot provide a sleep-item-excluded GWAS or independent genetic replication.

## Methods, checks and limitations

The [protocol](PROTOCOL.md) defines the analysis; [scoring](SCORING.md) lists all 49 working items, score construction and assumptions. [Code](score_and_analyze.R), [sample counts](counts.tsv), [all correlations](correlations.tsv), [paired contrasts](paired_contrasts.tsv), [run log](verified_run.log), and [input verification](input_verification.json) preserve the audit. All 32 public CDC files matched their recorded SHA-256; there was no independent publisher checksum. Survey-weighted correlation estimates were independently checked against the direct weighted formula, and standard errors were checked against a finite-difference gradient and the survey covariance matrix. Score-range and item-deletion arithmetic checks passed. These checks do not prove the original published scoring or the validity of the chosen missing-data assumptions.

The FI table's merged cells, grip trial choice, visit grouping, diagnosis coding, survey skip patterns and missingness remain material. The alternative score uses a quarter-mile walking question as an exhaustion proxy, not the project FFS exhaustion item. The primary sample excludes some people who need walking equipment and requires at least 40 FI items and complete adapted-Fried components; selection may bias associations. The approach has not been externally validated, and no item-excluded genetic signal or causal gene follows. First install attempt for the survey library failed for lack of space; [retry log](../item_mapping/survey_retry.log) records the successful package installation and the final run records versions. No derived participant-level rows were exported.
