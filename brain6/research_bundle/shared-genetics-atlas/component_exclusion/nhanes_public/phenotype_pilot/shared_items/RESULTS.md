# Shared physical-item sensitivity in the NHANES phenotype pilot

**Historical coding:** The tables below predate the [HUQ051 visit-code correction](../visit_code_audit/RESULTS.md). The corrected full-score estimate is 0.56416, and the corrected shared-item follow-up is linked there.

Among the same **2,433 respondents** in the [exploratory NHANES pilot](../RESULTS.md), the correlation between the approximate full FI and adapted Fried score was 0.563. When both scores were changed to remove two shared physical questions, the correlation was 0.237. This shows sensitivity to how the scores are constructed; it does **not** establish that the difference is caused solely by shared questions.

| Working comparison on the same sample | Weighted Pearson r | 95% survey interval |
| --- | ---: | ---: |
| Full approximate FI vs full adapted Fried score | 0.563 | 0.513–0.613 |
| FI without quarter-mile and lifting questions vs adapted score without those questions | 0.237 | 0.172–0.302 |
| FI additionally without sleep and fatigue vs that three-component adapted score | 0.236 | 0.172–0.300 |
| FI without sleep, fatigue, those physical questions, BMI and grip vs adapted score without the shared questions and BMI/weight component | 0.275 | 0.225–0.325 |

The last adapted score contains only room-walking difficulty and low activity. It is a **two-component research sensitivity**, not a validated Fried frailty score. The full-versus-last correlation difference is 0.288 (95% paired survey interval 0.251–0.325). Item removal changes content, denominators, score scale and reliability together, so the paired difference cannot isolate measurement overlap or establish a causal mechanism. The retained FI remains a broad, unvalidated deficit score.

The first removal targets PFQ061B and PFQ061E, which are present in both working scores. The three-component adapted score still contains a BMI/weight component; the FI contains BMI and uses BMI in its grip criterion. The last sensitivity also removes those FI components and the adapted score's BMI/weight component to reduce direct BMI-input overlap. Other related physical questions and common influences remain. The same primary cohort excludes some participants with mobility-related skipped responses, and the pilot's phenotype definitions and coding remain unvalidated. These intervals are post-hoc, uncorrected for the wider project, and describe phenotypes rather than genetic correlations or independent replication.

The [pre-result plan](PROTOCOL.md), [analysis](analyze.R), [all estimates](correlations.tsv), [paired contrast](contrast.tsv), [sample count](counts.tsv) and [successful run log](retry.log) preserve the calculation. Source-file hashes were checked against the existing receipts. Direct weighted correlations reproduced the survey estimates, and a finite-difference gradient reproduced the paired standard error. A first attempt stopped because a temporary variable name shadowed R's `c()` function ([log](run.log)); the corrected attempt produced the reported results. No person-level derived records were exported. See the parent [scoring crosswalk](../SCORING.md) and [mobility sensitivity](../RESULTS.md) for the working assumptions and limitations.
