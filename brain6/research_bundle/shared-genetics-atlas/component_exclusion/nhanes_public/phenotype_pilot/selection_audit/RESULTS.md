# Complete-score selection audit

The primary pilot includes 2,433 of 3,472 examined respondents aged 60 or older with positive examination weights: **75.96% weighted** (95% survey interval 73.10%–78.83%). The remaining 1,039 were excluded by the complete-score requirements. This post-hoc descriptive audit compares these groups; counts are unweighted and estimates use the survey design.

| Measure | Included | Excluded |
| --- | ---: | ---: |
| Mean age, top-coded at 80 | 68.84 (n=2,433) | 72.19 (n=1,039) |
| Female | 53.47% (n=2,433) | 59.53% (n=1,039) |
| Difficulty walking without special equipment because of a health problem | 0.00% (n=2,433) | 71.08% (n=1,039) |
| Any recent DPQ030 sleep symptom, among valid responses | 32.24% (n=2,433) | 42.64% (n=779) |
| Ever told a health professional about trouble sleeping, among valid responses | 29.40% (n=2,433) | 36.34% (n=1,037) |

PFQ054 asks about **difficulty walking without equipment**, not whether a person actually uses it. The CDC questionnaire skips the quarter-mile walking question for affirmative responses ([2011–12](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/PFQ_G.htm); [2013–14](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/PFQ_H.htm)). The primary scoring requires that question, so the zero in the included group reflects this eligibility rule; its zero linearized standard error is a boundary artifact, not proof of zero population uncertainty. The earlier mobility sensitivity assigns a deficit for the skipped item, which is an additional assumption and changes the people included.

Completion requires a valid DPQ030 response, while 260 excluded respondents lack one, so the sleep-symptom comparison is conditional on responding and cannot describe all excluded people. Disclosure is a different lifetime question, not an insomnia diagnosis. These comparisons demonstrate selection and differences in observed characteristics, but do not quantify bias in the correlation or establish a causal relationship. The FI and adapted Fried score are provisional reconstructions, and the analysis is restricted to the noninstitutionalized NHANES examination population.

The [pre-result plan](PROTOCOL.md), [analysis code](analyze.R), [inclusion estimate](inclusion.tsv), [group estimates and 95% intervals](groups.tsv), and [run log](run.log) record the analysis. Survey estimates were checked against direct weighted calculations; all 32 input files matched their saved hashes. This does not independently validate score construction. No participant-level output was saved.
