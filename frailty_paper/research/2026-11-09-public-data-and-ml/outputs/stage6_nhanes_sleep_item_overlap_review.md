# Illustration: putting sleep in the frailty outcome

Using public NHANES examined adults aged 60+ in the **same participants within each cycle**, we compared two exploratory adapted indices. The first has 24 deficit items and excludes sleep and fatigue. The second adds the sleep question and a fatigue question, for 26 items. Both use age/sex-adjusted survey-weighted models. The exposure is reporting having ever told a health professional about trouble sleeping.

| Cycle | n | 24-item index, sleep/fatigue excluded (95% CI) | 26-item index, sleep/fatigue included (95% CI) |
|---|---:|---:|---:|
| 2011–12 | 1,553 | 0.071 (0.051–0.092) | 0.112 (0.091–0.133) |
| 2013–14 | 1,651 | 0.090 (0.064–0.116) | 0.132 (0.106–0.157) |

**The included-item analysis is intentionally circular:** the same yes/no sleep response is both the exposure and a component of its outcome. An index containing it can show a larger association partly for arithmetic reasons. Fatigue was added at the same time, and varying denominators and other deficits prevent interpreting the observed difference as the isolated effect of sleep or fatigue. The sleep/fatigue-excluded analysis still has a positive association in this sample, but other cognitive, mood and physical items can relate to sleep. Neither adapted index has been validated here, and sleep discussion is not a clinical insomnia diagnosis.

This is a post hoc, cross-sectional illustration—not an analysis of the original Atkins genetic FI, a formal test of changes in genetic correlation, evidence of causation, or independent validation. A genetic comparison needs matched same-person GWAS for the original FI and prespecified item-excluded variants; the published Mak mental-wellbeing-stripped GWAS is a potential lead but its full summary statistics have not been located ([source audit](mak2023_sleep_excluded_fi_lead.md)). The [plan](stage6_nhanes_sleep_item_overlap_plan.md), [results](stage6_nhanes_sleep_item_overlap_results.csv) and [R code](../work/nhanes_sleep_item_overlap.R) record the definitions and MEC survey design.
