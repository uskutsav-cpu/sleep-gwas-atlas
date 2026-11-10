# LASI wave 1: candidate instrument items for sleep and physical frailty

The official [LASI wave 1 individual questionnaire](https://iipsindia.ac.in/sites/default/files/LASI_W1_Individual_Questionnaire.pdf) identifies the following **instrument question IDs**, not necessarily the variable names in an authorized raw or harmonized data release. Match them against the actual release and codebook before use. No participant data or score was analyzed.

| Instrument IDs | Measure | Relevance and limits |
|---|---|---|
| `HT219`, `HT220`, `HT221`, `HT222` | Prior-month sleep initiation, nighttime awakening, early awakening and daytime unrestedness | The published FI's sleep deficit uses all four; the published sleep–physical-frailty study uses the first three. Frequency cutoffs differ between those studies. These are symptom reports, not a clinical diagnosis. |
| `FS703`, `FS708` | Tired/low energy and everything was an effort during the past week | Candidate exhaustion items in the adapted physical phenotype. The instrument response categories refer to days in the past week; verify the published paper's precise combined rule and missingness. |
| `HB211` | Frequency of vigorous activity | Candidate low-activity item in the adapted physical phenotype. |
| `FO234` | Weight loss attributed to not enough household food | The published physical-phenotype proxy is narrower than unintentional weight loss from any cause. |
| `BM026`, `BM028`–`BM031` | Dominant hand, two left and two right grip trials | The questionnaire records the two hands separately. Its hand-injury safety skips and equally dominant hand response require explicit treatment. The two published frailty constructions use different reported grip summaries. |
| `BM056`, `BM057`, `BM058`, `BM059` | Two walking times, walking aid and compliance | Candidate measured 4-m walking fields. Check eligibility and actual data-release coding. The published FI and physical phenotype apply different thresholds. |

For the sleep questions the instrument uses 1=never, 2=rarely, 3=occasionally and 4=frequently. The [sleep/physical-frailty study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10339638/) describes coding occasional or frequent on any of `HT219`–`HT221` as a positive symptom composite. The [40-item FI supplement](https://pmc-oa-opendata.s3.amazonaws.com/PMC10391831.1/bmjopen-2023-071842supp001.pdf) describes code 4 on any of `HT219`–`HT222` as its sleep deficit. Handling partially missing responses is not fully specified by those descriptions; do not treat a missing question as a negative response without a prespecified, source-supported rule.

**Critical sentinel:** The instrument instructs interviewers to record `993` for a grip or walking attempt that the respondent could not complete **or** an unresolvable equipment problem. This value is not a plausible kilogram or seconds measurement and must not be averaged as an ordinary numeric result. It combines two potentially different reasons; the instrument alone cannot resolve whether inability should count as frailty or missingness in a particular released dataset. Check the actual release's recodes, accompanying safety/skip fields and published analytic implementation before defining a rule.

The 40-deficit FI includes many further diagnoses, symptoms and functional items that are not mapped here. The [LASI FI and physical-phenotype comparison](stage6_lasi_fi_ffs_sleep_overlap_review.md) describes the direct sleep overlap, differing grip and walk thresholds, and unresolved cutoff rules. This map is a concrete request and audit aid, not a validated executable recipe or evidence that the required fields can be linked to NIAGADS genotypes.

Local source `work/lasi/wave1_individual_questionnaire.pdf`, SHA-256 `4045852d209b537318d0283f1475d322ae03ec451498ea15d6aba1242738a0aa`.
