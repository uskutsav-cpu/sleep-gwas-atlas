# Framingham Offspring exam 8–9 variable crosswalk

Checked 29 October 2026. A feasibility crosswalk from the publicly available [2025a BioLINCC Offspring and Omni data dictionary](https://biolincc.nhlbi.nih.gov/media/studies/framoffspring/data_dictionary/Framingham_Offspring_2025a.pdf) and [exam 8](https://www.framinghamheartstudy.org/files/2017/08/Offspring-Exam-8.pdf) and [exam 9](https://biolincc.nhlbi.nih.gov/media/studies/framoffspring/Forms/Offspring_Cohort_Exam_9_Omni_1_Exam_4.pdf) forms. The downloaded dictionary's SHA-256 is `779fb00be2cfcc83d75daf1529d11845cb17df7a7681292baca026ebaac25f9f`. Its listing is public metadata, not access to participant data, and does not establish the exact dbGaP file, participant linkage, sample counts or marker quality.

| Assessment | Public variable(s) | Meaning and caveat |
|---|---|---|
| Exam 8 sleep | h730, h731, h732 | Trouble falling asleep, waking and difficulty returning to sleep, early waking. Form asks past-month frequency with codes 0–4; check release coding and missing values. |
| Exam 8 sleep | h726, h727, h733 | Workday sleep duration, time to fall asleep, excessive daytime sleepiness. Different constructs from insomnia diagnosis. |
| Exam 8 diagnosis | h745 | Ever told by health professional that participant has insomnia; form lists 0 no, 1 yes, 9 don't know. Do not treat 9 as a positive or negative response. |
| Exam 8 function | h605–h610, h611–h612 | Six grip trials and completion/reason variables. Need published phenotype's thresholds and treatment of missing tests. |
| Exam 8 function | h614–h623 | Walking aid, test completion, reasons and first, second and quick-walk times. Need distance and protocol before deriving speed. |
| Exam 9 diagnosis | j944 | Ever told by health professional of insomnia; form says 0 no, 1 yes, 9 don't know. This is the later exam, not the baseline exam 8 variable h745. |
| Exam 9 function | j765–j770, j771–j772 | Six grip trials in kilograms and test completion/reason variables. |
| Exam 9 function | j784–j793 | Walking aid, two walk times, laser times, completion/reason fields; exact measurement protocol and use in frailty definition require checking. |

The [SHARe dbGaP page](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000342.v23.p16) states that requests for individual-level phenotypes and molecular data for the top-level study and substudies go through `phs000007`; it lists three imputed panels and other genotyping resources. This narrows the request route but does not establish specific consent coverage or that these BioLINCC variable names or data are in the approved dbGaP release. Confirm current identifiers, sample overlap, imputation INFO, family structure and permissions before analysis.

For a prospective comparison, use baseline h730–h732 or h745 at exam 8 and subsequent frailty components at exam 9 if available. The exam 9 diagnosis j944 is contemporaneous with the later assessment, so do not substitute it for baseline sleep. Neither grip nor gait alone defines FFS; all five components, thresholds and missingness must be established. FI construction requires the complete deficit list and a prespecified sleep-free version. Sample selection and residual symptom overlap remain concerns. See [Framingham feasibility](framingham_no_fee_validation_route.md).
