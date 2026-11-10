# CLSA sleep and frailty: wave-specific feasibility crosswalk

Checked 29 October 2026 against the published instruments and papers. These are candidate variables, not an extract, and their availability and coding in the current released dataset must be checked against the data dictionary and access agreement.

| Measure | Baseline Comprehensive | First follow-up Comprehensive | Interpretation |
|---|---|---|---|
| Sleep-onset difficulty | `SLE_30MIN_COM` | `SLS_30MIN_COF1` | Frequency of taking over 30 minutes to fall asleep in the past month. |
| Sleep-maintenance difficulty | `SLE_MIDFQ_COM` | `SLE_MIDFQ_COF1` | Frequency of waking in the night or early morning and having difficulty returning to sleep. |
| Sleep satisfaction | `SLE_QLTY_COM` | `SLE_QLTY_COF1` | Separate symptom domain; not interchangeable with onset or maintenance. |
| Frailty index and Fried-type phenotype | Both built from baseline measures in [Pasquet et al.](https://link.springer.com/article/10.1007/s40520-024-02706-w) | Follow-up item availability and reproducibility not established by that paper | Baseline cross-sectional comparison is plausible; a longitudinal comparison requires confirming repeated frailty components. |
| Genotypes | [Release documentation](https://www.clsa-elcv.ca/wp-content/uploads/2023/06/clsa_gwas_v3.pdf) describes a sample-ID linkage to project-specific survey/physical data | Same inherited genotypes may be linked if the requested data and participant consent permit | Neither actual linkage nor adequate marker coverage has been demonstrated for this project. |

Source instruments: [baseline Comprehensive in-home](https://www.clsa-elcv.ca/wp-content/uploads/2023/06/clsa_baseline_comp_in-home_v4.0_2018_jun08.pdf) pages 88–90; [follow-up 1 data collection site](https://www.clsa-elcv.ca/wp-content/uploads/2023/06/clsa_dcs_fu1_v1.5_2018jun07.pdf) pages 195–197. They use the same broad response frequencies for these two questions: 1 never, 2 less than once weekly, 3 once or twice weekly, 4 three to five times weekly, 5 six to seven times weekly; 8 and 9 denote unknown/no answer and refusal. Confirm stored codes and skip patterns in the exact release before implementation. Follow-up's onset variable begins **SLS**, not SLE. `SLE_STAYFQ` concerns staying awake and is not sleep-maintenance difficulty.

A candidate symptom endpoint, consistent with the frequency threshold in [Cavallo et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12068708/), is positive if either onset or maintenance is 4 or 5; negative if both are observed and each is 1–3; unknown otherwise. This is a symptom measure, not a clinical insomnia diagnosis. If one response is missing but the other is positive, the logical OR is positive; if one is missing and the other negative, it is unknown. Specify the missing-data rule and check exact operationalization against the paper before claiming replication.

## Design implications

The baseline instrument directly documents sleep questions, despite the Cavallo study describing their introduction at follow-up. Pasquet et al. built both frailty approaches at baseline in the same Comprehensive cohort, making a same-wave joint analysis a plausible access request. However, different study samples do not establish how many genotyped participants have complete sleep and frailty data, or that the requested data will be released. Baseline genotype-to-symptom and genotype-to-frailty associations and prospective analyses using later outcomes are distinct questions; temporal ordering does not eliminate confounding or establish causal direction.

Pasquet et al.'s FI uses 93 deficits, while their phenotype counts five physical-frailty components. The [supplement audit](clsa_fi_sleep_fatigue_overlap_audit.md) identifies a restless-sleep FI item and shared exhaustion, grip and walking inputs, plus a walking-cutoff inconsistency; define exclusion sensitivities with explicit denominators and missingness rules. Account for repeated measures, selection and attrition. Neither published measure has validated the project's scores.
