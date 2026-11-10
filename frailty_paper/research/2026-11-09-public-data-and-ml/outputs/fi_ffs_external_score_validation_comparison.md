# External FI and FFS score evidence: what can be compared?

The strongest fair conclusion is **both broad scores show associations with their respective frailty outcomes in cohorts outside their discovery samples; these reports do not establish which score predicts better**. Their outcomes, cohorts, genetic-score methods and metrics differ. Neither evaluates the project's 15-marker candidate set or a shared insomnia–frailty score.

| | FI-derived score | FFS-derived score |
|---|---|---|
| Source | [Flint et al., 2024](https://www.nature.com/articles/s41598-024-63229-y) | [Ye et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/) |
| Discovery | UK Biobank component of Atkins FI GWAS, n=164,610, ages 60–70 | UK Biobank FFS GWAS, n=386,565, mean age 57 |
| Evaluation | LBC1936 and ELSA; FI, with 30 and 62 deficits respectively | HRS, n=9,720, mean age 73; FFS frailty outcome |
| Construction | Clumping and several p-value thresholds, with the best-fitting threshold selected in the evaluation analyses | PRS-CS score using 890,487 available SNPs in HRS |
| Reported evaluation | Extra FI variance explained above age, sex and PCs: 2.1% LBC Wave 1 and 1.8% ELSA age group 1 | Score alone AUC 0.57; score plus age, sex and PCs AUC 0.72 |
| Main limits | Target thresholds selected from multiple options; possible selection optimism. ELSA outcome includes restless sleep. | Different outcome/metric; standalone performance is modest; source's association and sign-test statistics have reporting inconsistencies documented in `ye_hrs_score_validation_audit.md`. |

The FI GWAS and external target cohorts are described as independent, but selecting the highest-performing score threshold among five in the evaluation data is **not equivalent to evaluating a single frozen threshold in a fresh test set**. This may inflate the reported best results. The paper identifies overlapping longitudinal participants and attrition, so the five ages should not be treated as five independent replications. The ELSA and LBC cohorts are also used in Foote et al.'s later validation; their results should not be counted as independent cohorts again.

I checked [Flint Supplementary Tables S7–S9](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-024-63229-y/MediaObjects/41598_2024_63229_MOESM1_ESM.docx): the discovery FI explicitly includes insomnia and recent tiredness/lethargy; ELSA includes restless sleep in the past week. No explicit sleep, insomnia, fatigue, tiredness, exhaustion or energy keyword appears in the LBC item table, although that does not rule out indirect sleep relationships. The supplement SHA-256 is `93d1b21a0bbb5c2dc648d699e06c46583b0367abb62e5b53b2a09bd478c6ea19`.

There is also an internal reporting ambiguity at the last age: the paper's results text describes both final-wave/group associations as p < 0.1, whereas its abstract and discussion describe significant associations at all five times. Exact last-wave estimates, test definitions and adjustment should be checked in the original data before claiming an association at the conventional 0.05 threshold for every time point. Do not infer a reliable age-dependent decline from these selected results alone.

For a fair score comparison, apply frozen FI and FFS scores, with harmonized alleles and documented missing-variant handling, in the same independent participants and evaluate the same prespecified outcome and time horizon. Use a paired uncertainty estimate for their performance difference and evaluate both calibration and discrimination. If developing the proposed 15–30-marker score, reserve independent samples for selection/training and final evaluation, and compare against these full scores and a clinical baseline model.
