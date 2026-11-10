# Open FRAILPOL data: exploratory physical-frailty benchmark

[FRAILPOL](https://doi.org/10.6084/m9.figshare.c.7874411) provides openly downloadable, CC BY 4.0 gait-sensor and clinical data on 668 community-dwelling older adults. The [dataset paper](https://doi.org/10.1038/s41597-026-06854-8) describes Fried-phenotype frailty classification. The published feature file contains one row per person: 342 robust, 279 prefrail and 47 frail, with 26 left/right gait-derived features. Figshare's MD5 for the three-class CSV matched the downloaded file (`72057ec58f8c648328344ff6e5eb67e5`). We did not find genome-wide genotype or insomnia measures in the gait feature file; the full participant workbook and all other files were not audited. This is **not genetic-score validation** or a sleep–frailty dataset.

An exploratory multinomial logistic-regression baseline used all 26 gait features, excluding the participant ID, with five stratified participant-level folds and scaling learned only in the corresponding training fold. The pooled out-of-fold results are:

| Measure | Result |
|---|---:|
| Overall accuracy | 0.599 |
| Balanced accuracy (mean class recall) | 0.453 |
| Majority-class accuracy | 0.512 |
| Recall: robust / prefrail / frail | 0.716 / 0.538 / 0.106 |

Only 5 of 47 frail participants were classified as frail; most were classified as prefrail. The majority classifier has balanced accuracy 1/3. Results are descriptive and the split was chosen after finding the dataset. This is one internal cross-validation run, not an independent external test, confidence interval, clinical model or evidence of utility. We did not compare against the authors' exact splits or establish a new state of the art. There was no hyperparameter search. The dataset is imbalanced and its data collection and exclusion procedures can affect generalization.

The outcome definition includes physical slowness, while the input is gait measured during a timed mobility task. These related measurements may share information; the result cannot establish independent prediction of broader frailty, causal pathways, sleep risk or genetics. Labels are created from the published phenotype, not the UK Biobank FFS. This is a useful public baseline for the mentor's ML direction only if later evaluated on an independent cohort with compatible measures and careful separation of participants and related data.

`Rscript work/frailpol_gait_baseline.R` reproduces the fixed seed, train-fold scaling, L2-penalized multinomial model and aggregate confusion matrix. [Metrics](stage6_frailpol_gait_baseline_metrics.csv) and [confusion matrix](stage6_frailpol_gait_baseline_confusion.csv) contain only aggregate counts; no participant predictions were retained. The dataset paper and Figshare attribution apply.

A subsequent [exploratory transfer](stage6_gstride_transfer_review.md) trained a different three-feature model on FRAILPOL and applied it to GSTRIDE. Accuracy was 0.387 and frail recall was zero; differences in task and features mean this is not a definitive external test of the original 26-feature model.

## Participant-workbook audit and label limitation

The public `participants.xlsx` has 771 numeric participant-ID rows and passed the publisher's MD5 check (`c438dc375281f2fb1e0afe08be89aed3`). Its headers describe demographics, cognition, five physical-frailty components, a field labeled `FI`, category target, gait and grip. That `FI` is the **sum of five physical phenotype components**, not the multi-deficit frailty index used in the project's FI GWAS. The paper's documentation describes the field named `Weakness` using the effort/getting-going question, an exhaustion component, and `Muscular strength` using grip; the labels should not be assumed to have the same meanings as in another cohort. No explicitly named insomnia or genotype columns were found. A limited English/Polish exact-term scan of text cells found no sleep-related matches; it is not a complete semantic review and cannot establish the absence of unlisted data.

Of 769 participant rows with numeric 0/1 values for all five components and numeric `FI`, the component sum agrees with `FI` for **763**, leaving six discrepancies. The source `Targets` category matches the standard 0, 1–2, and 3–5 categorization of the recorded `FI` for all **770** rows with both values. Matching published participant identifiers between the workbook and 668-row gait feature file yielded **668 matches**, all with the same recorded target category. Of 667 matched rows with complete numeric components and `FI`, five have component sums different from recorded `FI`; recomputing categories from the five component fields would differ from the supplied category for four. The reason is unverified: source coding or cleaning may explain it. Do not silently replace the published labels or treat the workbook's component sum as definitive.

The earlier baseline used the supplied gait-file labels, so its reported confusion matrix remains reproducible for those labels; this audit limits any claim that it was trained on categories independently reconstructed from every component. No case-level rows or participant identifiers are reproduced in the audit output. See the [aggregate audit](stage6_frailpol_participant_audit.json); reproduce using the bundled Python runtime and `work/frailpol_participant_audit.py`. Clarify the discordances with the dataset documentation or authors before using recomputed labels. The workbook still does not establish the genetics or sleep outcomes needed to validate the project scores.

## Sensitivity to the six ambiguous component records

As a post hoc label-quality sensitivity, I held the previously generated stratified five-fold assignment fixed, excluded the five gait records whose component sums disagree with the supplied `FI` and one record with an incomplete component, and compared two versions of the same exploratory model on the **same remaining 662 people**. Both use the publisher's target labels for those 662 people: the first trains each fold on all available other-fold rows, while the second also removes the six ambiguous records from its training folds. Feature scaling is fitted separately inside each training fold.

| Training set | Accuracy on common 662 | Balanced accuracy | Frail recall |
|---|---:|---:|---:|
| All other-fold rows | 0.600 | 0.455 | 5/46 (0.109) |
| Excluding six ambiguous rows | 0.597 | 0.452 | 5/46 (0.109) |

Eight of the 662 predictions changed, and neither run reported a fitting warning or nonconvergence. This narrow comparison suggests that excluding these records had little effect on the aggregate results in this specific split, but does not establish that the supplied categories are correct or that the model is robust in other data. The six rows were also excluded from the common evaluation, so this does not assess their individual predictions. The sensitivity is post hoc, not an independent test, clinical validation, sleep analysis, or genetic-score validation; the gait/physical-slowness overlap and poor frail-class recall remain. [Metrics](stage6_frailpol_label_sensitivity.csv), [original-training confusion matrix](stage6_frailpol_label_sensitivity_original_cm.csv), and [excluded-training confusion matrix](stage6_frailpol_label_sensitivity_excluded_cm.csv) contain aggregate data only. Reproduce with the bundled Python runtime and `work/frailpol_label_sensitivity.py`, which streams eligibility flags to the paired R analysis without saving participant identifiers.
