# Second public muscle-expression dataset: candidate follow-up

We examined the baseline portion of [GSE117525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117525), an Affymetrix RMA expression dataset of vastus-lateralis biopsies. The 61 (pre-)frail and 73 healthy older baseline samples have 134 distinct subject IDs. We excluded all young and post-training measurements, avoiding repeated-measure inflation. We evaluated the same four postselected genes as in the [GSE144304 RNA-seq check](stage6_gse144304_muscle_candidate_expression_review.md), using the GEO GPL20880 Entrez-based probe annotation. An exploratory OLS model adjusted for age and sex with HC3 standard errors; the BH column covers only these four selected genes, not the wider project search.

| Gene | Adjusted difference in processed expression, frail minus fit | HC3 SE | Nominal p | BH q within four genes |
|---|---:|---:|---:|---:|
| FOXP2 | +0.001 | 0.066 | 0.985 | 0.985 |
| SMAD5 | +0.137 | 0.040 | 0.00093 | 0.0037 |
| RBM6 | −0.096 | 0.035 | 0.0064 | 0.0128 |
| MEIS1 | −0.063 | 0.054 | 0.239 | 0.319 |

**These are fragile observational signals, not genetic validation.** Sample titles for all 61 frail baseline samples contain `PROM`, while none of the 73 healthy older baseline titles does. The meaning of these title labels has not been fully verified, but their perfect association with group raises a serious possibility of study or recruitment confounding. A [published reanalysis](https://academic.oup.com/qjmed/article/118/9/657/8127924) also reports baseline batch differences. Frailty, study provenance, collection and processing effects cannot be disentangled from these data on the current evidence. The microarray measures bulk expression, not causal variant action or a cell type.

As a limited sensitivity, we additionally adjusted for sample-title prefixes G084–G087; their exact relation to technical batch is unverified and they are **not** the PROM indicator. SMAD5 and RBM6 associations remained nominally detectable (respectively estimates +0.169, p=0.00020; −0.088, p=0.020), but this does not resolve the perfectly group-aligned title label or unknown confounding. Results in GSE144304 had no detected difference for any of the four genes; they do not independently replicate these associations, though their small sample and different RNA-seq platform limit conclusions. Participant or close-relative overlap between the studies has not been verified. Do not combine estimates across platforms or claim the absence of effects.

Code: `work/geo_muscle/analyze_gse117525_candidates.py` and `work/geo_muscle/analyze_gse117525_plate_sensitivity.py`; aggregate results: `outputs/stage6_gse117525_candidate_expression.json` and `outputs/stage6_gse117525_plate_sensitivity.json`. The [source audit](stage6_gse144304_muscle_candidate_expression_review.md) notes sample counts and source files. These analyses neither test insomnia nor evaluate the genetic scores.

## Within-participant follow-up

We also identified 31 frail/pre-frail and 41 healthy older individuals with both baseline and post-training biopsies. For each participant, we calculated post-minus-baseline expression and compared changes between groups with baseline age and sex adjustment and HC3 standard errors. This avoids treating the two measurements from the same participant as independent, but does not remove the study/recruitment differences between groups.

| Gene | Adjusted difference in change, frail minus fit | HC3 SE | Nominal p | BH q within four genes |
|---|---:|---:|---:|---:|
| FOXP2 | +0.216 | 0.135 | 0.114 | 0.227 |
| SMAD5 | +0.002 | 0.080 | 0.980 | 0.980 |
| RBM6 | +0.014 | 0.059 | 0.811 | 0.980 |
| MEIS1 | +0.195 | 0.080 | 0.0175 | 0.070 |

MEIS1's nominal result does not pass BH correction even within these four postselected genes; the wider project search is not accounted for. The other estimates are uncertain and do not prove equal changes. This is an observational group comparison of paired changes, not a randomized test of exercise: the original study describes different training frequencies for frail and healthy participants (two versus three sessions per week), and follow-up selection, baseline study provenance, sample processing and cell composition may differ. Both groups received training, so a time effect cannot be separated from exercise. There is no insomnia phenotype, genetic-variant analysis or score validation here.

Reproduce with `work/geo_muscle/analyze_gse117525_paired_exercise.py`; results: `outputs/stage6_gse117525_paired_exercise.json`.
