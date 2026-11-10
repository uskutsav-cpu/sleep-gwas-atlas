# Post hoc transfer of the neuronal sequence baseline across immune conditions

We kept the existing baseline-trained sequence-only logistic model fixed and checked its ranking on the same held-out chromosome 1–2 variants that had usable paired results under IL-6 or interferon-alpha. The same variant must have a complete matching pair and a clear reported assay-detection classification in both conditions. No cue-condition results were used to train or tune the model.

| Matched comparison | Variants | Baseline AUROC | Cue AUROC | Cue average precision | Cue detected-call prevalence |
|---|---:|---:|---:|---:|---:|
| IL-6 | 136 | 0.552 | 0.516 | 0.231 | 0.176 |
| Interferon-alpha | 138 | 0.543 | 0.509 | 0.108 | 0.101 |

These are post hoc, small-sample descriptive results with no confidence intervals. Both cue AUROCs are close to 0.5; average precision depends on the class prevalence. They do not establish a difference in model performance between conditions. FDR-based calls differ for 14 of the 136 IL-6 pairs and 17 of the 138 interferon pairs, but a change in significance is **not** evidence of an interaction: a direct interaction test and uncertainty are required.

The same variants, study and donors underpin these comparisons, so they are not independent study or donor validation. A nonsignificant call is not proof of no biological activity, and these psychiatric-risk-selected reporter constructs do not provide sleep–frailty outcomes. No AlphaGenome outputs or derivatives were used, and source article reuse restrictions remain as described in the [original baseline review](stage6_neuronal_mpra_ml_baseline_review.md). The [plan](stage6_neuronal_mpra_cue_transfer_plan.md), [results](stage6_neuronal_mpra_cue_transfer.json), and [code](../work/neuronal_mpra_cue_transfer.py) record the matching, fixed training and evaluation. Independent genetic-score validation remains outstanding.
