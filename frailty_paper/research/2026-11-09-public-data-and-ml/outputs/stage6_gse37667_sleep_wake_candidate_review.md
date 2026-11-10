# Prolonged wakefulness: paired candidate-expression sensitivity

The public [GSE37667 study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE37667) measured whole-blood microarray expression in nine volunteers at baseline, after prolonged wakefulness, and after one night of recovery. It is a within-person sleep-loss experiment, **not** a chronic-insomnia case-control or frailty study. A [published data inventory](https://pmc.ncbi.nlm.nih.gov/articles/PMC13242190/) identifies the nine participants as male. We cannot use it to validate either genetic score or claim independent replication of the older-adult insomnia analysis.

We took the same four project-selected candidate genes previously screened in [GSE208668](stage6_gse208668_insomnia_candidate_review.md). Using the public GPL570 annotation, we retained probe sets mapping exactly to a single named gene and selected the probe with the largest overall mean across all 27 samples, without looking at group differences. For each gene we calculated the nine matched prolonged-wakefulness-minus-baseline differences and a paired t interval; BH adjustment covers these four tests only. Sample names were matched by volunteer number, not treated as 27 independent people. The submitted values were RMA-processed and ComBat-adjusted for array-processing day by the original submitters; this analysis does not reconstruct those processing steps. Recovery values below are descriptive means, without a separate hypothesis test.

| Gene | Selected probe | Wake minus baseline, log₂ | 95% CI | p | Four-gene BH q | Recovery minus baseline mean |
|---|---|---:|---:|---:|---:|---:|
| FOXP2 | 1564876_s_at | +0.063 | −0.086 to +0.213 | 0.359 | 0.507 | −0.053 |
| SMAD5 | 225219_at | −0.581 | −1.190 to +0.028 | 0.0588 | 0.235 | −0.056 |
| RBM6 | 201967_at | −0.109 | −0.470 to +0.252 | 0.507 | 0.507 | +0.153 |
| MEIS1 | 204069_at | −0.145 | −0.511 to +0.220 | 0.386 | 0.507 | −0.0076 |

None passed BH 0.05 in this small postselected screen. This does not show absence of an effect. There is no randomized parallel time control in this comparison, so time, intervention and other study effects cannot be separated; uncertainty includes assumptions of the small-sample paired t model and prior normalization. The nominal SMAD5 decrease differs in sign from the previous insomnia case-control contrast, but the experiments, tissues, arrays, probes and participants differ; do not call this a failed replication or infer statistically established heterogeneity. The selected probes have not been independently checked for transcript specificity or reliable expression above background.

Reproduce with `work/stage4_env/bin/python work/geo_sleep_wake/analyze_candidates.py`; [results](stage6_gse37667_candidate_paired.csv) and [checksum/method receipt](stage6_gse37667_candidate_paired.json) record the exact inputs. The analysis produces aggregate results only.
