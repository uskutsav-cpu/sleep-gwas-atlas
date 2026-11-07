# Clinical-insomnia outcome substitution: ML proxy stress check

The previously saved group-held-out bounded-ridge predictions for questionnaire insomnia were compared with the new clinical MVP insomnia LDSC outcomes, without retraining:

| Frailty target | Saved bounded prediction | MVP observed genetic correlation (SE) | Absolute difference |
| --- | ---: | ---: | ---: |
| FI | 0.581 | 0.553 (0.039) | 0.028 |
| FFS | 0.536 | 0.400 (0.035) | 0.136 |
| HFRS | 0.469 | 0.570 (0.031) | 0.100 |

The mean absolute difference across these three outcomes was 0.088 for bounded ridge, 0.106 for unbounded ridge, 0.139 for 3-nearest-neighbor and 0.373 for the training mean. These three dependent, already-inspected comparisons cannot establish a model ranking or independent predictive performance. They compare point estimates with uncertainty, and the model never received MVP-specific sleep–brain features: the questionnaire representation is reused as a proxy. The FFS and HFRS differences show why substituting a sleep source and treating it as a validated prediction is not justified. The source-specific GWAS, LD support, ascertainment and possible cohort overlaps differ. Both MVP and HFRS LDSC intercepts were elevated, and the frailty definitions retain sleep-related or exhaustion content.

The [protocol](PROTOCOL.md) was recorded before running the comparison. The [full results](comparison.tsv), [code](run.py), and [input hashes](receipt.json) preserve all four methods and the source checks; the new targets are documented in the [MVP LDSC sensitivity](../../../../../../replication/mvp_frailty_genomewide/RESULTS.md). This analysis is a post-hoc outcome-substitution stress check, not a new external model fit, clinical prediction, causal inference, or demonstration of physical frailty independent of sleep. No AlphaGenome output was used for model training.
