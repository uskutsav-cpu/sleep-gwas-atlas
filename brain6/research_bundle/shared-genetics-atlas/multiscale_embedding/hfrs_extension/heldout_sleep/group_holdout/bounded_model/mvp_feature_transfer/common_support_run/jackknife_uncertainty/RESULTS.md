# Paired block uncertainty for clinical-insomnia prediction residuals

On the 991,707 shared markers, we reran the nine LDSC pairs with 200 delete-block estimates and propagated the six feature estimates and each target jointly through the fixed bounded model. The residual is model prediction minus observed genetic correlation.

| Target | Residual | Conditional block-jackknife SE | Approximate 95% normal interval |
| --- | ---: | ---: | ---: |
| Frailty index | 0.112 | 0.045 | 0.024 to 0.199 |
| Fried frailty score | 0.213 | 0.039 | 0.137 to 0.290 |
| Hospital frailty risk score | −0.041 | 0.055 | −0.148 to 0.066 |

The calculation accounts for the covariance represented by jointly deleting the same genetic blocks from the six features and each target, **conditional on the fixed training data and model**. The FFS residual remains positive over this conditional normal interval. These are exploratory conditional uncertainty summaries, not full predictive intervals or independent validation: they omit training-set and hyperparameter uncertainty, source selection, cohort overlap uncertainty, possible model misspecification, and the broader post-hoc research search. They rely on LDSC's block-jackknife and normal approximations and do not establish statistical significance across the project or causation. The target FI/FFS/HFRS definitions still contain sleep- or exhaustion-related content, and elevated MVP and HFRS intercepts constrain interpretation.

All nine rerun r_g estimates match the earlier common-support run to printed precision. Each reconstructed jackknife SE matches its LDSC log to rounding precision, all delete-block heritabilities were positive and finite, and all nine MVP h² delete vectors were identical. The independently verified inputs have the same ordered SNP and allele-pair signatures; the runs used the same reference, common first trait and 991,707 valid variants. These checks support paired block alignment but do not validate the upstream GWAS or the analysis assumptions. The 22 LD-score files matched their preparation hashes immediately after the run began.

The [pre-result protocol](PROTOCOL.md), [runner](run.sh), [full LDSC log](ldsc-20261029T060147-2867.log), [delete values](.), [calculation](analyze.py), [results table](residual_uncertainty.tsv) and [receipt](receipt.json) document the calculation. The analysis uses a fixed model and no AlphaGenome output for model training.
