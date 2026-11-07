# Paired comparison of one- and two-component genetic projections

Using the same 991,707 SNPs and 200 matched LDSC delete blocks, we compared each insomnia–frailty partial genetic correlation given ADHD with the corresponding value given ADHD and MDD. Differences are **two-component minus ADHD-only**:

| Frailty measure | ADHD-only | ADHD + MDD | Difference | Paired SE | Approximate 95% interval | Bonferroni p, three contrasts |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| FI | 0.378 | 0.200 | −0.178 | 0.030 | −0.237 to −0.119 | 1.07×10⁻⁸ |
| FFS | 0.165 | 0.073 | −0.092 | 0.020 | −0.131 to −0.054 | 7.09×10⁻⁶ |
| HFRS | 0.426 | 0.328 | −0.098 | 0.018 | −0.133 to −0.064 | 7.08×10⁻⁸ |

All three differences cross the stated correction for this **post-hoc three-contrast family**; this does not account for the broader research search. The point estimates use the rounded LDSC logs, so intervals and p-values are approximate. The two projections are different mathematical estimands, including different residual scaling. Their differences are not evidence that depression causes or mediates the sleep–frailty relationship, nor that adjusting for a person's depression diagnosis would have the same effect. MDD and the frailty measures can include sleep or fatigue symptoms, and unknown participant overlap, effective-N approximations, ancestry/reference fit and elevated intercepts limit inference. No difference identifies a shared causal variant, a gene-by-sleep interaction, or physical frailty independent of sleep and exhaustion.

The [protocol](PROTOCOL.md) was recorded before the contrasts were calculated. The [analysis](analyze.py), [results](results.tsv) and [receipt](receipt.json) retain the calculation and hashes. Both component projections and their individual jackknife SEs were independently reproduced from their matched delete vectors, and every four-trait point and delete-block matrix passed the positive-definiteness check (minimum eigenvalue 0.320). Identical corresponding trait h² vectors and common input signatures support block alignment. These checks verify the calculation, not the original GWAS, model calibration or causal assumptions. No AlphaGenome output was used as a model-training feature.
