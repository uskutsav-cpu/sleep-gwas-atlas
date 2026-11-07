# Why the overall correlation does not identify the item-removed result

The existing clinical-insomnia–FFS estimate is approximately 0.395. That number alone cannot determine the corresponding genetic correlation after exhaustion is removed. This is a mathematical illustration, not a reanalysis of those participants or an estimate of any real item effect.

Suppose a standardized hypothetical physical component P and exhaustion component E have equal genetic variance and zero genetic correlation. Let the composite F = (P + E)/√2 and S represent sleep. The [illustration script](example.py) constructs two valid correlation matrices:

| Hypothetical r(S,E) | Implied r(S,P) | Same r(S,F) |
|---:|---:|---:|
| 0.2000 | 0.3586 | 0.3950 |
| 0.6000 | −0.0414 | 0.3950 |

Both matrices are positive definite, but their item-removed correlations have different signs. Real frailty score construction is more complicated: the components have unequal variances and correlations, may be binary or ordinal, and the FI GWAS used transformations. Knowing separate component summary statistics may help only if their sample, scale, covariates, missingness and cross-component genetic covariance are adequately characterized. Simply subtracting two reported correlations or their standard errors is not justified.

The appropriate direct test remains a same-participant, consistently processed GWAS of the original and item-deleted definitions, followed by a paired genetic-correlation comparison with uncertainty and appropriate multiple-testing control. The available physical-function composite and component GWAS are complementary sensitivities, not that matched comparison. See the [construction protocol](../CONSTRUCTION_PROTOCOL.md).
