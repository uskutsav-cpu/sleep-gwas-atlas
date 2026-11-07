# Exploratory differences between physical-function correlations

The source physical-function comparison did not save paired block-jackknife estimates. As a **post-hoc sensitivity**, we calculated a deliberately conservative standard-error bound for differences between its three frailty correlations. For two estimators A and B, Cauchy–Schwarz gives SD(A−B) ≤ SD(A)+SD(B), whatever their covariance. We use the two reported LDSC SEs as approximations, add a small allowance for four-decimal rounding, and use a normal approximation. The SEs are estimates, not guaranteed upper bounds on the true sampling uncertainty; this calculation cannot replace a paired jackknife or correct model violations.

| Contrast in r_g with physical function | Difference | Approximate 95% interval using maximum SE | Conservative approximate p |
|---|---:|---:|---:|
| FI − FFS | 0.2214 | 0.1386 to 0.3042 | 1.57×10⁻⁷ |
| FI − HFRS | −0.0930 | −0.1758 to −0.0102 | 0.0277 |
| FFS − HFRS | −0.3144 | −0.3888 to −0.2400 | 1.11×10⁻¹⁶ |

The FI–FFS and FFS–HFRS differences pass a post-hoc three-contrast benchmark of p < 0.0167 even under this maximal-SE calculation; FI–HFRS does not. The broader project-wide testing burden is not controlled. All contrasts share the physical-function source, and possible sample overlap, inflated or large intercepts, different phenotype definitions, cohort and ancestry differences and estimation error constrain interpretation. Do not attribute these differences to removing exhaustion, causality, or a formal comparison of individual risk. The logged SEs are rounded, and exact retained marker lists and paired delete-block estimates were not saved.

Inputs are the [audited completed LDSC results](../audit-20261029T161051-7498.json); the [script](calculate.py) and [output](summary.json) preserve the arithmetic. The original per-trait diagnostic details are in the [parent results](../RESULTS.md). A properly designed comparison should retain aligned delete-block estimates and predefine the full comparison family.
