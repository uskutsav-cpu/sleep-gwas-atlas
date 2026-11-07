# FI–ADHD posterior sensitivity: effect and uncertainty perturbations

At the model's leading FI–ADHD shared-signal candidate rs9969232, the aligned FI association statistic is Z=4.71 in the meta-analysis input and Z=5.51 in the UKB-only input. By comparison, at meta-FI's top-weighted rs936145 these are 6.09 and 6.28; at UKB-FI's top-weighted rs2396766 they are 6.09 and 6.34. Meta-analysis FI betas and standard errors are rounded to four decimal places. Z values are descriptive and the studies are not independent; no formal test of their differences is implied.

To inspect mathematical sensitivity only, we replaced one input column at a time across the fixed 1,055 SNPs. At p12=10⁻⁵ the FI–ADHD H4 values were:

| FI beta input | FI standard-error input | H4 |
|---|---|---:|
| Meta-analysis | Meta-analysis | 0.0441 |
| UKB-only | Meta-analysis | 0.2719 |
| Meta-analysis | UKB-only | 0.0698 |
| UKB-only | UKB-only | 0.3689 |

The middle two rows are **artificial hybrids**, not valid GWAS fits: combining betas and standard errors estimated in different samples has no established sampling interpretation. They show that the mathematical posterior can move when either input changes, especially when betas are swapped in this particular exercise. They do not separate causal contributions, estimate uncertainty, resolve cohort orientation, or indicate which source is correct. These are post hoc changes on a region selected using related data; effects and standard errors are coupled, and the nonlinear model makes effects of swaps context-dependent. Both real fits retain all the parent limitations, including overlapping samples, FI sleep items, uncertain study dependence, single-signal assumptions, and untested region boundaries.

The [reproducible calculation](analyze.py), [full table](counterfactual.tsv) and [input receipt](receipt.json) also cover p12=10⁻⁶ and FI–insomnia. The implementation reproduces H3 and H4 for all corresponding non-hybrid source fits; this is an arithmetic check, not a calibration or validation test. Restricted per-variant input remains on the SSD.
