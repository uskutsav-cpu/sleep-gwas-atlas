# Sensitivity to indicator-specific genetic covariance

The factor model assumes insomnia has no genetic covariance with physical-indicator residuals. Three post-hoc relaxed models each allow one indicator's residual to correlate with the residual component of both insomnia sources. Loadings and both factor correlations are refitted jointly using the same nine saved moments and 200-block sampling covariance.

| Residual link allowed | Registry factor correlation (jackknife SE) | Questionnaire factor correlation (jackknife SE) |
|---|---:|---:|
| None: restricted joint model | −0.0870 (0.0449) | −0.4011 (0.0320) |
| Walking | −0.1439 (0.0708) | −0.3742 (0.0801) |
| Grip | −0.0623 (0.0376) | −0.3193 (0.0689) |
| Activity | −0.1189 (0.0497) | −0.4104 (0.0326) |

All three scenarios completed all 200 delete refits without accepted boundary estimates. Five starts per full model agreed, and a second optimizer reproduced each best objective within 1.3×10⁻¹². Analytic gradients passed finite-difference checks. These scenarios show dependence on the zero-residual restrictions, especially for registry insomnia; they do not identify the correct biological model. Their fit diagnostics and source contrasts are in [the output table](joint_fit/residual_diagnostics.txt), with [parameters](joint_fit/residual_parameters.txt), [protocol](RESIDUAL_PROTOCOL.md) and [code](residual_fit.R). The sampling-covariance weights remain fixed, and all original source, overlap, reference-fit and post-selection limitations apply.

## Allowing all indicator-specific links: lack of identification

Fixing the physical loadings at their original physical-only values, let rho be the factor–insomnia genetic correlation and r_i the observed insomnia–indicator correlations. A covariance model allowing all three indicator-specific components to share genetics with insomnia is positive semidefinite when

`rho² + sum_i [(r_i − loading_i × rho)² / (1 − loading_i²)] ≤ 1`.

Solving this quadratic gives conditional mathematical compatibility ranges of **[−0.747, +0.628] for registry insomnia** and **[−0.864, +0.456] for questionnaire insomnia**. These are not confidence intervals, biological plausibility ranges or estimated effects. They use rounded point moments and omit sampling uncertainty. Their purpose is to show that positive, zero and negative latent-factor correlations can reproduce the same physical and insomnia correlations once all residual links are free. The observed summaries alone therefore cannot identify that latent-factor association without additional restrictions.

The [calculation](identification_range.py) constructs explicit positive-semidefinite Gram matrices for both endpoints, zero and positive/negative interior examples, verifies all diagonal variances and observed correlations, and retains [the witnesses](joint_fit/identification_ranges.json). This is an identification argument, not a claim that the alternative models are equally biologically plausible.

A more transparent alternate estimand is the [explicit equal-weight genetic physical-function index](INDEX_RESULTS.md). Its definition fixes the combination rather than attributing covariance to an unobserved factor; it still requires valid source GWAS and cannot replace the requested same-person item-deleted frailty GWAS.
