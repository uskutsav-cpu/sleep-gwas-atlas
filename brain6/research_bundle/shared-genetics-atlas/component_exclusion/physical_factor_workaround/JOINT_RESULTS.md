# Joint physical-factor sensitivity

On 6 November 2026, we fitted all nine saved genetic correlations jointly, estimating three shared positive physical loadings and separate correlations with registry and questionnaire insomnia. The nine-moment sampling covariance was reconstructed from the same 200 genomic-block delete estimates; it was positive definite with condition number 6.13. Source log hashes matched the saved evaluation receipts, all logs retained 1,053,614 markers, and per-trait heritability delete vectors matched across pairs. The SSD was disconnected, so this verifies saved summaries and prior source receipts rather than rechecking the original GWAS files.

| Parameter | Joint estimate | Block-jackknife SE |
|---|---:|---:|
| Walking loading | 0.6885 | 0.0538 |
| Grip loading | 0.4396 | 0.0417 |
| Activity loading | 0.1658 | 0.0346 |
| Registry insomnia–physical factor | −0.0870 | 0.0449 |
| Questionnaire insomnia–physical factor | −0.4011 | 0.0320 |

Unlike the earlier two-stage projection, each of 200 delete fits re-estimated all five parameters. The registry estimate's unadjusted normal-approximation 95% interval is approximately **−0.175 to +0.001**, so it does not exclude zero. This weakens the registry evidence relative to the preliminary −0.112 projection. The questionnaire estimate remains strongly negative under the fitted model, but its sample shares UK Biobank with the physical indicators. The paired questionnaire-minus-registry contrast was −0.3141 (SE 0.0483; unadjusted normal p≈8.1×10⁻¹¹). Different insomnia definitions, populations, ascertainment, shared samples and reference fit change together; the contrast cannot identify a cause and is not independent replication.

The fitted quadratic statistic was 6.95 on approximately four degrees of freedom (p≈0.139). Non-rejection does not validate the model, its residual-covariance assumptions or the physical phenotype. Activity's loading remains small. These are nonlinear GLS fits to a vector of genetic correlations, using its estimated sampling covariance; they are not full multivariate LDSC/Genomic SEM or a factor GWAS. The weighting covariance is fixed across delete refits, so its estimation and model-selection uncertainty are not fully propagated. Point correlations were rounded to four decimal places in the original logs. All earlier source, overlap, LD-reference and phenotype limitations remain.

The model assumes insomnia is uncorrelated with each indicator's specific residual genetic component. All loadings and factor correlations stayed within bounds in the accepted full and 200 delete fits. Four of five L-BFGS-B starts reached the same best objective; one converged to a poor boundary solution. An independent `nlminb` optimizer reached the better solution from all five starts and agreed on its objective within 1.1×10⁻¹² and parameters within 4.3×10⁻⁸. Analytic gradients matched finite differences within 1.5×10⁻⁹. These are numerical checks, not statistical or biological validation.

The [protocol](JOINT_PROTOCOL.md), [input reconstruction](joint_prepare.py), [R fitting code](joint_fit.R), [parameters](joint_fit/parameters.txt), [diagnostics](joint_fit/diagnostics.txt), [optimizer starts](joint_fit/multistart.txt), [moment residuals](joint_fit/moment_residuals.txt), [200 parameter delete estimates](joint_fit/parameter_deletes.txt) and [source receipt](joint_fit/receipt.json) retain the analysis. The physical factor contains no explicit sleep or exhaustion input, but this is a distinct physical-function construct, not a matched exhaustion-deleted frailty score or evidence of a causal relationship.
