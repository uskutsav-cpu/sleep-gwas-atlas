# Post-failure LD repair sensitivity

**The larger reference does not resolve the ADHD LD–Z mismatch.** This is an exploratory, post-failure sensitivity analysis of the GCTB block, not a validation of the unmodified matrix or the biological interpretation. After the original matrix failed the PSD check, the [protocol](REPAIR_SENSITIVITY_PROTOCOL.md) specified two alternatives before examining these results: clipping negative eigenvalues and renormalizing, and shrinking correlations by 0.001 toward identity. Both adjusted matrices pass numerical symmetry, diagonal and PSD checks. Eigenvalue clipping changed an entry by at most 0.000057 (relative Frobenius change 0.0000069); the identity shrinkage had a maximum change of 0.001 and relative Frobenius change 0.000995. Neither establishes that the resulting correlations equal study-specific LD. The [receipt](repair_sensitivity_receipt.json) records the transformations.

On identical support of 1,547 variants, the table gives the largest absolute standardized LD–Z residual, followed by the number exceeding 5. These cutoffs are descriptive and post hoc, not calibrated pass/fail rules.

| Association dataset | 503-person 1000 Genomes | GCTB, eigenvalue clip | GCTB, identity shrinkage |
|---|---:|---:|---:|
| Oligodendrocyte LINC03051 eQTL | 6.42 (4) | 5.05 (1) | 5.05 (1) |
| Excitatory-neuron LINC03051 eQTL | 5.04 (1) | 4.47 (0) | 4.46 (0) |
| Insomnia GWAS | 6.61 (2) | 4.76 (0) | 4.76 (0) |
| ADHD GWAS | 11.27 (6) | 14.25 (10) | 14.24 (10) |

Some residual summaries are smaller under GCTB, while the most extreme ADHD residual and its count above 5 are larger. The estimated inconsistency parameter is smaller for all four datasets, but that statistic and the residual flags assess different aspects and do not prove compatibility. The two GCTB adjustments are similar here; they do not independently replicate each other or the disease association. Residuals may reflect LD reference mismatch, sample size, allele handling, association-model assumptions, sampling, or other causes; the experiment does not identify which.

The comparison uses the prior scalar sample-size assumptions of 737 for each QTL, 313,750 for insomnia and 128,214 for ADHD. Their effective-sample-size interpretation remains uncertain, especially because ADHD sample sizes vary across variants. GCTB's UK Biobank reference is not the source cohort for the brain QTL, ADHD or Finnish insomnia data and may overlap UK Biobank discovery samples. This experiment was selected after seeing the original matrix's failure, and its transformations introduce assumptions. No new SuSiE fine-mapping, colocalization or LAVA claim was made.

The complete [comparison](repair_ldz_comparison.tsv), [variant-level values](repair_ldz_details.tsv), [diagnostic code](diagnose_repaired.R), [numerical adjustment code](repair_sensitivity.py), [runtime](repair_runtime.txt) and [verification](repair_verification.json) are available. Independent arithmetic reproduced the residuals and summary counts; it does not validate the model or correct sample sizes.

A subsequent [ADHD sample-size and source audit](ADHD_N_RESULTS.md) found that the ten strongest GCTB residuals have lower audited N than the regional median. Using a lower scalar N barely changed the diagnostics; the post-hoc pattern does not establish its cause.
