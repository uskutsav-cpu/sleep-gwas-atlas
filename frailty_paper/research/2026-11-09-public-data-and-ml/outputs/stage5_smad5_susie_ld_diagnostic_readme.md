# SMAD5 skeletal-muscle splice: LD consistency diagnostic

This is a quality-control analysis, not evidence of shared causal risk. For 2,345 aligned variants, susieR 0.14.2 `estimate_s_rss` estimated lambda = 0.06360 for the SMAD5 skeletal-muscle splice event (assumed N = 702) and 0.01627 for the Mak frailty index (median per-variant N = 402,110). Lambda measures inconsistency between summary Z scores and reference LD under the diagnostic model; a higher value indicates more inconsistency. There is no universal cutoff here, and these values do not identify the cause or alone establish whether any particular fine-mapping result is valid.

The LD reference comprises only 503 1000 Genomes European samples, including Finnish individuals; it is not the original genotype data for either association study. Its centered genotype rank is at most 502, versus 2,345 variants. It has incomplete coverage of the published molecular credible sets and contains missing genotypes mean-imputed per variant. Sampling noise, differing ancestry or study design, allele or variant alignment errors, and other model mismatches could contribute to the diagnostics. The Mak sample size varies by variant and the median is an approximation. The splice sample size is taken from the original molecular study metadata and may not capture every event-specific exclusion.

The stored float32 matrix has a minimum eigenvalue of -1.47e-6. The diagnostic software warned that the matrix was not positive semidefinite and set negative eigenvalues to zero internally. This small negative value may arise from numerical rounding, but the warning and the reference's rank deficiency should be retained when interpreting the results. No SuSiE credible sets or multi-signal colocalization probabilities are reported: the combination of LD uncertainty and incomplete molecular credible-set coverage does not support a reliable shared-signal conclusion. The earlier single-signal analysis cannot resolve this molecular event's two published signals.

A separate preliminary high-LD check found no absolute signed-adjusted Z-score differences above 5 among 13,245 reference variant pairs with |r| > 0.999. This arbitrary screening check is not a formal test and does not rule out more subtle alignment problems. Results: `stage5_smad5_high_ld_z_check.json` and `stage5_smad5_high_ld_z_check_top_pairs.csv`.

## Reproduction

From the project root, use Python with NumPy to reconstruct the temporary double-precision matrix and install the stated R package version and dependencies. The diagnostic script expects `work/stage5_Rlib` if using a private R library; otherwise a normally installed susieR is also discoverable.

```sh
mkdir -p work/stage5_susie_diagnostic
python3 -c "import numpy as np; np.load('outputs/stage5_smad5_muscle_eur_reference_ld.npz')['r'].astype('<f8').tofile('work/stage5_susie_diagnostic/eur_ld_float64.bin')"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript outputs/stage5_smad5_susie_ld_diagnostics.R
```

Inputs are `stage5_smad5_susie_ld_diagnostic_inputs.csv`, in precisely the matrix's variant order, and `stage5_smad5_muscle_eur_reference_ld.npz`. Output is `stage5_smad5_susie_ld_consistency.csv`; the original run log is `stage5_smad5_susie_ld_diagnostics.log`. See `stage5_smad5_muscle_eur_reference_ld_readme.md` and `stage5_smad5_muscle_splice_indel_colocalization_review.md` for LD construction, published signal coverage, and the allele-matching limitations. The script computes one eigendecomposition and reuses it in both diagnostics; neither diagnostic should be interpreted as a causal test.

Methods: [susieR summary-statistics diagnostic](https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html); [coloc SuSiE guidance](https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html).

## Variant-level follow-up

We additionally ran `susieR::kriging_rss` on the same matrix and estimated lambda values. For each variant it compares its sample-size-adjusted Z score with the conditional expectation from the other variants and computes a log likelihood ratio for an allele-flip model. The output contains all 2,345 variants per trait, including the original and adjusted Z scores, conditional mean and variance, standardized residual, and the package's plot-highlighting rule (`logLR > 2` and absolute adjusted Z > 2). Both traits have zero variants meeting that rule. This is not proof of correct allele orientation or adequate LD: for example, the largest absolute standardized residual is 6.89 for the splice event and 9.74 for Mak FI. Large residuals can arise from reference-panel sampling, mismatched LD, or other discrepancies; they do not on their own establish an allele error or biological signal.

Results: `stage5_smad5_SMAD5_muscle_splice_variant_ld_diagnostic.csv` and `stage5_smad5_Mak_FI_variant_ld_diagnostic.csv`; code: `stage5_smad5_susie_variant_diagnostics.R`; run log: `stage5_smad5_susie_variant_diagnostics.log`. After reconstructing the binary matrix as above and running the initial diagnostic, run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript outputs/stage5_smad5_susie_variant_diagnostics.R`. As above, negative eigenvalues were truncated by susieR; this is a reference-LD diagnostic, not a clinical test or a causal inference.
