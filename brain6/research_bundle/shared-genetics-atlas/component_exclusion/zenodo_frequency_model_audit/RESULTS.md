# Component prefixes: allele frequency and standard-error diagnostics

These are **source-quality checks, not genetic association results**. They use 20 previously selected variants from the beginning of chromosome 1 and incomplete file prefixes, so they cannot validate a whole cohort, full variant file or study model.

## Reference frequency check

For 99 of 100 component–variant combinations, the prefix contained a unique matching variant. Counting the allele named `ALLELE1` in 503 European-labeled 1000 Genomes Phase 3 participants gave median absolute differences from the reported `A1FREQ` of **0.016–0.017** across components; the largest checked difference was **0.0372**. One weight-loss variant was outside its stored prefix and was not assessed. These are differences in selected allele frequencies, not an ancestry estimate or independent confirmation that the reported beta is signed for that allele. Small reference-sample size, reference population composition, correlated variants, source missingness and multi-allelic representation limit interpretation.

The reported frequencies for the same variants across the five components differ by at most **0.000197** in this selected set, despite different published component case/control totals. This could have several explanations, including common underlying samples or frequency calculation conventions; it does not establish identical participants, overlap, or how frequencies were calculated.

## Standard-error scale check

The [depositor's README](https://zenodo.org/records/14011550/files/0_Readme.txt?download=1) reports case and control totals for each component and names Regenie v3.5. For the complete rows in each stored prefix, the median of

`SE² × 2 × A1FREQ × (1 − A1FREQ) × cases × controls / (cases + controls)`

was **1.021–1.069**, with the 5th–95th percentiles spanning **1.008–1.145** across the five components. Under a simplified independent-sample, Hardy–Weinberg logistic score approximation with no covariate, relatedness or imputation complications, this quantity is near one. The result is therefore consistent with that scale, but it is **not proof** that the depositor used a logistic model, that the sample counts apply to every variant, or that the standard errors and p-values are calibrated. The prefix observations are LD-correlated and their percentile ranges are descriptive, not confidence intervals or hypothesis tests.

[Frequency comparisons](frequency_comparisons.tsv), [scale summaries](se_scale_summary.tsv), [script](audit.py), [protocol](PROTOCOL.md) and [receipt](receipt.json) retain inputs, counts and hashes. The reference VCF was accessed remotely, and no participant genotypes were saved. The repository's [Regenie documentation](https://rgcgithub.github.io/regenie/options/) describes standard output conventions, but these processed files still require source-specific confirmation of build throughout, phenotype definitions, model/effect allele, ancestry, covariates, overlap and per-variant sample sizes. The previously documented duplicate/multiallelic issues remain. No inference hold has been lifted, and no sleep–component genetic correlation or component-excluded FFS GWAS was calculated.
