# Genomic SEM chromosome-split validation

Genomic SEM model discovery used odd chromosomes and confirmation used even
chromosomes. The input set was fixed before either split fit: 42 of 45 traits
passed the multivariable LDSC diagonal gate. T2D, melanoma, and HDL were
excluded for the prespecified h2-Z/intercept criteria. The two independent
GenomicSEM LDSC structures each contain 42 traits, 903 covariance elements, and
947 jackknife blocks.

The split genetic-correlation estimates were strongly but imperfectly
reproducible: off-diagonal Pearson correlation 0.9111, Spearman correlation
0.8129, median absolute difference 0.0444, and sign agreement 80.1%.

Exploratory factor analyses compared one through ten factors with ML and
minimum-residual extraction after explicitly recorded positive-definite
smoothing of the odd-chromosome correlation matrix. Confirmatory candidates
were derived only from the odd-chromosome ML loadings. The tested family
included primary-loading models and models retaining all absolute loadings of
at least 0.30. A six-factor boundary-residual sensitivity model fixed only the
three Heywood residuals detected in odd-chromosome discovery before inspecting
even-chromosome confirmation.

No candidate passed the held-out criteria CFI at least 0.90, SRMR at most 0.10,
and no negative observed-trait residuals. The closest admissible-fit candidate
was the six-factor salient-cross-loading model (CFI 0.9665, SRMR 0.0888), but it
had negative sleep-duration and Crohn residuals. The discovery-boundary version
still developed a new negative Crohn residual on even chromosomes (CFI 0.9678,
SRMR 0.0888). Every candidate also triggered GenomicSEM's warning that matrix
smoothing changed at least one genetic-covariance Z statistic by more than
0.025.

Therefore the latent model is **not validated**. Factor GWAS and Q_SNP are not
run, because doing so would promote a failed model into an apparently final
result. `scripts/30_finalize_genomicsem.py` instead validates the exact
ten-candidate family and publishes schema-bearing, header-only factor-GWAS and
Q_SNP tables plus `factor_gwas.provenance.json` with terminal status
`NOT_APPLICABLE_NO_VALIDATED_MODEL`. The real model comparison and all loadings
remain available in the ignored local result tables for audit:

- `results/tables/genomic_sem_efa_models.tsv`
- `results/tables/genomic_sem_model_fit.tsv`
- `results/tables/genomic_sem_factor_loadings.tsv`
- `results/tables/genomic_sem_model_syntax.tsv`
- `results/tables/genomic_sem_split_diagnostics.tsv`

Rerun with:

```bash
python3 scripts/28_prepare_chromosome_split.py
Rscript scripts/28_chromosome_split_covariance.R --split odd
Rscript scripts/28_chromosome_split_covariance.R \
  --split even \
  --munged-dir data/munged_chromosome_split/even \
  --ld-dir ref/eur_w_ld_chr_even \
  --out results/tables/genomicsem_validation_even.rds \
  --metadata results/tables/genomicsem_validation_even_metadata.tsv \
  --log-prefix results/logs/genomicsem/validation_even
Rscript scripts/29_genomicsem_model.R
python3 scripts/30_finalize_genomicsem.py
```
