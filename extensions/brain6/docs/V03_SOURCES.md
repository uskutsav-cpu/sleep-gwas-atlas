# Primary sources checked for v0.3

Checked September 8, 2026. Retrieval can return newer documentation than the
versions previously used in an archived run; do not relabel old runs with new versions.

- PLACO/PLACO+ official repository and notes: https://github.com/RayDebashree/PLACO
  Genome-wide variance/correlation estimation; same effect alleles; extreme-Z
  policy; original PLACO versus PLACO+; pleiotropy is association, not mediation.
- SuSiE-RSS diagnostics: https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html
  Matched signed LD/summary-statistic consistency, no automatic causal claims.
- SuSiE-RSS interface: https://stephenslab.github.io/susieR/reference/susie_rss.html
- coloc.susie: https://chr1swallace.github.io/coloc/reference/coloc.susie.html
  All signal-pair hypotheses retained; SNP H4 probabilities conditional on H4.
- coloc prior sensitivity: https://chr1swallace.github.io/coloc/reference/sensitivity.html
- TwoSampleMR harmonization: https://mrcieu.github.io/TwoSampleMR/reference/harmonise_data.html
- LAVA official README: https://github.com/josefin-werme/LAVA
- GenomicSEM official repository: https://github.com/GenomicSEM/GenomicSEM
- TwoSampleMR official repository: https://github.com/MRCIEU/TwoSampleMR

Exact checked GitHub package commits are in configs/native_source_pins.json.
The archive-derived effect/N conventions are read from the user's
config/gwas_schemas.tsv and source registry, not imputed from these method docs.

The exact competitive pathway calculation is an elementary conditional
randomization test: X_s ~ Hypergeometric(N_s,M_s,K_s) and the independent-strata
sum is computed by convolution. Tests compare it against exhaustive enumeration.
It is explicitly not labeled as a third-party software method or a validated
biological mechanism. Within-stratum exchangeability is a substantive assumption.
