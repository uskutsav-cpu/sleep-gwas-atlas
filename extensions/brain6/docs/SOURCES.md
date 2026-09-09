# Primary implementation references

Consulted September 8, 2026. These document method interfaces and assumptions; they do **not** certify this new adapter code or substitute for production execution tests.

| Component | Primary documentation/source |
|---|---|
| LDSC original CLI | https://github.com/bulik/ldsc |
| PLACO / PLACO+ | https://github.com/RayDebashree/PLACO |
| PLACO v0.2.0 API | https://raw.githubusercontent.com/RayDebashree/PLACO/master/PLACO_v0.2.0.R |
| PLACO+ paper | Park and Ray (2025), Human Genetics and Genomics Advances, doi:10.1016/j.xhgg.2025.100501 |
| Original PLACO paper | Ray and Chatterjee (2020), PLOS Genetics, doi:10.1371/journal.pgen.1009218 |
| LAVA API | https://github.com/josefin-werme/LAVA |
| PLINK 1.9 clumping | https://www.cog-genomics.org/plink/1.9/postproc#clump |
| PLINK signed LD | https://www.cog-genomics.org/plink/1.9/ld |
| SuSiE summary-statistics interface | https://stephenslab.github.io/susieR/reference/susie_rss.html |
| SuSiE-RSS paper | Zou et al. (2022), PLOS Genetics, doi:10.1371/journal.pgen.1010299 |
| coloc input data structures | https://chr1swallace.github.io/coloc/articles/a02_data.html |
| coloc multi-signal analysis | https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html |
| coloc runsusie | https://chr1swallace.github.io/coloc/reference/runsusie.html |
| GenomicSEM multivariate LDSC | https://github.com/GenomicSEM/GenomicSEM/blob/master/man/ldsc.Rd |
| GenomicSEM usermodel | https://github.com/GenomicSEM/GenomicSEM/blob/master/man/usermodel.Rd |
| TwoSampleMR harmonization | https://mrcieu.github.io/TwoSampleMR/articles/harmonise.html |
| TwoSampleMR methods | https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html |
| Official pleioFDR | https://github.com/precimed/pleiofdr |
| Official MATLAB runme | https://raw.githubusercontent.com/precimed/pleiofdr/master/runme.m |
| Official MATLAB configuration | https://raw.githubusercontent.com/precimed/pleiofdr/master/config_default.txt |
| Snakemake rules/resources | https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html |

## Source-specific decisions retained in the implementation

PLACO+ fits genome-wide variance/correlation parameters, handles potentially correlated traits, and recommends explicit filtering of extreme Z scores. Its developers caution against candidate-only testing and interpreting pleiotropic association as a causal mechanism. The native adapter calls their R functions rather than rewriting tail integrals.

SuSiE/coloc use signed, order-matched LD and report credible sets/signals under the specified model. The adapter does not automatically estimate residual variance as if external LD were in-sample LD. Both GWAS effects are aligned to the alleles counted in the LD reference.

pleioFDR is run through the official MATLAB program and its configuration structure. Its method-specific random-pruning settings are not equivalent to arbitrarily thinning a SuSiE fine-mapping locus. Source-specific settings need review rather than a universal no-pruning switch.

No upstream implementation, genotype reference, controlled dataset, MATLAB license or R package binary is redistributed in this delivery.


## v0.2 additions (checked September 8, 2026)

- TwoSampleMR native workflow and estimator output schema: https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html
- MR-specific LD-clumping defaults: https://mrcieu.github.io/TwoSampleMR/reference/clump_data.html
- Full QTL summary-data access and effect-allele conventions: https://www.ebi.ac.uk/eqtl/Data_access/
- SuSiE summary/LD mismatch diagnostics: https://stephenslab.github.io/susieR/articles/susierss_diagnostic.html
- Multiple-signal colocalization: https://chr1swallace.github.io/coloc/articles/a06_SuSiE.html
- PLACO/PLACO+ genome-wide parameter estimation and restrictions: https://github.com/RayDebashree/PLACO

The external pages document methods; they are not evidence that this user's analysis has executed. No upstream R implementation is redistributed here. Pin and review native source files/packages on the actual data host before use.
