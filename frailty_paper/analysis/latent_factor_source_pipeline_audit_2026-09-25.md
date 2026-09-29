# Latent-factor Q_SNP source-pipeline audit — 2026-09-25

## Author-pipeline details verified

The primary [Foote et al. paper](https://www.nature.com/articles/s41588-025-02269-0) links its custom analysis code; the corresponding [Zenodo code release](https://doi.org/10.5281/zenodo.15654249) archives a versioned repository snapshot. The publicly accessible author-repository `main` versions of [`6_Calculate_Q_and_Nhat.R`](https://github.com/IsyFoote/Frailty-Multivariate-GWAS/blob/main/6_Calculate_Q_and_Nhat.R) and [`12_LDSC_Frailty_latent_factors.R`](https://github.com/IsyFoote/Frailty-Multivariate-GWAS/blob/main/12_LDSC_Frailty_latent_factors.R) document the following pipeline:

1. Remove rows with missing `Z_Estimate`, nonzero `warning`, or nonzero `Z_smooth`.
2. Define Q_SNP hits with the inclusive rule `Q_SNP_pval <= 5e-8/7`.
3. Expand each Q_SNP hit by 1 Mb on both sides and remove the resulting SNPs from the post-GWAS set.
4. Estimate effective sample size from the retained rows with MAF 0.10–0.40, then write the `noQ_postGWAS` file.

Script 12 consumes those seven `noQ_postGWAS` files, munges against HapMap3 with INFO 0.9 and MAF 0.01, and supplies the factor-specific sample-size vector to LDSC. The Zenodo record lists both scripts in its archive, but the web-accessible archive did not expose their contents for a byte-level comparison with GitHub `main`; exact archive/source-version parity remains unverified.

## Comparison with this project’s source exports

The downloaded GWAS Catalog exports expose coordinates and `Q_metric_pvalue`, so the local audit can test whether source Q_SNP regions would remove any Catalog-export rows. They do not include the author intermediate fields `Z_Estimate`, `warning`, or `Z_smooth`; consequently this check cannot reproduce all prefilters or certify byte-for-byte equivalence to the author-generated `noQ_postGWAS` inputs. The audit remains a bounded source-Q-region check, not a full reproduction of the latent-factor GWAS preparation pipeline.

The boundary operator has been aligned to the author script (`<=`) and covered by a regression test. A streaming equality check over all 32,035,589 rows found zero values below, exactly equal to, or below-or-equal to `5e-8/7` across the seven factors. Thus the former strict comparison and corrected inclusive comparison produce the same Q_SNP hit set for these exports. The full coordinate-window audit also finds zero significant Q_SNP variants and zero rows excluded; no analysis input, estimate, or result changes.

## Local verification

The corrected audit was rerun with `make -C frailty_paper PYTHON=<pinned Python 3.11.11> audit-latent-qsnp`. Output is written with explicit LF line endings for stable TSV bytes and clean diff checks. It completed with `source_rows=32035589`, `significant_q_variants=0`, and `excluded_rows=0`. The refreshed table SHA-256 is `762eb1e273e4f13f0c994de2cf6e809194fc11a9a126965789892992fd5c6c62`; provenance SHA-256 is `636c9905cb9a2f98c24f7a0f9fd139567b7a5b660d9600bad2f28b4cfedc7abb`; audit-script SHA-256 is `7955fe935e85d0b69684b91f393e11391076e7b47aed612dfd77b057974d37c7`. The updated boundary regression test passes, and the pinned package suite passes 131/131 tests with the untracked Brain6-only test module excluded because that interpreter lacks `pytest`; captured suite log: `latent_qsnp_audit_unit_suite_2026-09-25.log` (SHA-256 `e71447eb2fc72008615785f027f3987e696734cfda976fafaa22cde7a903be8d`). The full integrated `preflight` has not yet been rerun after this narrowly scoped threshold correction.

## Interpretation and remaining requirement

The seven-factor results remain a secondary sensitivity family. The zero-window finding only means that the author’s Q_SNP-window removal would remove no rows from these particular Catalog exports. It does not resolve exact cohort overlap, establish that these exports are identical to the source `noQ_postGWAS` outputs, or establish full upstream Genomic SEM reproducibility. Preserve those limits and do not promote the latent results to independent replication.
