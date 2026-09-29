# Fried FFS signed-effect convention recheck — 2026-09-26 23:58 UTC

## Evidence

- The article-linked Figshare file is named `Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv` and has columns `CHROM POS A1 A2 BETA P SNP SE OR`.
- Ye et al. report BOLT-LMM single-variant testing and use FUMA in the post-GWAS analysis.
- The official FUMA input guide says its header auto-detection treats `A1` as the effect allele, `A2` as the non-effect allele, and `Beta` as the signed beta field. The FUMA tutorial cautions users to observe this default.
- The official BOLT-LMM manual says its standard `ALLELE1` is the effect allele for reported `BETA` coefficients.
- Independently, the acquired file's `OR` equals `exp(BETA)` row-by-row; this does not establish an independent logistic odds-ratio estimate.

## Interpretation and scope

The best-supported interpretation is that `A1` is the effect allele, `A2` is the other allele, and `BETA` is the per-copy signed effect on the ordinal 0–5 FFS as analyzed with BOLT-LMM. This is a documented-format inference from the exact FUMA-labelled file and both tools' standard conventions; no author-provided per-file sidecar explicitly maps the renamed columns. Preserve that distinction in provenance and manuscript wording. Do not use the `OR` column as a separate effect measure.

This closes the practical signed-allele convention needed for standard summary-statistic formatting, subject to the stated inference. It does not establish exact participant overlap with any sleep GWAS. Do not call Fried results independent replication; keep overlap-dependent RG/replication interpretations qualified or blocked under the existing plan. No harmonization, LDSC, LAVA, or other inferential analysis was run during this check.

## Primary documentation

- Ye et al. (2023), source study and data availability: https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/
- FUMA documentation, automatic allele-header mapping: https://fuma-docs.readthedocs.io/en/latest/snp2gene/prepare_input_files.html
- BOLT-LMM manual, `ALLELE1`/`BETA` convention: https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html
- Article-linked Figshare source: https://figshare.com/s/6683396c68807fe4e729
