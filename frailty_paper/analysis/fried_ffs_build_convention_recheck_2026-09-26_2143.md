# Fried Frailty Score source-convention recheck — 2026-09-26 21:43 UTC

## Question

Does newly checked source metadata resolve the build and signed-effect convention for the exact Ye et al. (2023) Figshare file acquired at `frailty_paper/data/gwas/fried_frailty_score/Frailty_GWAS_maf0.01_info0.9_hwe1e-6_FUMA.tsv`?

## Evidence checked

- The primary Ye et al. paper identifies the UK Biobank Fried Frailty Score GWAS and links its full summary statistics to the Figshare share page. It reports use of FUMA for annotation, but the retrieved article text does not define A1/BETA direction or state a build for the exact file. [Ye et al. (2023)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10651618/)
- The official PGS Catalog record for PGS005229 identifies its source GWAS as `GCST90295968` / PMID 36928559, with 386,565 European-ancestry UK Biobank participants, and records the original genome build as GRCh37. This corroborates the study-level build for the same accession, but does not explicitly identify the exact article-linked Figshare file as the file from which the accession metadata were derived. [PGS Catalog PGS005229](https://www.pgscatalog.org/score/PGS005229/)
- The current GWAS Catalog page for `GCST90295968` says full summary statistics are not available on the Catalog page. It therefore does not supply an accession-side file or metadata record that can be compared byte-for-byte with the Figshare object. [GWAS Catalog GCST90295968](https://www.ebi.ac.uk/gwas/studies/GCST90295968)
- FUMA's official input documentation says automatic header recognition treats `A1` as the effect allele and `A2` as the non-effect allele by default, and recognizes `Beta` as the beta column. The acquired filename contains `_FUMA`, and its header is `CHROM POS A1 A2 BETA P SNP SE OR`; this is consistent with that convention. However, the filename and column names do not independently prove that the authors supplied this exact file to FUMA with those labels or that no allele transformation occurred before upload. [FUMA input documentation](https://fuma-docs.readthedocs.io/en/latest/snp2gene/prepare_input_files.html)
- A 2024 cross-trait paper also reports GRCh37 for the Ye et al. frailty summary statistics, but does not identify the acquired Figshare filename. [El-Saadi et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11069310/)

## Decision

The study-level evidence now supports GRCh37 more strongly: both the official PGS Catalog entry for the same GWAS accession and a published downstream study report GRCh37. The Figshare file is directly linked by the primary article and matches the reported 8,883,488-variant count, but the available Catalog page has no full-statistics file and no exact-object crosswalk was found. Retain the registered exact-file build as **UNKNOWN** until that link is authoritative or a documented reference-coordinate concordance audit is completed.

The FUMA default convention is consistent with interpreting `A1` as effect allele and `A2` as non-effect allele. It is not a source-specific signed-effect statement for this exact file. Retain effect-allele semantics as unresolved; do not flip, harmonize, or analyze the file based only on header conventions or `OR == exp(BETA)`.

Acquisition integrity remains valid. No analysis eligibility or locked threshold changed. Exact participant intersections with the UKB-based sleep sources remain unknown.
