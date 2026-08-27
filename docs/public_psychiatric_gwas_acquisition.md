# Public psychiatric GWAS acquisition

This log records the exact public psychiatric releases selected for
`atlas-v1.0`. Source verification means identity, byte count, checksum,
container traversal, literal schema, phenotype, sample, ancestry, and build
claims were checked. It does not mean the ignored raw file is present or that
harmonization, prevalence, LDSC, or downstream analysis has passed.

| Trait | Selected release | Verified file evidence | Remaining blocker |
|---|---|---|---|
| ADHD | [Demontis 2023 PGC Figshare v1](https://doi.org/10.6084/m9.figshare.22564390.v1) | `ADHD2022_iPSYCH_deCODE_PGC.meta.gz`; 200,588,408 bytes; MD5 `ea57120c720e880e34f67f57483dcbb4`; SHA-256 `c58a96031ec44b1edba81c91d603037a2efd03fcd946531424866e3f5f1d40c5`; 6,774,224 rows; hg19; 38,691 cases and 186,843 controls. | Raw archive is not materialized. PGC/iPSYCH terms limit use to scientific or educational non-commercial research and prohibit reposting. |
| Bipolar disorder | [Mullins 2021 PGC Figshare v2](https://doi.org/10.6084/m9.figshare.14102594.v2) | `pgc-bip2021-all.vcf.tsv.gz`; 385,014,083 bytes; MD5 `02e610aaf630e0c869a22fe179d37067`; SHA-256 `3bf42353fd8fcd832986c3ad122747fba458e9712556659fe73b16a91b3081dd`; 7,608,183 rows; GRCh37; 41,917 cases and 371,549 controls. | Raw archive is not materialized. Canonicalization must strip `##` metadata and convert `NEFFDIV2` to total effective N; the harmonizer has a tested named-field path. |
| Schizophrenia | [Trubetskoy 2022 PGC Figshare v7](https://doi.org/10.6084/m9.figshare.19426775.v7) | European autosomal v3 file; 239,710,564 bytes; MD5 `6ebe2376f5cda972d37efa0f214c4df0`; SHA-256 `dbba3a85575c99fd1c2e3497d0c7a44539ccfdbf2e69742f8bdbda879230bcd7`; 7,659,767 rows; GRCh37; 53,386 cases and 77,258 controls including the documented trio convention. | Raw archive is not materialized. Multi-ancestry core/primary files remain secondary-analysis candidates; the locked EUR LDSC input must not silently change to them. |
| Major depression | [Howard 2019 Edinburgh DataShare](https://doi.org/10.7488/ds/2458) | `PGC_UKB_depression_genome-wide.txt`; 364,301,885 bytes; repository MD5 `ed4597a4e7fa168fb96970e3286a0b31`; SHA-256 `082a1602e405fa5e8b40917521f063444131db68e2f435e7aacfd80084852110`; exactly 8,483,301 seven-field rows; 170,756 cases and 329,443 controls; 23andMe excluded. | The source has rsIDs but no chromosome/position. An audited rsID-to-GRCh37 mapping is required before harmonization. The repository labels the dataset CC BY 4.0. |

The schizophrenia publication headline of 76,755 cases and 243,649 controls is
a two-stage maximum. Its deCODE contribution was evaluated at 1,249 selected
variants, so those counts are not attached to the genome-wide European file.
Likewise, the Howard paper's 246,363/561,190 headline includes 23andMe; the
selected public genome-wide no-23andMe release is the documented
170,756/329,443 UK Biobank plus PGC_139k meta-analysis.
