# Dense GWAS QC: bipolar

This report summarizes source normalization only; it is not an association result.

- Input: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/bipolar.harmonized.tsv.gz`
- Input SHA256: `f55dc16a1c45f90189e7469065409e92ac68b0c47eb34de4a307a81fc518e15a`
- Receipt SHA256: `d5584d5b8cfbfab34b198e8ee39355c7ca8a3ae86232aa81e1eed5a02e62a235`
- Genome build: GRCh37; ancestry: EUR
- Rows read / retained: 6339466 / 6339466
- Duplicate IDs: 0; ambiguous positions: 0
- Source columns: `SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N`
- EAF mapping: FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD; INFO: NO_INFO_FIELD_IN_SOURCE_HEADER.
- Brain6 normalization did not map an EAF field (6339466 rows) or INFO (6339466 rows), so its local MAF/INFO filters were not applied.
- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.
- This dataset-level QC does not establish participant-level cohort overlap.
