# Dense GWAS QC: longsleep

This report summarizes source normalization only; it is not an association result.

- Input: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/longsleep.harmonized.tsv.gz`
- Input SHA256: `ff0f50390411fa0d0b236f5b38311ee5ec0df423707adfd5ab18b42ba8a59d76`
- Receipt SHA256: `995df4b56b076075cecf7ece8f1fee3436ada7e2c14e686a50fc4a22e5038eb7`
- Genome build: GRCh37; ancestry: EUR
- Rows read / retained: 6549769 / 6549769
- Duplicate IDs: 0; ambiguous positions: 0
- Source columns: `SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N`
- EAF mapping: FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD; INFO: NO_INFO_FIELD_IN_SOURCE_HEADER.
- Brain6 normalization did not map an EAF field (6549769 rows) or INFO (6549769 rows), so its local MAF/INFO filters were not applied.
- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.
- This dataset-level QC does not establish participant-level cohort overlap.
