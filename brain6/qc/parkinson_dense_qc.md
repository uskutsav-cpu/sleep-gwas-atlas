# Dense GWAS QC: parkinson

This report summarizes source normalization only; it is not an association result.

- Input: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/parkinson.harmonized.tsv.gz`
- Input SHA256: `d6fef80964c812c66cadbaa99a38c040dac22383f84f81de5558edada067b132`
- Receipt SHA256: `b7c95dd9ed7e6e2f457cbd3e56cd3cda443a9a6711c13fd6328a94fa05a87f71`
- Genome build: GRCh37; ancestry: EUR
- Rows read / retained: 1132078 / 1132078
- Duplicate IDs: 0; ambiguous positions: 0
- Source columns: `SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N`
- EAF mapping: FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD; INFO: NO_INFO_FIELD_IN_SOURCE_HEADER.
- Brain6 normalization did not map an EAF field (1132078 rows) or INFO (1132078 rows), so its local MAF/INFO filters were not applied.
- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.
- This dataset-level QC does not establish participant-level cohort overlap.
