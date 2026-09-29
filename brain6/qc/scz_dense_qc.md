# Dense GWAS QC: scz

This report summarizes source normalization only; it is not an association result.

- Input: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/scz.harmonized.tsv.gz`
- Input SHA256: `02a168a4693e2b1f589401f345626c3f4b966d2080356d25c1419c4fb601c6df`
- Receipt SHA256: `778f1b31a7893fde1f1e975b92c092e6a49059ec18e68a74d5686ec434c0bfba`
- Genome build: GRCh37; ancestry: EUR
- Rows read / retained: 6341702 / 6341702
- Duplicate IDs: 0; ambiguous positions: 0
- Source columns: `SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N`
- EAF mapping: FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD; INFO: NO_INFO_FIELD_IN_SOURCE_HEADER.
- Brain6 normalization did not map an EAF field (6341702 rows) or INFO (6341702 rows), so its local MAF/INFO filters were not applied.
- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.
- This dataset-level QC does not establish participant-level cohort overlap.
