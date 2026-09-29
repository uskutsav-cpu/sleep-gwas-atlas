# Dense GWAS QC: mdd

This report summarizes source normalization only; it is not an association result.

- Input: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/mdd.harmonized.tsv.gz`
- Input SHA256: `7aa0c5f4903d951ece36b64d6349c45c0b134a26dbe8dccd7ae3b9f77d030f58`
- Receipt SHA256: `7765b0cf99604518a57af8e84b3c393003c3c19f29ca957b40e19d480a391c37`
- Genome build: GRCh37; ancestry: EUR
- Rows read / retained: 1156949 / 1156949
- Duplicate IDs: 0; ambiguous positions: 0
- Source columns: `SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N`
- EAF mapping: FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD; INFO: NO_INFO_FIELD_IN_SOURCE_HEADER.
- Brain6 normalization did not map an EAF field (1156949 rows) or INFO (1156949 rows), so its local MAF/INFO filters were not applied.
- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.
- This dataset-level QC does not establish participant-level cohort overlap.
