# Dense GWAS QC: adhd

This report summarizes source normalization only; it is not an association result.

- Input: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/adhd.harmonized.tsv.gz`
- Input SHA256: `25cdc39268559321465583a445fadc02b49efc1019cfa9a3a07a33e2b24dfcd3`
- Receipt SHA256: `50f2be56bfec0c597953c7f850b45f54ae08b37f16626a84195d19875cb7cea8`
- Genome build: GRCh37; ancestry: EUR
- Rows read / retained: 5692669 / 5692669
- Duplicate IDs: 0; ambiguous positions: 0
- Source columns: `SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N`
- EAF mapping: FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD; INFO: NO_INFO_FIELD_IN_SOURCE_HEADER.
- Brain6 normalization did not map an EAF field (5692669 rows) or INFO (5692669 rows), so its local MAF/INFO filters were not applied.
- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.
- This dataset-level QC does not establish participant-level cohort overlap.
