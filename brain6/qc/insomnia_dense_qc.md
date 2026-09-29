# Dense GWAS QC: insomnia

This report summarizes source normalization only; it is not an association result.

- Input: `/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/harmonized/insomnia.harmonized.tsv.gz`
- Input SHA256: `23e06f4a682477064e73551eeab78c292af4d2058fe636ddf328c4bccf560ea4`
- Receipt SHA256: `48077ffa0222131ea7b9035ea5acc1f6ed4a2d27fd4c74571bdc1b799f632534`
- Genome build: GRCh37; ancestry: EUR
- Rows read / retained: 6077635 / 6077635
- Duplicate IDs: 0; ambiguous positions: 0
- Source columns: `SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N`
- EAF mapping: FRQ_PRESENT; EAF_SEMANTICS_NOT_MAPPED_IN_BRAIN6_SOURCE_CARD; INFO: NO_INFO_FIELD_IN_SOURCE_HEADER.
- Brain6 normalization did not map an EAF field (6077635 rows) or INFO (6077635 rows), so its local MAF/INFO filters were not applied.
- Source-level allele-frequency filtering is inherited only where an upstream, source-bound QC record documents it; do not infer it from the unused FRQ column.
- This dataset-level QC does not establish participant-level cohort overlap.
