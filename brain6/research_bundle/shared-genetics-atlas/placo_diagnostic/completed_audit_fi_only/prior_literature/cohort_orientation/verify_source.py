import csv, gzip, json, pathlib
p=pathlib.Path(__file__).resolve().parent
source=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/v9-public-data/atkins-2021/Atkins_2021_FrailtyIndex_GWAS_EuropeansAged60to70years_GRCh37_UKBiobank_TwinGene_associations.txt.gz')
with (p/'cohort_rows.tsv').open() as f:
    saved={r['SNP']:r for r in csv.DictReader(f,delimiter='\t')}
counts={k:0 for k in saved}; rows=0
with gzip.open(source,'rt') as f:
    reader=csv.DictReader(f,delimiter='\t')
    for row in reader:
        rows+=1
        snp=row['SNP']
        if snp in saved:
            counts[snp]+=1
            assert row==saved[snp], (snp,'saved row differs from source')
assert all(v==1 for v in counts.values()),counts
(p/'source_verify.json').write_text(json.dumps({'source':str(source),'source_rows':rows,'selected_snp_occurrences':counts,'saved_rows_exactly_match':True,'gzip_stream_completed':True},indent=2)+'\n')
print('Verified',rows,'source rows; each target appears exactly once and matches the extract.')
