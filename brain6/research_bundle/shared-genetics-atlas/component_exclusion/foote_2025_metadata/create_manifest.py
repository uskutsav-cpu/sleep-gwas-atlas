import hashlib,pathlib,re,csv
p=pathlib.Path(__file__).resolve().parent
labels=['General frailty','Limited social support','Unhealthy lifestyle','Multimorbidity','Metabolic problems','Poorer cognition','Disability','Pneumonia']
rows=[]
for i,label in enumerate(labels,90624046):
 name=f'GCST{i}.yaml'; f=p/name; s=f.read_text();md5=re.search(r'data_file_md5sum:\s*(\w+)',s); size=re.search(r'sample_size:\s*(\d+)',s)
 assert f'gwas_id: GCST{i}' in s and md5 and size
 rows.append({'accession':f'GCST{i}','reported_trait':label,'metadata_file':name,'metadata_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'source_metadata_url':f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90624001-GCST90625000/GCST{i}/GCST{i}.tsv-meta.yaml','summary_filename':f'GCST{i}.tsv','publisher_summary_md5':md5.group(1),'catalog_sample_size_field':size.group(1),'size_interpretation':'Metadata sample-size field; unique participant count unverified; per-variant n is Nhat' if i<=90624052 else 'reported pneumonia case+control total','downloaded_full_statistics':'NO'})
with (p/'catalog_manifest.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print('verified',len(rows),'metadata files and created catalog manifest')
