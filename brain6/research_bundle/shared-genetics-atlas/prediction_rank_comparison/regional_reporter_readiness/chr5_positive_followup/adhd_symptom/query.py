import hashlib,json,math,pathlib,pysam
p=pathlib.Path(__file__).resolve().parent;root=p.parents[3];url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90568001-GCST90569000/GCST90568440/harmonised/GCST90568440.h.tsv.gz';fields=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','rsid','hm_coordinate_conversion','hm_code','variant_id']
raw=list(pysam.TabixFile(url).fetch('5',104568697,104568698));rows=[]
for line in raw:
 v=line.split('\t');assert len(v)==len(fields);rows.append(dict(zip(fields,v)))
matched=[r for r in rows if r['chromosome']=='5' and r['base_pair_location']=='104568698' and r['rsid']=='rs12659965' and {r['effect_allele'],r['other_allele']}=={'G','A'}]
status='NO_MATCH' if not matched else ('AMBIGUOUS' if len(matched)>1 else 'MATCHED')
if status=='MATCHED':
 x=matched[0];b,se,pv,fr=map(float,(x['beta'],x['standard_error'],x['p_value'],x['effect_allele_frequency']));assert all(map(math.isfinite,(b,se,pv,fr))) and se>0 and 0<=pv<=1 and 0<=fr<=1
metadata=root/'replication/childhood_symptoms/metadata.yaml';receipt=dict(url=url,query='5:104568698 GRCh38',status=status,raw_rows=rows,matched=matched,raw_sha256=hashlib.sha256(('\n'.join(raw)+'\n').encode()).hexdigest(),metadata_sha256=hashlib.sha256(metadata.read_bytes()).hexdigest(),whole_file_or_index_hashed=False)
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(status,'raw rows',len(rows),'exact matches',len(matched))
