import csv,gzip,pathlib,hashlib,json,math
base=pathlib.Path('/Volumes/Extreme SSD/brain6-placo-ffs-rerun-20261027'); inputs=[base/'FFS_1_569131.tsv.gz',base/'FFS_569132_1138262.tsv.gz'];m={x['variant']:x for x in csv.DictReader(open('outputs/shared22_eight_source_evidence_matrix.csv'))};ex=next(x for x in csv.DictReader(open('outputs/integrated_15_candidate_review.csv')) if x['variant']=='rs323509');wanted=set(m)|{'rs323509'};hits={};receipts=[]
for p in inputs:
 h=hashlib.sha256();
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 n=0;status={}
 with gzip.open(p,'rt') as f:
  reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames==['SNP','A1','A2','Z_insomnia','Z_frailty','P_insomnia','P_frailty','P_PLACO_plus','status']
  for x in reader:
   n+=1;status[x['status']]=status.get(x['status'],0)+1
   if x['SNP'] in wanted:hits.setdefault(x['SNP'],[]).append(x)
 receipts.append({'path':str(p),'sha256':h.hexdigest(),'rows':n,'status':status})
print('files',receipts)
out=[]
for v in sorted(wanted):
 rows=hits.get(v,[]); rec={'variant':v,'matching_id_rows':len(rows)}
 if len(rows)==1:
  x=rows[0];a,b=(m[v]['score_counted_allele'],m[v]['other_allele']) if v in m else (ex['counted_allele'],ex['other_allele']);rec.update(x);rec['candidate_allele_pair_match']={a,b}=={x['A1'],x['A2']};rec['candidate_allele_pair']=a+'/'+b
  if rec['candidate_allele_pair_match'] and x['status']=='TESTED':
   rec['Z_insomnia_counted_allele']=float(x['Z_insomnia'])*(1 if x['A1']==a else -1);rec['Z_frailty_counted_allele']=float(x['Z_frailty'])*(1 if x['A1']==a else -1)
 out.append(rec)
fields=list(dict.fromkeys(k for x in out for k in x));p=pathlib.Path('outputs/stage3_existing_placo_ffs_candidate_lookup.csv')
with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
pathlib.Path('outputs/stage3_existing_placo_ffs_receipt.json').write_text(json.dumps({'inputs':receipts,'candidate_hits':len(hits),'n_candidates':len(wanted),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'limitations':'Source and PLACO+ configuration not yet authenticated; rsID and allele pair only, coordinates unavailable in these output shards; do not promote to discovery or rank until provenance verified.'},indent=2)+'\n')
print('candidate ids',len(hits),'missing',sorted(wanted-set(hits)));print('all output p results are provisional until source and method provenance checked')
