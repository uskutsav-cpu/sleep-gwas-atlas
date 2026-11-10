import csv,gzip,hashlib,json,pathlib,collections
scores=pathlib.Path('work/PGS005228_complete.txt.gz'); source=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz')
expected={'score_md5':'584c4839ba4d8caa299ff558ed161259','fi_sha256':'241ed294ccc367b22787a0c051f8d377c020e16089e82920a555286f3148e6de'}
def hashfile(path,algo):
 h=hashlib.new(algo)
 with path.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
assert hashfile(scores,'md5')==expected['score_md5'];assert hashfile(source,'sha256')==expected['fi_sha256']
selected=[];rsids=set()
with gzip.open(scores,'rt',newline='') as f:
 reader=csv.DictReader((s for s in f if not s.startswith('#')),delimiter='\t')
 for x in reader:
  key=x['rsID']+'|'+x['chr_name']+'|'+x['chr_position']
  if int.from_bytes(hashlib.sha256(key.encode()).digest()[:2],'big')<64:
   selected.append(x);rsids.add(x['rsID'])
print('sampled',len(selected),'unique rsid',len(rsids),flush=True)
found=collections.defaultdict(list)
with gzip.open(source,'rt',newline='') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  if x['SNP'] in rsids:found[x['SNP']].append(x)
c=collections.Counter(); rows=[]
for s in selected:
 a=found.get(s['rsID'],[]);c['sampled']+=1
 if len(a)==0:state='not_found_in_filtered_FI'
 elif len(a)>1:state='multiple_rsid_matches_unresolved'
 else:
  y=a[0];coord=s['chr_name']==y['CHR'] and s['chr_position']==y['BP'];effect=s['effect_allele'];
  if not coord:state='rsid_found_coordinate_mismatch'
  elif y['A1']==y['A2']:state='invalid_fi_alleles'
  elif effect==y['A1']:state='exact_coordinate_effect_equals_FI_A1'
  elif effect==y['A2']:state='exact_coordinate_effect_equals_FI_A2'
  else:state='exact_coordinate_effect_not_in_FI_pair'
 c[state]+=1
 rows.append({'PGS_variant':s['rsID'],'chr':s['chr_name'],'pos':s['chr_position'],'pgs_effect_allele':s['effect_allele'],'pgs_other_allele':s['other_allele'],'matched_FI_rows':len(a),'FI_A1':a[0]['A1'] if len(a)==1 else '','FI_A2':a[0]['A2'] if len(a)==1 else '','FI_beta_A1':a[0]['BETA'] if len(a)==1 else '','FI_p':a[0]['P'] if len(a)==1 else '','PGS_weight':s['effect_weight'],'status':state})
assert sum(v for k,v in c.items() if k!='sampled')==len(selected)
with open('outputs/fi_pgs_allele_sample_audit.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
out={'sample_method':'SHA256(rsID|chr|pos) first 16 bits < 64; deterministic approximately 1/1024; score file order independent','counts':dict(c),'sampled_unique_rsid':len(rsids),'score_md5':expected['score_md5'],'FI_input_sha256':expected['fi_sha256'],'limitations':'Matched, filtered FI subset is not the full score; no weight or phenotype validation; all score rows have duplicated other allele; no automated repair.'}
pathlib.Path('outputs/fi_pgs_allele_sample_summary.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
