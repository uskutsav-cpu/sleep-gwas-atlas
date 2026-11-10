import gzip,csv,hashlib,json,math
from pathlib import Path
p=Path('work/genoa_mqtl/meQTL_summarystat_chr15.txt.gz'); h=hashlib.md5()
with p.open('rb') as f:
 for b in iter(lambda:f.read(8388608),b''):h.update(b)
assert h.hexdigest()=='679b3497ae07e5694c41c9e6aaadaeb8'
rows=[];n=0;by_id=0
with gzip.open(p,'rt') as f:
 for x in csv.DictReader(f,delimiter='\t'):
  n+=1
  if x['rs']=='rs11636381':by_id+=1
  if x['chr']=='15' and x['ps']=='74134883' and set((x['allele1'],x['allele0']))=={'C','T'}:rows.append(x)
rows.sort(key=lambda x:float(x['p_wald']))
with open('outputs/stage6_genoa_rs11636381_mqtl.csv','w',newline='') as f:
 if rows:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
rc={'source':'https://doi.org/10.5281/zenodo.7697509','sha_md5':h.hexdigest(),'rows_scanned':n,'rsid_matches_all_alleles':by_id,'position_allele_pair_matches':len(rows),'window_bp':50000,'nominal_p_lt_0_05':sum(float(r['p_wald'])<.05 for r in rows),'bonferroni_p_lt_0_05':sum(float(r['p_wald'])*len(rows)<.05 for r in rows),'top':rows[:5]}
Path('outputs/stage6_genoa_rs11636381_mqtl_receipt.json').write_text(json.dumps(rc,indent=2)+'\n');print(json.dumps(rc,indent=2))
