import gzip,csv,hashlib,json
from pathlib import Path
p=Path('work/genoa_mqtl/meQTL_summarystat_chr15.txt.gz');h=hashlib.md5()
with p.open('rb') as f:
 for b in iter(lambda:f.read(8388608),b''):h.update(b)
assert h.hexdigest()=='679b3497ae07e5694c41c9e6aaadaeb8'
rows=[];target=[];n=0
with gzip.open(p,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  n+=1
  if r['CpG']=='cg12594244':rows.append(r)
  if r['CpG']=='cg12594244' and (r['rs']=='rs11636381' or r['ps']=='74134883'):target.append(r)
with open('outputs/stage6_genoa_cg12594244.csv','w',newline='') as f:
 if rows:
  w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
for rs in ['rs4077284','rs2165241']:
 print(rs,[x for x in rows if x['rs']==rs])
print('rows scanned',n,'cpg_rows',len(rows),'target',len(target));print('smallest',min(rows,key=lambda x:float(x['p_wald'])) if rows else None)
Path('outputs/stage6_genoa_cg12594244_receipt.json').write_text(json.dumps({'md5':h.hexdigest(),'rows_scanned':n,'cpg_rows':len(rows),'candidate_exact_rows':len(target),'source':'https://doi.org/10.5281/zenodo.7697509','study_cis_window_bp':50000,'target_distance_bp':83871},indent=2)+'\n')
