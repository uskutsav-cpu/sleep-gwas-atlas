#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parent
sources=[root/'mvp_six_leads/six_leads.tsv',root/'finngen_adhd_raw_lookup/indexed_results.tsv']
def read(x):
 with x.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a=read(sources[0]);b={r['variant']:r for r in read(sources[1])};assert len(a)==6 and len({x['variant'] for x in a})==6
out=[]
for x in a:
 y=b[x['variant']]; assert x['chr']==y['chr'] and x['pos38']==y['position38'] and x['effect_allele']==y['effect_allele']
 r={'variant':x['variant'],'chr':x['chr'],'position38':x['pos38'],'allele':x['effect_allele'],'mvp_insomnia_p':x['mvp_p'],'finngen_adhd_p':y['source_p'],'max_p':'','max_p_times_six':'','status':'not estimable'}
 if x['status']=='MATCHED' and y['status']=='MATCHED':
  p1=float(x['mvp_p']);p2=float(y['source_p']);assert all(math.isfinite(z) and 0<=z<=1 for z in (p1,p2))
  m=max(p1,p2);r.update(max_p=m,max_p_times_six=min(1,6*m),status='matched')
 out.append(r)
with (p/'results.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
(p/'receipt.json').write_text(json.dumps({'input_sha256':{str(x.relative_to(root)):hashlib.sha256(x.read_bytes()).hexdigest() for x in sources},'family':6,'matched':sum(x['status']=='matched' for x in out),'interpretation':'postselected numerical benchmark; no adjustment for prior discovery'},indent=2)+'\n')
print((p/'results.tsv').read_text())
