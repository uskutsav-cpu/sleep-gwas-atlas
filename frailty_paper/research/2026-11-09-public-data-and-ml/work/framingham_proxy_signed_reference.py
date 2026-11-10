import csv,json
from pathlib import Path
import numpy as np
P=Path('outputs')
tags={r['variant']:r for r in csv.DictReader((P/'stage6_framingham_omni5_reference_tags.csv').open())}
direct={r['marker'] for r in csv.DictReader((P/'stage6_framingham_omni5_nexus_coordinate_matches.csv').open())}
scores=[]
for label,path in [('score15','unvalidated_15_variant_score_prototype.csv'),('score22','unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader((P/path).open()):
  if r['variant'] in direct:continue
  t=tags[r['variant']];scores.append({'score':label,'variant':r['variant'],'chr':t['chr_GRCh37'],'pos':int(t['position_GRCh37']),'counted':r['counted_allele'],'other':r['other_allele'],'tag':t['best_tag_id'],'tag_pos':int(t['best_tag_pos']),'expected_r2':float(t['best_r2'])})
base=Path('work/stage4_magma_reference/g1000_eur'); found={}
with open(str(base)+'.bim') as f:
 for idx,line in enumerate(f):
  c,rs,_,p,a1,a2=line.split();p=int(p)
  for x in scores:
   key=None
   if c==x['chr'] and p==x['pos'] and {a1,a2}=={x['counted'],x['other']}:key=(x['variant'],'target')
   if c==x['chr'] and p==x['tag_pos'] and rs==x['tag']:key=(x['variant'],'tag')
   if key and (idx,a1,a2) not in found.setdefault(key,[]):found[key].append((idx,a1,a2))
 nvar=idx+1
n=sum(1 for _ in open(str(base)+'.fam'));step=(n+3)//4;bed=np.memmap(str(base)+'.bed',mode='r',dtype=np.uint8);assert len(bed)==3+nvar*step and bytes(bed[:3])==b'\x6c\x1b\x01'
lut=np.array([2,np.nan,1,0],dtype=float)
def decode(i):
 b=bed[3+i*step:3+(i+1)*step];return lut[((b[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:n]]
result=[]
for x in scores:
 a=found.get((x['variant'],'target'),[]);b=found.get((x['variant'],'tag'),[])
 assert len(a)==1 and len(b)==1,(x['variant'],a,b)
 ai,aa,ab=a[0];bi,ba,bb=b[0]; y=decode(ai);z=decode(bi)
 if aa!=x['counted']:y=2-y
 ok=np.isfinite(y)&np.isfinite(z);r=float(np.corrcoef(y[ok],z[ok])[0,1]);assert abs(r*r-x['expected_r2'])<2e-6,(x['variant'],r*r,x['expected_r2'])
 result.append({'score':x['score'],'target':x['variant'],'target_counted_allele':x['counted'],'reference_tag':x['tag'],'tag_reference_A1':ba,'tag_reference_A2':bb,'signed_r_target_counted_vs_tag_A1':round(r,6),'r_squared':round(r*r,6),'paired_n':int(ok.sum())})
with (P/'stage6_framingham_proxy_reference_direction.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=result[0].keys());w.writeheader();w.writerows(result)
print(json.dumps({'rows':len(result),'positive':sum(r['signed_r_target_counted_vs_tag_A1']>0 for r in result),'negative':sum(r['signed_r_target_counted_vs_tag_A1']<0 for r in result),'n_reference':n}))
