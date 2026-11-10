import csv, json
from pathlib import Path
import numpy as np
P=Path('outputs'); direct={r['marker'] for r in csv.DictReader((P/'stage6_framingham_omni5_nexus_coordinate_matches.csv').open())}
signed={(r['score'],r['target']):r for r in csv.DictReader((P/'stage6_framingham_proxy_reference_direction.csv').open())}
models=[]
for score,path in [('score15','unvalidated_15_variant_score_prototype.csv'),('score22','unvalidated_22_shared_direction_score.csv')]:
 for r in csv.DictReader((P/path).open()):
  marker=r['variant'];chr=r.get('chr') or r['chr_GRCh37_assumed'];pos=int(r.get('position_GRCh37_inferred') or r.get('position_GRCh37_assumed'))
  a={'score':score,'marker':marker,'chr':chr,'pos':pos,'counted':r['counted_allele'],'other':r['other_allele'],'direct':marker in direct}
  if score=='score15':a['weights']={'FI':float(r['FI_beta_weight']),'FFS':float(r['FFS_beta_weight_inferred']),'insomnia':float(r['insomnia_beta_weight'])}
  else:a['weights']={'count':1.0}
  if not a['direct']:
   s=signed[(score,marker)];a.update(tag=s['reference_tag'],tag_a1=s['tag_reference_A1'],tag_a2=s['tag_reference_A2'],r=float(s['signed_r_target_counted_vs_tag_A1']))
  models.append(a)
base=Path('work/stage4_magma_reference/g1000_eur'); found={}
target_idx={};tag_idx={}
for x in models:
 target_idx.setdefault((x['chr'],x['pos']),[]).append(x)
 if not x['direct']:tag_idx.setdefault((x['chr'],x['tag']),[]).append(x)
with open(str(base)+'.bim') as f:
 for idx,line in enumerate(f):
  c,rs,_,pos,a1,a2=line.split();pos=int(pos)
  for x in target_idx.get((c,pos),[]):
   if {a1,a2}=={x['counted'],x['other']}:found.setdefault((x['score'],x['marker'],'target'),[]).append((idx,a1,a2))
  for x in tag_idx.get((c,rs),[]):
   if {a1,a2}=={x['tag_a1'],x['tag_a2']}:found.setdefault((x['score'],x['marker'],'tag'),[]).append((idx,a1,a2))
 nvar=idx+1
n=sum(1 for _ in open(str(base)+'.fam'));step=(n+3)//4;bed=np.memmap(str(base)+'.bed',dtype=np.uint8,mode='r');assert bytes(bed[:3])==b'\x6c\x1b\x01' and len(bed)==3+nvar*step
lut=np.array([2,np.nan,1,0],dtype=float)
def decode(i):
 b=bed[3+i*step:3+(i+1)*step];return lut[((b[:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(-1)[:n]]
for x in models:
 a=found.get((x['score'],x['marker'],'target'),[]);assert len(a)==1,(x['marker'],a)
 i,a1,a2=a[0];v=decode(i);x['target']=v if a1==x['counted'] else 2-v
 if x['direct']:x['proxy']=x['target'].copy()
 else:
  b=found.get((x['score'],x['marker'],'tag'),[]);assert len(b)==1,(x['marker'],b)
  i,a1,a2=b[0]; v=decode(i);v=v if a1==x['tag_a1'] else 2-v;x['proxy']=v if x['r']>0 else 2-v
  ok=np.isfinite(x['target']) & np.isfinite(x['proxy']);r=np.corrcoef(x['target'][ok],x['proxy'][ok])[0,1];assert abs(r-abs(x['r']))<2e-6,(x['marker'],r,x['r'])
fam=[line.split()[1] for line in open(str(base)+'.fam')]
panel={r['sample']:r['pop'] for r in csv.DictReader(open('work/1kg_phase3_panel.tsv'),delimiter='\t')}
assert all(x in panel for x in fam)
groups={k:np.array([i for i,sample in enumerate(fam) if panel[sample]==k]) for k in ('CEU','FIN','GBR','IBS','TSI')}
assert sum(len(v) for v in groups.values())==503
results=[]
for score,trait in [('score15','FI'),('score15','FFS'),('score15','insomnia'),('score22','count')]:
 group=[x for x in models if x['score']==score];full=np.sum([x['weights'][trait]*x['target'] for x in group],axis=0);proxy=np.sum([x['weights'][trait]*x['proxy'] for x in group],axis=0);partial=np.sum([x['weights'][trait]*x['target'] for x in group if x['direct']],axis=0)
 for name,indices in groups.items():
  ok=np.isfinite(full[indices])&np.isfinite(proxy[indices])&np.isfinite(partial[indices]);ix=indices[ok];assert len(ix)>50
  results.append({'score':score,'weight':trait,'reference_subgroup':name,'n_complete':len(ix),'proxy_full_correlation':round(float(np.corrcoef(full[ix],proxy[ix])[0,1]),6),'direct_only_full_correlation':round(float(np.corrcoef(full[ix],partial[ix])[0,1]),6),'proxy_reference_sd_ratio':round(float(np.std(proxy[ix],ddof=1)/np.std(full[ix],ddof=1)),6)})
with (P/'stage6_framingham_proxy_score_eur_subgroups.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
print(json.dumps(results,indent=2))
