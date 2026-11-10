import runpy,collections,json,openpyxl
import numpy as np
m=runpy.run_path('work/matoba_fixed_baseline_transfer.py'); b=m['b']; effects=m['effects'];seq=m['seq'];loc=m['loc'];comp=m['comp'];pairs=[];ids=[]
for rs,p in effects.items():
 if len(loc[rs])!=1 or next(iter(loc[rs]))[0] not in ('1','2'):continue
 if len(seq[rs]['A1'])!=1 or len(seq[rs]['A2'])!=1:continue
 a=next(iter(seq[rs]['A1']));c=next(iter(seq[rs]['A2']))
 if len(a)!=150 or len(c)!=150 or sum(x!=y for x,y in zip(a,c))!=1:continue
 pairs.append((a,c,int(p<.05)));ids.append(rs)
if len(ids)!=108:raise ValueError('unexpected sample')
w=openpyxl.load_workbook('work/matoba2026/table3.xlsx',read_only=True,data_only=True);it=iter(w['a. Regulatory activity (Veh)'].iter_rows(values_only=True));next(it);loci=collections.defaultdict(set)
for r in it:
 if len(r)>1 and r[1] in ids:loci[str(r[1])].add(str(r[0]))
if any(len(loci[x])!=1 for x in ids):raise ValueError('ambiguous locus assignment')
y=np.array([v[2] for v in pairs]);X=[]
for a,c,t in pairs:
 a,c=a.translate(comp)[::-1],c.translate(comp)[::-1];X.append(np.abs(b['kcounts'](a)-b['kcounts'](c)))
X=np.clip((np.stack(X)-b['mean'])/b['sd'],-20,20);scores=1/(1+np.exp(-np.clip(X@b['wgt']+b['intercept'],-30,30)))
clusters=collections.defaultdict(list)
for i,rs in enumerate(ids):clusters[next(iter(loci[rs]))].append(i)
def auc(x,y):
 positives=x[y==1];negatives=x[y==0]
 if not len(positives) or not len(negatives):return np.nan
 return float(((positives[:,None]>negatives).sum()+.5*(positives[:,None]==negatives).sum())/(len(positives)*len(negatives)))
point=auc(scores,y)
if abs(point-m['r']['reference_plus_orientation']['auroc'])>1e-10:raise ValueError('point mismatch')
rng=np.random.default_rng(20261028);keys=list(clusters);draws=[]
for _ in range(5000):
 sampled=rng.choice(keys,len(keys),replace=True);idx=np.array([i for k in sampled for i in clusters[k]],int);v=auc(scores[idx],y[idx]);
 if np.isfinite(v):draws.append(v)
r={'variants':108,'reported_locus_clusters':len(keys),'largest_cluster':max(map(len,clusters.values())),'draws_requested':5000,'draws_valid':len(draws),'reference_plus_auroc':point,'bootstrap_percentile_95':list(map(float,np.quantile(draws,[.025,.975]))),'method':'resample reported locus labels with replacement; fixed model; preserves within-cluster variants; no class stratification','seed':20261028}
open('outputs/stage6_matoba_locus_bootstrap.json','w').write(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
