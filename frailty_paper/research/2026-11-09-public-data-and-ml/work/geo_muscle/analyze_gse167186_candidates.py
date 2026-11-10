import csv,gzip,json,hashlib,math
import numpy as np
from scipy import stats
from pathlib import Path
meta={}
with gzip.open('work/geo_muscle/GSE167186-GPL20301.txt.gz','rt') as f:
 for r in csv.reader(f,delimiter='\t'):
  if not r:continue
  if r[0]=='!Sample_title':titles=r[1:]
  if r[0]=='!Sample_characteristics_ch1':
   for i,v in enumerate(r[1:]):
    if ': ' in v:
     k,x=v.split(': ',1);meta.setdefault(k.lower(),[None]*len(titles))[i]=x
  if r[0]=='!series_matrix_table_begin':break
assert len(titles)==76 and len(set(titles))==76
lookup={v:i for i,v in enumerate(titles)};genes=['FOXP2','SMAD5','RBM6','MEIS1'];vals={};lib=None
with gzip.open('work/geo_muscle/GSE167186_counts.csv.gz','rt') as f:
 reader=csv.reader(f);head=next(reader)[1:];assert len(head)==76 and set(head)==set(titles);idx=[lookup[x] for x in head]
 for row in reader:
  if not row:continue
  a=np.array([float(x) for x in row[1:]]);assert len(a)==76 and np.all(a>=0);lib=a.copy() if lib is None else lib+a
  if row[0] in genes:
   assert row[0] not in vals;vals[row[0]]=a
assert len(vals)==4
sel=[i for i,j in enumerate(idx) if meta['group'][j] in ['Sarcopenia','Old Healthy']];n=len(sel);assert n==53
group=np.array([int(meta['group'][idx[i]]=='Sarcopenia') for i in sel]);assert group.sum()==24
age=np.array([float(meta['age'][idx[i]]) for i in sel]);X=np.column_stack([np.ones(n),group,age-age.mean()]);inv=np.linalg.inv(X.T@X);hat=np.sum(X@inv*X,axis=1);out=[]
for g in genes:
 a=vals[g][sel];y=np.log2((a+.5)/(lib[sel]+1)*1e6);beta=inv@X.T@y;res=y-X@beta;cov=inv@(X.T@((res/(1-hat))[:,None]**2*X))@inv;se=math.sqrt(cov[1,1]);pv=float(2*stats.t.sf(abs(beta[1]/se),df=n-3));out.append({'gene':g,'adjusted_difference':float(beta[1]),'hc3_se':se,'p_unadjusted':pv,'median_count':float(np.median(a)),'zero_count':int(sum(a==0))})
order=sorted(range(4),key=lambda i:out[i]['p_unadjusted']);prev=1
for k in range(3,-1,-1):
 i=order[k];prev=min(prev,out[i]['p_unadjusted']*4/(k+1));out[i]['bh_q_four_genes']=prev
Path('outputs/stage6_gse167186_candidate_expression.json').write_text(json.dumps({'n_sarcopenia':24,'n_healthy_older':29,'results':out,'method':'exploratory OLS log2 CPM on group and age, HC3 robust SE and BH across four selected genes','limitations':'Sex unavailable in sample metadata, batch and cell composition unadjusted, small postselected, raw count model; sarcopenia not equivalent to FI/FFS; no insomnia or score validation'},indent=2)+'\n');print(json.dumps(out,indent=2))
