import csv,gzip,json,hashlib,math,re
import numpy as np
from scipy import stats
from pathlib import Path
meta={}
with gzip.open('work/geo_muscle/GSE144304.txt.gz','rt') as f:
 for r in csv.reader(f,delimiter='\t'):
  if not r:continue
  if r[0]=='!Sample_title': titles=r[1:]
  if r[0]=='!Sample_characteristics_ch1':
   for i,v in enumerate(r[1:]):
    if ': ' in v:
     k,val=v.split(': ',1);meta.setdefault(k,[None]*len(titles))[i]=val
  if r[0]=='!series_matrix_table_begin':break
lookup={re.match(r'^(s\d+)',t).group(1):i for i,t in enumerate(titles)};assert len(lookup)==80
ens={'FOXP2':'ENSG00000128573','SMAD5':'ENSG00000113658','RBM6':'ENSG00000004534','MEIS1':'ENSG00000143995'};counts={};lib=None
with gzip.open('work/geo_muscle/GSE144304_counts.txt.gz','rt') as f:
 reader=csv.reader(f,delimiter='\t');heads=next(reader)[1:];idx=[lookup[re.match(r'^(s\d+)',h).group(1)] for h in heads];assert len(idx)==80 and len(set(idx))==80
 for r in reader:
  if not r:continue
  a=np.array([float(x) for x in r[1:]],dtype=float);assert len(a)==80 and np.all(a>=0)
  lib=a.copy() if lib is None else lib+a
  gene=r[0].split('.')[0]
  if gene in ens.values():counts[gene]=a
sel=np.array([i for i,j in enumerate(idx) if meta['group'][j] in ('Frail','Fit')]);n=len(sel);assert n==54
frail=np.array([int(meta['group'][idx[i]]=='Frail') for i in sel]);age=np.array([float(meta['age'][idx[i]]) for i in sel]);sex=np.array([int(meta['gender'][idx[i]]=='Male') for i in sel]);X=np.column_stack([np.ones(n),frail,age-age.mean(),sex]);inv=np.linalg.inv(X.T@X);hat=np.sum((X@inv)*X,axis=1);out=[]
for name,g in ens.items():
 if g not in counts:out.append({'gene':name,'status':'not_in_count_matrix'});continue
 a=counts[g][sel];cpm=(a+0.5)/(lib[sel]+1)*1e6;y=np.log2(cpm);beta=inv@X.T@y;res=y-X@beta;meat=X.T@((res/(1-hat))[:,None]**2*X);cov=inv@meat@inv;se=math.sqrt(cov[1,1]);t=beta[1]/se;pv=2*stats.t.sf(abs(t),df=n-X.shape[1]);out.append({'gene':name,'ensembl':g,'n_frail_or_prefrail':int(frail.sum()),'n_fit':int(n-frail.sum()),'median_raw_count':float(np.median(a)),'n_raw_count_zero':int((a==0).sum()),'log2cpm_adjusted_difference':float(beta[1]),'hc3_se':se,'p_unadjusted':float(pv),'low_count_caution':bool((a==0).sum()>n/2 or np.median(a)<10)})
valid=[v for v in out if 'p_unadjusted' in v];rank=sorted(range(len(valid)),key=lambda i:valid[i]['p_unadjusted']);prev=1
for k in range(len(rank)-1,-1,-1):
 i=rank[k];prev=min(prev,valid[i]['p_unadjusted']*len(valid)/(k+1));valid[i]['bh_q_four_genes']=prev
Path('outputs/stage6_gse144304_candidate_expression.json').write_text(json.dumps({'source':'GSE144304 JoinedCounts and GEO sample metadata','results':out,'method':'exploratory log2 CPM, HC3 OLS adjusted age and sex, 24 (pre-)frail vs 30 fit older; BH across the four selected genes','limitations':'small, postselected, no count-variance modeling, batch/cell composition unadjusted, no genetic linkage or sleep phenotype; association not mediation, causality or score validation'},indent=2)+'\n');print(json.dumps(out,indent=2))
