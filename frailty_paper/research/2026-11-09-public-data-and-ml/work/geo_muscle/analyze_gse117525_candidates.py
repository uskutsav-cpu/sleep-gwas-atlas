import csv,gzip,json,math,hashlib
import numpy as np
from scipy import stats
from pathlib import Path
p=Path('work/geo_muscle/GSE117525.txt.gz');platform=Path('work/geo_muscle/GPL20880.soft.gz')
genes={'FOXP2':('93986_at','forkhead box P2'),'SMAD5':('4090_at','SMAD family member 5'),'RBM6':('10180_at','RNA binding motif protein 6'),'MEIS1':('4211_at','Meis homeobox 1')};validated={}
with gzip.open(platform,'rt') as f:
 for l in f:
  if l.startswith('!platform_table_begin'):break
 header=next(f).rstrip('\n').split('\t');assert header[:2]==['ID','Description']
 for l in f:
  if l.startswith('!platform_table_end'):break
  r=l.rstrip('\n').split('\t')
  for name,(probe,desc) in genes.items():
   if r[0]==probe:validated[name]=r[1];assert r[1].lower()==desc.lower()
assert len(validated)==4
characteristics=[];values={}
with gzip.open(p,'rt') as f:
 reader=csv.reader(f,delimiter='\t')
 for r in reader:
  if not r:continue
  if r[0]=='!Sample_geo_accession':samples=r[1:]
  if r[0]=='!Sample_characteristics_ch1':characteristics.append(r[1:])
  if r[0]=='!series_matrix_table_begin':break
 head=next(reader);assert head[1:]==samples
 for r in reader:
  if r and r[0]=='!series_matrix_table_end':break
  if r and r[0] in {v[0] for v in genes.values()}:values[r[0]]=np.array([float(v) for v in r[1:]])
assert len(values)==4
meta=[{} for _ in samples]
for row in characteristics:
 assert len(row)==len(samples)
 for i,v in enumerate(row):
  if ': ' in v:
   k,x=v.split(': ',1);meta[i][k.lower()]=x
idx=[i for i,m in enumerate(meta) if m.get('group') in ['Frail (baseline)','Healthy older (baseline)']];assert len(idx)==134
assert len({meta[i]['subjectid'] for i in idx})==134
frail=np.array([int(meta[i]['group']=='Frail (baseline)') for i in idx]);assert frail.sum()==61
age=np.array([float(meta[i]['age (yrs)']) for i in idx]);sex=np.array([int(meta[i]['sex']=='M') for i in idx]);assert set(meta[i]['sex'] for i in idx)=={'M','F'}
n=len(idx);X=np.column_stack([np.ones(n),frail,age-age.mean(),sex]);inv=np.linalg.inv(X.T@X);hat=np.sum((X@inv)*X,axis=1);out=[]
for name,(probe,_) in genes.items():
 y=values[probe][idx];beta=inv@X.T@y;res=y-X@beta;meat=X.T@((res/(1-hat))[:,None]**2*X);cov=inv@meat@inv;se=math.sqrt(cov[1,1]);pv=float(2*stats.t.sf(abs(beta[1]/se),df=n-X.shape[1]));out.append({'gene':name,'probe_id':probe,'adjusted_difference':float(beta[1]),'hc3_se':se,'p_unadjusted':pv})
order=sorted(range(4),key=lambda i:out[i]['p_unadjusted']);prev=1
for k in range(3,-1,-1):
 i=order[k];prev=min(prev,out[i]['p_unadjusted']*4/(k+1));out[i]['bh_q_four_genes']=prev
Path('outputs/stage6_gse117525_candidate_expression.json').write_text(json.dumps({'results':out,'baseline_n':n,'n_frail':int(frail.sum()),'n_fit':int(n-frail.sum()),'method':'GEO processed RMA expression, OLS group+age+sex with HC3 SE and t(130), BH across four postselected genes','source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'annotation_sha256':hashlib.sha256(platform.read_bytes()).hexdigest(),'limitations':'postselected candidate analysis; no batch or cell composition adjustments, age/sex model assumptions, potential cohort overlap, not genetic association, cell-specific or score validation'},indent=2)+'\n');print(json.dumps(out,indent=2))
