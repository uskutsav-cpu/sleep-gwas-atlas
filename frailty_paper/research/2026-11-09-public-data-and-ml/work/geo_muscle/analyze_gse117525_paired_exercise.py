import csv,gzip,json,collections,math
import numpy as np
from scipy import stats
from pathlib import Path
probes={'FOXP2':'93986_at','SMAD5':'4090_at','RBM6':'10180_at','MEIS1':'4211_at'}
with gzip.open('work/geo_muscle/GSE117525.txt.gz','rt') as f:
 reader=csv.reader(f,delimiter='\t');chars=[];vals={}
 for r in reader:
  if not r:continue
  if r[0]=='!Sample_geo_accession':samples=r[1:]
  if r[0]=='!Sample_characteristics_ch1':chars.append(r[1:])
  if r[0]=='!series_matrix_table_begin':break
 assert next(reader)[1:]==samples
 for r in reader:
  if r and r[0]=='!series_matrix_table_end':break
  if r and r[0] in probes.values():vals[r[0]]=np.array([float(x) for x in r[1:]])
assert len(vals)==4
m=[{} for _ in samples]
for row in chars:
 for i,v in enumerate(row):
  if ': ' in v:
   k,x=v.split(': ',1);m[i][k.lower()]=x
by=collections.defaultdict(list)
for i,x in enumerate(m):by[x['subjectid']].append(i)
pairs=[]
for indices in by.values():
 base=[i for i in indices if m[i].get('group') in ['Frail (baseline)','Healthy older (baseline)']]
 if len(base)!=1 or len(indices)!=2:continue
 before=base[0];after=next(i for i in indices if i!=before);group=int(m[before]['group']=='Frail (baseline)');expected='Frail (after training)' if group else 'Healthy older (after training)';assert m[after].get('group')==expected
 pairs.append((before,after,group,float(m[before]['age (yrs)']),int(m[before]['sex']=='M')))
assert len(pairs)==72 and sum(x[2] for x in pairs)==31
n=len(pairs);group=np.array([r[2] for r in pairs]);age=np.array([r[3] for r in pairs]);sex=np.array([r[4] for r in pairs]);X=np.column_stack([np.ones(n),group,age-age.mean(),sex]);inv=np.linalg.inv(X.T@X);hat=np.sum((X@inv)*X,axis=1);out=[]
for gene,probe in probes.items():
 y=np.array([vals[probe][post]-vals[probe][pre] for pre,post,*_ in pairs]);beta=inv@X.T@y;res=y-X@beta;meat=X.T@((res/(1-hat))[:,None]**2*X);cov=inv@meat@inv;se=math.sqrt(cov[1,1]);p=float(2*stats.t.sf(abs(beta[1]/se),df=n-4));out.append({'gene':gene,'mean_paired_change_frail':float(y[group==1].mean()),'mean_paired_change_fit':float(y[group==0].mean()),'difference_in_change_adjusted':float(beta[1]),'hc3_se':se,'p_unadjusted':p})
order=sorted(range(4),key=lambda i:out[i]['p_unadjusted']);prev=1
for k in range(3,-1,-1):
 i=order[k];prev=min(prev,out[i]['p_unadjusted']*4/(k+1));out[i]['bh_q_four_genes']=prev
Path('outputs/stage6_gse117525_paired_exercise.json').write_text(json.dumps({'paired_n_frail':31,'paired_n_fit':41,'results':out,'method':'Within-subject post minus baseline RMA expression; compare group differences in paired changes with baseline age and sex HC3 OLS; BH across four postselected genes','limits':'observational differential exercise response, different exercise frequency and cohort/study label, no randomization of groups, selection to follow-up, unverified technical batch, no direct causal or genetic score inference'},indent=2)+'\n');print(json.dumps(out,indent=2))
