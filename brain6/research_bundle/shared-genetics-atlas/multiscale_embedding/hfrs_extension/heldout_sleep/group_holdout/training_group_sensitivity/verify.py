import csv,math,json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
src=read(p.parents[2]/'source_matrix.tsv'); d=read(p/'deletions.tsv'); names=list(dict.fromkeys(x['sleep_trait'] for x in src)); features=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS']; look={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in src};X=np.array([[look[n,k] for k in features] for n in names]);Y=np.array([[look[n,k] for k in targets] for n in names]); groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names}
def fit(tr,va,j,lam):
 a=X[tr];mu=a.mean(axis=0);sd=a.std(axis=0);sd[sd<1e-12]=1; a=(a-mu)/sd; v=(X[va]-mu)/sd; ym=Y[tr,j].mean(); b=np.linalg.lstsq(np.vstack((a,np.sqrt(lam)*np.eye(6))),np.r_[Y[tr,j]-ym,np.zeros(6)],rcond=None)[0];return ym+v@b
maxdiff=0
for x in d:
 tr=np.array([i for i,n in enumerate(names) if groups[n] not in (x['held_group'],x['excluded_group'])]);va=np.array([names.index(x['sleep_trait'])]);j=targets.index(x['target']);gs=list(dict.fromkeys(groups[names[i]] for i in tr)); scores=[]
 for lam in (.1,1.,10.,100.):
  e=[]
  for g in gs:
   a=np.array([i for i in tr if groups[names[i]]!=g]);b=np.array([i for i in tr if groups[names[i]]==g]);e.extend((fit(a,b,j,lam)-Y[b,j])**2)
  scores.append(float(np.mean(e)))
 grid=[.1,1.,10.,100.]; k=min(range(4),key=lambda i:(scores[i],-grid[i]));assert grid[k]==float(x['lambda'])
 val=float(fit(tr,va,j,grid[k])[0]);maxdiff=max(maxdiff,abs(val-float(x['deleted_prediction'])))
assert len(d)==252 and maxdiff<1e-9
(p/'verification.json').write_text(json.dumps({'rows':len(d),'max_absolute_prediction_difference':maxdiff,'method':'independent augmented least-squares solver'},indent=2)+'\n');print('Verified',len(d),'fits; max difference',maxdiff)
