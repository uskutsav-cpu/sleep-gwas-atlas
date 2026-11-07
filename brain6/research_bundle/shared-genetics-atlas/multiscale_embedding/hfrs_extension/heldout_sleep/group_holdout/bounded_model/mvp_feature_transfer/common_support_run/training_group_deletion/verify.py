import csv,hashlib,json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas');repo=atlas.parent.parent/'work/sleep-gwas-atlas/brain6/research_prototypes/frailty_prediction';matrix=repo/'source_matrix.tsv';features=repo/'mvp_common_support_features.json';output=p/'predictions.tsv'
with matrix.open(newline='') as f:r=list(csv.DictReader(f,delimiter='\t'))
lookup={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in r};names=list(dict.fromkeys(x['sleep_trait'] for x in r));feat=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS'];X=np.array([[lookup[n,k] for k in feat] for n in names]);Y=np.array([[lookup[n,k] for k in targets] for n in names]);v=json.loads(features.read_text());x=np.array([v[k] for k in feat]);groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names};grid=[.1,1.,10.,100.]
def fit(tr,values,target,lam):
 a=X[tr];mu=a.mean(0);sd=a.std(0);sd=np.where(sd>1e-12,sd,1);a=(a-mu)/sd;z=np.arctanh(Y[tr,target]);coef=np.linalg.lstsq(np.vstack((a,np.sqrt(lam)*np.eye(6))),np.r_[z-z.mean(),np.zeros(6)],rcond=None)[0];return np.tanh(z.mean()+((values-mu)/sd)@coef)
with output.open(newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==21;maxerr=0
for row in rows:
 n=row['deleted_training_group'];j=targets.index(row['target']);tr=np.array([i for i,name in enumerate(names) if groups[name] not in ('insomnia',n)]);gs=list(dict.fromkeys(groups[names[i]] for i in tr));scores=[]
 for lam in grid:
  e=[]
  for g in gs:
   sub=np.array([i for i in tr if groups[names[i]]!=g]);valid=np.array([i for i in tr if groups[names[i]]==g]);e.extend((fit(sub,X[valid],j,lam)-Y[valid,j])**2)
  scores.append(np.mean(e))
 best=min(range(4),key=lambda i:(scores[i],-grid[i]));assert grid[best]==float(row['lambda']);estimate=float(fit(tr,x,j,grid[best]));maxerr=max(maxerr,abs(estimate-float(row['prediction'])));assert abs(estimate-float(row['prediction']))<1e-10
(p/'verification.json').write_text(json.dumps({'rows':len(rows),'max_absolute_prediction_difference':maxerr,'method':'independent augmented least squares and group-held-out selection','output_sha256':hashlib.sha256(output.read_bytes()).hexdigest()},indent=2)+'\n');print('verified',len(rows),'refits; max difference',maxerr)
