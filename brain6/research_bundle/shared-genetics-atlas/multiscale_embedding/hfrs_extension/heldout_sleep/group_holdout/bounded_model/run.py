import csv,json,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent; source=p.parents[2]/'source_matrix.tsv';orig=p.parent/'predictions.tsv'
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
r=read(source); old=read(orig); names=list(dict.fromkeys(x['sleep_trait'] for x in r)); feats=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS']; lookup={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in r}; assert len(names)==12 and len(lookup)==108
X=np.array([[lookup[n,k] for k in feats] for n in names]);Y=np.array([[lookup[n,k] for k in targets] for n in names]);assert np.isfinite(X).all() and np.isfinite(Y).all() and np.all(abs(Y)<1)
groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names}; grid=[.1,1.,10.,100.]; oldby={(x['sleep_trait'],x['target']):x for x in old}; assert len(oldby)==36
def predict(train,valid,j,lam):
 a=X[train]; mu=a.mean(axis=0);sd=a.std(axis=0);sd=np.where(sd>1e-12,sd,1);a=(a-mu)/sd;v=(X[valid]-mu)/sd;z=np.arctanh(Y[train,j]);zm=z.mean();b=np.linalg.solve(a.T@a+lam*np.eye(6),a.T@(z-zm));return np.tanh(zm+v@b)
out=[]
for g in dict.fromkeys(groups.values()):
 tr=np.array([i for i,n in enumerate(names) if groups[n]!=g]);te=np.array([i for i,n in enumerate(names) if groups[n]==g]);inner=list(dict.fromkeys(groups[names[i]] for i in tr))
 for j,target in enumerate(targets):
  mse=[]
  for lam in grid:
   e=[]
   for gi in inner:
    sub=np.array([i for i in tr if groups[names[i]]!=gi]); val=np.array([i for i in tr if groups[names[i]]==gi]);e.extend((predict(sub,val,j,lam)-Y[val,j])**2)
   mse.append(float(np.mean(e)))
  k=min(range(4),key=lambda i:(mse[i],-grid[i])); pred=predict(tr,te,j,grid[k])
  for i,v in zip(te,pred):
   o=oldby[names[i],target];assert abs(float(o['observed'])-Y[i,j])<1e-12
   out.append({'sleep_trait':names[i],'held_group':g,'target':target,'observed':Y[i,j],'training_mean':float(o['training_mean']),'knn3':float(o['knn3']),'ridge':float(o['ridge']),'bounded_ridge':v,'lambda':grid[k]})
summary=[]
for target in ['ALL']+targets:
 a=[x for x in out if target=='ALL' or x['target']==target]
 for model in ['training_mean','knn3','ridge','bounded_ridge']:
  e=np.array([x[model]-x['observed'] for x in a]); summary.append({'target':target,'model':model,'n':len(a),'rmse':float(np.sqrt(np.mean(e*e))),'mae':float(np.mean(abs(e))),'out_of_range':int(sum(abs(x[model])>1 for x in a))})
def write(name,a):
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=a[0].keys(),delimiter='\t');w.writeheader();w.writerows(a)
write('predictions.tsv',out);write('summary.tsv',summary)
(p/'receipt.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'baseline_sha256':hashlib.sha256(orig.read_bytes()).hexdigest(),'n_predictions':len(out),'selection_metric':'nested group-held-out MSE on r scale'},indent=2)+'\n')
print(json.dumps(summary,indent=2))
