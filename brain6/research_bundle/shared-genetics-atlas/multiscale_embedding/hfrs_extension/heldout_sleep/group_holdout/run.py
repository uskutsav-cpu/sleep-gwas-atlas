import csv,hashlib,json,math
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;source=p.parents[1]/'source_matrix.tsv';cross=p.parents[2]/'sleep_source_crosswalk.tsv';old=p.parent/'predictions.tsv'
with source.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
with cross.open() as f: cr=list(csv.DictReader(f,delimiter='\t'))
names=list(dict.fromkeys(r['sleep_trait'] for r in rows)); features=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS'];assert len(names)==12 and len(rows)==108
lookup={(r['sleep_trait'],r['outcome']):float(r['rg']) for r in rows};assert len(lookup)==108
X=np.array([[lookup[(n,k)] for k in features] for n in names]);Y=np.array([[lookup[(n,k)] for k in targets] for n in names]);assert np.isfinite(X).all() and np.isfinite(Y).all()
assert all(next(r for r in cr if r['sleep_trait']==n)['brain6_source_id'].startswith('dashti_2019_') for n in ('sleepdur','shortsleep','longsleep'))
assert len({next(r for r in cr if r['sleep_trait']==n)['brain6_source_id'] for n in ('sleep_efficiency','accel_sleep_duration','sleep_timing')})==1
groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names}
lamgrid=[.1,1.,10.,100.]
def predict(tx,ty,vx,lam):
 mu=tx.mean(axis=0);sd=tx.std(axis=0);sd=np.where(sd>1e-12,sd,1);a=(tx-mu)/sd;v=(vx-mu)/sd;ym=ty.mean();coef=np.linalg.solve(a.T@a+lam*np.eye(a.shape[1]),a.T@(ty-ym));return ym+v@coef
previous={(r['held_sleep'],r['target']):r for r in csv.DictReader(old.open(),delimiter='\t')};pred=[];scores=[]
for group in dict.fromkeys(groups.values()):
 test=np.array([i for i,n in enumerate(names) if groups[n]==group]);train=np.array([i for i,n in enumerate(names) if groups[n]!=group]);ingroup=list(dict.fromkeys(groups[names[i]] for i in train))
 mu=X[train].mean(axis=0);sd=X[train].std(axis=0);sd=np.where(sd>1e-12,sd,1)
 neigh={i:train[np.argsort(np.sum(((X[train]-X[i])/sd)**2,axis=1),kind='stable')[:3]] for i in test}
 for j,target in enumerate(targets):
  err=[]
  for lam in lamgrid:
   inner=[]
   for g in ingroup:
    tr=np.array([i for i in train if groups[names[i]]!=g]);val=np.array([i for i in train if groups[names[i]]==g]);v=predict(X[tr],Y[tr,j],X[val],lam);inner.extend((v-Y[val,j])**2)
   mse=float(np.mean(inner));err.append(mse);scores.append({'held_group':group,'target':target,'lambda':lam,'inner_mse':mse})
  k=min(range(len(lamgrid)),key=lambda i:(err[i],-lamgrid[i]));v=predict(X[train],Y[train,j],X[test],lamgrid[k])
  for idx,guess in zip(test,v):
   prev=previous[(names[idx],target)];assert math.isclose(float(prev['observed']),Y[idx,j],abs_tol=1e-12)
   pred.append({'sleep_trait':names[idx],'held_group':group,'target':target,'observed':Y[idx,j],'training_mean':Y[train,j].mean(),'knn3':Y[neigh[idx],j].mean(),'knn3_neighbors':','.join(names[i] for i in neigh[idx]),'ridge':guess,'lambda':lamgrid[k],'old_row_holdout_ridge':float(prev['ridge_prediction']),'n_train':len(train)})
assert len(pred)==36
summary=[]
for target in ['ALL']+targets:
 a=[r for r in pred if target=='ALL' or r['target']==target]
 for model in ('training_mean','knn3','ridge','old_row_holdout_ridge'):
  e=np.array([r[model]-r['observed'] for r in a]);summary.append({'target':target,'model':model,'n':len(a),'rmse':float(np.sqrt(np.mean(e**2))),'mae':float(np.mean(abs(e))),'out_of_range':int(sum(abs(r[model])>1 for r in a))})
def write(name,data):
 with (p/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
write('predictions.tsv',pred);write('nested_scores.tsv',scores);write('summary.tsv',summary)
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
(p/'receipt.json').write_text(json.dumps({'input_sha256':{str(f):sha(f) for f in (source,cross,old)},'groups':groups,'nested_cv':'inner leave-one-group-out; outer leave-one-group-out','lambda_grid':lamgrid,'rows':len(pred),'warning':'same-cohort leakage between groups is possible; six known r_g features required; post-hoc grouping'},indent=2)+'\n')
print(json.dumps(summary,indent=2))
