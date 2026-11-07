#!/usr/bin/env python3
"""Descriptive SVD and repeated masked-cell experiment; no AlphaGenome scores."""
import csv,json,hashlib,pathlib,sys,platform,numpy as np
root=pathlib.Path(__file__).resolve().parent; p=root/'source_brain6_72.tsv'
rows=list(csv.DictReader(p.open(),delimiter='\t')); sleeps=list(dict.fromkeys(r['sleep_trait'] for r in rows)); diseases=list(dict.fromkeys(r['brain_disorder'] for r in rows));assert len(rows)==72 and len(sleeps)==12 and len(diseases)==6
M=np.empty((12,6));S=np.empty((12,6));seen=set()
for r in rows:
 i=sleeps.index(r['sleep_trait']);j=diseases.index(r['brain_disorder']);assert (i,j) not in seen;seen.add((i,j));M[i,j]=float(r['rg']);S[i,j]=float(r['se']);assert np.isfinite(M[i,j]) and np.isfinite(S[i,j]) and S[i,j]>0
u,d,vt=np.linalg.svd(M,full_matrices=False)
# Fix arbitrary signs for stable tables, choosing positive loading for largest absolute disease coordinate.
for k in range(len(d)):
 j=int(np.argmax(np.abs(vt[k])));sign=1 if vt[k,j]>=0 else -1;u[:,k]*=sign;vt[k]*=sign
sleep=u[:,:2]*np.sqrt(d[:2]);disease=vt[:2,:].T*np.sqrt(d[:2]);recon=sleep@disease.T
with (root/'coordinates.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['type','trait','dimension_1','dimension_2']);
 for label,names,arr in [('sleep',sleeps,sleep),('disorder',diseases,disease)]:
  for name,x in zip(names,arr):w.writerow([label,name,*x])
with (root/'reconstruction.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['sleep_trait','brain_disorder','observed_rg','rg_se','reconstructed_rank2_rg','residual']);
 for i,s in enumerate(sleeps):
  for j,t in enumerate(diseases):w.writerow([s,t,M[i,j],S[i,j],recon[i,j],M[i,j]-recon[i,j]])
def fit_als(obs,train,seed,rank=2,lam=.02,steps=100):
 rng=np.random.RandomState(seed);m,n=obs.shape;U=rng.normal(scale=.15,size=(m,rank));V=rng.normal(scale=.15,size=(n,rank));I=np.eye(rank)*lam
 for it in range(steps):
  for i in range(m):
   ids=np.flatnonzero(train[i]);A=V[ids];U[i]=np.linalg.solve(A.T@A+I,A.T@obs[i,ids])
  for j in range(n):
   ids=np.flatnonzero(train[:,j]);A=U[ids];V[j]=np.linalg.solve(A.T@A+I,A.T@obs[ids,j])
 return U@V.T
# Freeze 10 shuffled five-fold assignments before fitting; no tuning on hold-outs.
rng=np.random.RandomState(20261027);folds=[]
for rep in range(10):
 while True:
  order=rng.permutation(M.size);fold=np.empty(M.size,dtype=int);fold[order]=np.arange(M.size)%5;fold=fold.reshape(M.shape)
  if all(np.all((fold!=k).sum(axis=0)>=3) and np.all((fold!=k).sum(axis=1)>=3) for k in range(5)):break
 folds.append(fold)
cv=[]
for rep,fold in enumerate(folds):
 for k in range(5):
  train=fold!=k;test=~train;baseline=np.array([np.mean(M[train[:,j],j]) for j in range(6)])[None,:]*np.ones((12,1));pred=fit_als(M,train,seed=3000+rep*5+k);pred1=fit_als(M,train,seed=4000+rep*5+k,rank=1);design=np.zeros((M.size,1+11+5));design[:,0]=1
  for r in range(12):
   for c in range(6):
    if r>0: design[r*6+c,r]=1
    if c>0: design[r*6+c,11+c]=1
  coeff=np.linalg.lstsq(design[train.ravel()],M.ravel()[train.ravel()],rcond=None)[0];additive=(design@coeff).reshape(M.shape);ii,jj=np.where(test)
  for i,j in zip(ii,jj):cv.append({'repeat':rep,'fold':k,'sleep_trait':sleeps[i],'brain_disorder':diseases[j],'observed_rg':M[i,j],'column_mean_prediction':baseline[i,j],'rank2_prediction':pred[i,j],'rank1_prediction':pred1[i,j],'additive_prediction':additive[i,j],'rank1_sq_error':(M[i,j]-pred1[i,j])**2,'additive_sq_error':(M[i,j]-additive[i,j])**2,'baseline_sq_error':(M[i,j]-baseline[i,j])**2,'rank2_sq_error':(M[i,j]-pred[i,j])**2})
with (root/'masked_cell_predictions.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=cv[0],delimiter='\t');w.writeheader();w.writerows(cv)
summary={'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'n_sleep':12,'n_disease':6,'n_observed':72,'svd_singular_values':d.tolist(),'rank2_in_sample_squared_frobenius_fraction':float(sum(d[:2]**2)/sum(d**2)),'rank2_in_sample_rmse':float(np.sqrt(np.mean((M-recon)**2))),'masked_cell_test':{'repeats':10,'folds':5,'n_predictions':len(cv),'seed':20261027,'model':'2-component unweighted ridge ALS, lambda 0.02, 100 alternating steps, one fixed initialization per fold; no tuning','baseline':'mean of observed training rows in the same disease column','baseline_pooled_rmse':float(np.sqrt(np.mean([x['baseline_sq_error'] for x in cv]))),'rank2_pooled_rmse':float(np.sqrt(np.mean([x['rank2_sq_error'] for x in cv]))),'rank1_pooled_rmse':float(np.sqrt(np.mean([x['rank1_sq_error'] for x in cv]))),'additive_pooled_rmse':float(np.sqrt(np.mean([x['additive_sq_error'] for x in cv]))),'warning':'same 72 observed source estimates reused; repeated folds correlated; this is a missing-cell reconstruction diagnostic, not independent-cohort replication or risk prediction'},'python':platform.python_version(),'numpy':np.__version__}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
