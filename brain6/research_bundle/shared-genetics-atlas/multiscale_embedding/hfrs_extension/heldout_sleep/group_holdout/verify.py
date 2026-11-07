import csv,json,math,hashlib
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
read=lambda f:list(csv.DictReader(f.open(),delimiter='\t'))
s=read(p.parents[1]/'source_matrix.tsv');pred=read(p/'predictions.tsv');scores=read(p/'nested_scores.tsv');summ=read(p/'summary.tsv');receipt=json.loads((p/'receipt.json').read_text())
names=list(dict.fromkeys(r['sleep_trait'] for r in s));features=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS'];data={(r['sleep_trait'],r['outcome']):float(r['rg']) for r in s};X=np.array([[data[n,f] for f in features] for n in names]);Y=np.array([[data[n,t] for t in targets] for n in names]);grp=receipt['groups'];assert len(pred)==36 and len(scores)==96
def estimate(xx,y,v,lam):
 mu=xx.mean(0);sd=xx.std(0);sd=np.where(sd>1e-12,sd,1);a=np.column_stack((np.ones(len(xx)),(xx-mu)/sd));pen=np.column_stack((np.zeros(6),np.sqrt(lam)*np.eye(6)));b=np.r_[y,np.zeros(6)];coef=np.linalg.lstsq(np.vstack((a,pen)),b,rcond=None)[0];return float(np.r_[1,(v-mu)/sd]@coef)
for r in pred:
 i=names.index(r['sleep_trait']);j=targets.index(r['target']);g=grp[names[i]];assert r['held_group']==g
 train=[k for k,n in enumerate(names) if grp[n]!=g];assert int(r['n_train'])==len(train)
 assert math.isclose(float(r['observed']),Y[i,j],abs_tol=1e-12)
 candidates=[x for x in scores if x['held_group']==g and x['target']==r['target']];assert len(candidates)==4
 for c in candidates:
  err=[]
  for h in set(grp[names[k]] for k in train):
   tt=[k for k in train if grp[names[k]]!=h];vv=[k for k in train if grp[names[k]]==h]
   for k in vv:err.append((estimate(X[tt],Y[tt,j],X[k],float(c['lambda']))-Y[k,j])**2)
  assert math.isclose(sum(err)/len(err),float(c['inner_mse']),rel_tol=1e-9,abs_tol=1e-12)
 best=min(candidates,key=lambda v:(float(v['inner_mse']),-float(v['lambda'])));assert float(r['lambda'])==float(best['lambda'])
 assert math.isclose(estimate(X[train],Y[train,j],X[i],float(best['lambda'])),float(r['ridge']),rel_tol=1e-9,abs_tol=1e-12)
 assert math.isclose(float(r['training_mean']),float(Y[train,j].mean()),abs_tol=1e-12)
 mu=X[train].mean(0);sd=X[train].std(0);sd=np.where(sd>1e-12,sd,1);idx=np.array(train)[np.argsort(np.sum(((X[train]-X[i])/sd)**2,axis=1),kind='stable')[:3]]
 assert r['knn3_neighbors']==','.join(names[k] for k in idx) and math.isclose(float(r['knn3']),float(Y[idx,j].mean()),abs_tol=1e-12)
for r in summ:
 subset=[x for x in pred if r['target']=='ALL' or x['target']==r['target']];err=[float(x[r['model']])-float(x['observed']) for x in subset]; assert math.isclose(float(r['rmse']),math.sqrt(sum(x*x for x in err)/len(err)),abs_tol=1e-12)
report={'status':'PASS','predictions':len(pred),'inner_scores':len(scores),'independent_solver':'augmented least squares','checks':'groups, inner and outer predictions, neighbors, RMSE','hashes':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ['run.py','predictions.tsv','nested_scores.tsv','summary.tsv']},'limitation':'No validation of independence or source homogeneity.'}
(p/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
