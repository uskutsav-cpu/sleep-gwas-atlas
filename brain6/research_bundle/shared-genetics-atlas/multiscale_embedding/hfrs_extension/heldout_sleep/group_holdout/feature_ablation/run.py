#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent; parent=p.parent
source=parent.parents[1]/'source_matrix.tsv'; previous=parent/'predictions.tsv';groupfile=parent/'receipt.json'
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
rows=read(source);old={(r['sleep_trait'],r['target']):r for r in read(previous)};groups=json.loads(groupfile.read_text())['groups']
names=list(dict.fromkeys(r['sleep_trait'] for r in rows));features=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS']
assert len(rows)==108 and len(names)==12 and set(names)==set(groups)
lookup={(r['sleep_trait'],r['outcome']):float(r['rg']) for r in rows};assert len(lookup)==108
X=np.array([[lookup[(n,k)] for k in features] for n in names]);Y=np.array([[lookup[(n,k)] for k in targets] for n in names]);assert np.isfinite(X).all() and np.isfinite(Y).all()
grid=[.1,1.,10.,100.]
def predict(x,y,test,lam):
 mu=x.mean(axis=0);sd=x.std(axis=0);sd=np.where(sd>1e-12,sd,1);a=(x-mu)/sd;b=(test-mu)/sd;ym=y.mean();return ym+b@np.linalg.solve(a.T@a+lam*np.eye(a.shape[1]),a.T@(y-ym))
out=[]
for dropped in ['none']+features:
 cols=[i for i,n in enumerate(features) if n!=dropped];Z=X[:,cols]
 for group in dict.fromkeys(groups[n] for n in names):
  test=np.array([i for i,n in enumerate(names) if groups[n]==group]);train=np.array([i for i,n in enumerate(names) if groups[n]!=group]);innergroups=list(dict.fromkeys(groups[names[i]] for i in train))
  for j,target in enumerate(targets):
   losses=[]
   for lam in grid:
    errors=[]
    for g in innergroups:
     tr=np.array([i for i in train if groups[names[i]]!=g]);val=np.array([i for i in train if groups[names[i]]==g]);errors.extend((predict(Z[tr],Y[tr,j],Z[val],lam)-Y[val,j])**2)
    losses.append(float(np.mean(errors)))
   chosen=min(range(len(grid)),key=lambda i:(losses[i],-grid[i]));pred=predict(Z[train],Y[train,j],Z[test],grid[chosen])
   for idx,guess in zip(test,pred):
    if dropped=='none':assert abs(float(guess)-float(old[(names[idx],target)]['ridge']))<1e-10
    out.append({'dropped_feature':dropped,'sleep_trait':names[idx],'held_group':group,'target':target,'observed':Y[idx,j],'prediction':float(guess),'lambda':grid[chosen]})
assert len(out)==7*36
def metrics(sub):
 e=np.array([x['prediction']-x['observed'] for x in sub]);return float(np.sqrt(np.mean(e*e))),float(np.mean(np.abs(e)))
summ=[]
for target in ['ALL']+targets:
 baseline=[r for r in out if r['dropped_feature']=='none' and (target=='ALL' or r['target']==target)];base_rmse,_=metrics(baseline)
 for dropped in ['none']+features:
  subset=[r for r in out if r['dropped_feature']==dropped and (target=='ALL' or r['target']==target)];rmse,mae=metrics(subset)
  summ.append({'target':target,'dropped_feature':dropped,'n':len(subset),'rmse':rmse,'mae':mae,'rmse_change_vs_all_features':rmse-base_rmse,'out_of_range':sum(abs(r['prediction'])>1 for r in subset)})
def write(path,data):
 with path.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
write(p/'predictions.tsv',out);write(p/'summary.tsv',summ)
sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest()
receipt={'input_hashes':{str(x):sha(x) for x in (source,previous,groupfile)},'rows':len(out),'full_feature_predictions_match_previous':True,'lambda_grid':grid,'scope':'Post-hoc ablation on 12 correlated sleep measures. No proof of feature causality, external generalization, or individual risk.'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(summ,indent=2))
