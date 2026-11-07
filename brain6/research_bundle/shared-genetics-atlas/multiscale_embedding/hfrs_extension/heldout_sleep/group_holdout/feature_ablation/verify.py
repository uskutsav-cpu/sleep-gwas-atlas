#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;parent=p.parent
def read(x):
 with x.open() as f:return list(csv.DictReader(f,delimiter='\t'))
source=read(parent.parents[1]/'source_matrix.tsv');pred=read(p/'predictions.tsv');summary=read(p/'summary.tsv');groups=json.loads((parent/'receipt.json').read_text())['groups']
names=list(dict.fromkeys(x['sleep_trait'] for x in source));features=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS'];grid=[.1,1.,10.,100.]
d={(r['sleep_trait'],r['outcome']):float(r['rg']) for r in source};X=np.array([[d[(n,k)] for k in features] for n in names]);Y=np.array([[d[(n,k)] for k in targets] for n in names]);assert len(pred)==252
def solve(train_x,train_y,test_x,lam):
 mu=train_x.mean(0);sd=train_x.std(0);sd=np.where(sd>1e-12,sd,1);X0=(train_x-mu)/sd;xt=(test_x-mu)/sd;ym=train_y.mean()
 A=np.vstack((X0,np.sqrt(lam)*np.eye(X0.shape[1])));b=np.r_[train_y-ym,np.zeros(X0.shape[1])];return ym+xt@np.linalg.lstsq(A,b,rcond=None)[0]
maxerr=0
for row in pred:
 keep=[i for i,f in enumerate(features) if f!=row['dropped_feature']];Z=X[:,keep];target=targets.index(row['target']);ix=names.index(row['sleep_trait']);held=row['held_group'];assert groups[names[ix]]==held
 train=[i for i,n in enumerate(names) if groups[n]!=held];ig=list(dict.fromkeys(groups[names[i]] for i in train));loss=[]
 for lam in grid:
  err=[]
  for group in ig:
   tr=[i for i in train if groups[names[i]]!=group];val=[i for i in train if groups[names[i]]==group];err.extend((solve(Z[tr],Y[tr,target],Z[val],lam)-Y[val,target])**2)
  loss.append(float(np.mean(err)))
 chosen=min(range(len(grid)),key=lambda i:(loss[i],-grid[i]));assert grid[chosen]==float(row['lambda'])
 value=float(solve(Z[train],Y[train,target],Z[[ix]],grid[chosen])[0]);maxerr=max(maxerr,abs(value-float(row['prediction'])))
 assert math.isclose(float(row['observed']),Y[ix,target],abs_tol=1e-12)
assert maxerr<1e-9
for r in summary:
 x=[a for a in pred if a['dropped_feature']==r['dropped_feature'] and (r['target']=='ALL' or r['target']==a['target'])];e=[float(a['prediction'])-float(a['observed']) for a in x]
 assert len(e)==int(r['n']) and math.isclose(math.sqrt(sum(y*y for y in e)/len(e)),float(r['rmse']),abs_tol=1e-12)
receipt={'predictions_verified':len(pred),'max_prediction_error':maxerr,'lambda_choices_verified':True,'summary_rows_verified':len(summary),'sha256':{x:hashlib.sha256((p/x).read_bytes()).hexdigest() for x in ('run.py','predictions.tsv','summary.tsv')},'limitations':'Verifies numerical calculations and fold selection, not independence, biological meaning, generalization or study covariance.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
