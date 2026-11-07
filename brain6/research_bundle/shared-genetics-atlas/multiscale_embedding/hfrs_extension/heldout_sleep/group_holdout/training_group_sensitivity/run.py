import csv, hashlib, json
from pathlib import Path
import numpy as np
here=Path(__file__).resolve().parent
source=here.parents[2]/'source_matrix.tsv'; orig=here.parent/'predictions.tsv'
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
rows=read(source); old=read(orig); names=list(dict.fromkeys(x['sleep_trait'] for x in rows)); features=['adhd','mdd','scz','bipolar','alz','parkinson']; targets=['FI','FFS','HFRS']; lookup={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in rows}; assert len(names)==12 and len(lookup)==108
X=np.array([[lookup[n,k] for k in features] for n in names]);Y=np.array([[lookup[n,k] for k in targets] for n in names]);groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names}; allgroups=list(dict.fromkeys(groups.values())); assert len(allgroups)==8
oldby={(x['sleep_trait'],x['target']):x for x in old}; assert len(oldby)==36
lams=[.1,1.,10.,100.]
def predict(tx,ty,vx,lam):
 mean=tx.mean(axis=0); sd=tx.std(axis=0);sd=np.where(sd>1e-12,sd,1); a=(tx-mean)/sd;v=(vx-mean)/sd;ym=ty.mean(); return ym+v@np.linalg.solve(a.T@a+lam*np.eye(a.shape[1]),a.T@(ty-ym))
def fit(txids,testids,j):
 gs=list(dict.fromkeys(groups[names[i]] for i in txids)); errors=[]
 for lam in lams:
  e=[]
  for g in gs:
   tr=np.array([i for i in txids if groups[names[i]]!=g]);va=np.array([i for i in txids if groups[names[i]]==g]);e.extend((predict(X[tr],Y[tr,j],X[va],lam)-Y[va,j])**2)
  errors.append(np.mean(e))
 k=min(range(len(lams)),key=lambda i:(errors[i],-lams[i]));return predict(X[txids],Y[txids,j],X[testids],lams[k]),lams[k]
out=[]
for g in allgroups:
 test=np.array([i for i,n in enumerate(names) if groups[n]==g]); train=np.array([i for i,n in enumerate(names) if groups[n]!=g]); ingroups=list(dict.fromkeys(groups[names[i]] for i in train))
 for j,target in enumerate(targets):
  check,lam=fit(train,test,j)
  for i,guess in zip(test,check):assert abs(guess-float(oldby[names[i],target]['ridge']))<1e-10 and lam==float(oldby[names[i],target]['lambda'])
  for excluded in ingroups:
   sub=np.array([i for i in train if groups[names[i]]!=excluded]); guesses,reg=fit(sub,test,j)
   for i,guess in zip(test,guesses):
    prev=float(oldby[names[i],target]['ridge']); out.append({'sleep_trait':names[i],'held_group':g,'target':target,'excluded_group':excluded,'n_train':len(sub),'observed':Y[i,j],'full_prediction':prev,'deleted_prediction':guess,'shift':guess-prev,'lambda':reg})
summary=[]
for n in names:
 for target in targets:
  rr=[x for x in out if x['sleep_trait']==n and x['target']==target]; worst=max(rr,key=lambda x:abs(x['shift'])); vals=[x['deleted_prediction'] for x in rr];summary.append({'sleep_trait':n,'target':target,'observed':rr[0]['observed'],'full_prediction':rr[0]['full_prediction'],'min_deleted':min(vals),'max_deleted':max(vals),'largest_absolute_shift':abs(worst['shift']),'most_influential_group':worst['excluded_group'],'deleted_fits':len(rr)})
def write(file,rr):
 with (here/file).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rr[0].keys(),delimiter='\t');w.writeheader();w.writerows(rr)
write('deletions.tsv',out);write('summary.tsv',summary)
(here/'receipt.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'baseline_sha256':hashlib.sha256(orig.read_bytes()).hexdigest(),'deletion_fits':len(out),'predictions':len(summary),'max_shift':max(x['largest_absolute_shift'] for x in summary),'baseline_reproduced':36},indent=2)+'\n')
print(json.dumps({'deletion_fits':len(out),'largest_shifts':sorted(summary,key=lambda x:-x['largest_absolute_shift'])[:5]},indent=2))
