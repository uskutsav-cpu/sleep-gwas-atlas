import csv,json,hashlib,math
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent; src=p.parents[3]/'source_matrix.tsv';old=p.parent/'predictions.tsv'
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
r=read(src);o=read(old);names=list(dict.fromkeys(x['sleep_trait'] for x in r));f=['adhd','mdd','scz','bipolar','alz','parkinson'];targets=['FI','FFS','HFRS'];d={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in r}; X=np.array([[d[n,k] for k in f] for n in names]);Y=np.array([[d[n,k] for k in targets] for n in names]);oldby={(x['sleep_trait'],x['target']):x for x in o};assert len(names)==12 and len(d)==108 and len(oldby)==36 and np.all(np.abs(Y)<1)
groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names};grid=[.1,1.,10.,100.]
def predict(tr,va,j,lam):
 a=X[tr];mu=a.mean(0);sd=a.std(0);sd=np.where(sd>1e-12,sd,1);a=(a-mu)/sd;v=(X[va]-mu)/sd;z=np.arctanh(Y[tr,j]);zm=z.mean();b=np.linalg.solve(a.T@a+lam*np.eye(6),a.T@(z-zm));return np.tanh(zm+v@b)
out=[]
for excluded in ['longsleep','sleep_apnea']:
 keep=[i for i,n in enumerate(names) if n!=excluded];outer=list(dict.fromkeys(groups[names[i]] for i in keep))
 for held in outer:
  test=np.array([i for i in keep if groups[names[i]]==held]);train=np.array([i for i in keep if groups[names[i]]!=held]);inner=list(dict.fromkeys(groups[names[i]] for i in train))
  for j,t in enumerate(targets):
   loss=[]
   for lam in grid:
    errs=[]
    for g in inner:
     tr=np.array([i for i in train if groups[names[i]]!=g]);va=np.array([i for i in train if groups[names[i]]==g]);errs.extend((predict(tr,va,j,lam)-Y[va,j])**2)
    loss.append(float(np.mean(errs)))
   best=min(range(len(grid)),key=lambda k:(loss[k],-grid[k]));ys=predict(train,test,j,grid[best]);ym=float(Y[train,j].mean())
   for i,v in zip(test,ys):
    oldr=oldby[names[i],t];assert math.isclose(float(oldr['observed']),Y[i,j],abs_tol=1e-12)
    out.append({'excluded':excluded,'sleep_trait':names[i],'held_group':held,'target':t,'observed':Y[i,j],'exclusion_prediction':float(v),'original_prediction':float(oldr['bounded_ridge']),'training_mean':ym,'lambda':grid[best],'n_train':len(train)})
assert len(out)==66
summary=[]
for ex in ['longsleep','sleep_apnea']:
 for t in ['ALL']+targets:
  z=[x for x in out if x['excluded']==ex and (t=='ALL' or x['target']==t)]
  for model in ['original_prediction','exclusion_prediction','training_mean']:
   err=[x[model]-x['observed'] for x in z];summary.append({'excluded':ex,'target':t,'model':model,'n':len(z),'rmse':math.sqrt(sum(e*e for e in err)/len(err)),'mae':sum(abs(e) for e in err)/len(err)})
def write(name,rows):
 with (p/name).open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
write('predictions.tsv',out);write('summary.tsv',summary)
(p/'receipt.json').write_text(json.dumps({'inputs_sha256':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in (src,old)},'rows':len(out),'grid':grid},indent=2)+'\n')
for x in summary: print(x)
