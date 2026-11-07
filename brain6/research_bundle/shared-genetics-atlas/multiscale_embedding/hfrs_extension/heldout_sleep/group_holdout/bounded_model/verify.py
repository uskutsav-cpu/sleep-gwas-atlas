import csv,json
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
src=read(p.parents[2]/'source_matrix.tsv'); rows=read(p/'predictions.tsv'); names=list(dict.fromkeys(x['sleep_trait'] for x in src)); feats=['adhd','mdd','scz','bipolar','alz','parkinson']; targets=['FI','FFS','HFRS']; lookup={(x['sleep_trait'],x['outcome']):float(x['rg']) for x in src}; X=np.array([[lookup[n,k] for k in feats] for n in names]);Y=np.array([[lookup[n,k] for k in targets] for n in names]);groups={n:('dashti_2019_duration' if n in ('sleepdur','shortsleep','longsleep') else 'jones_accelerometer' if n in ('sleep_efficiency','accel_sleep_duration','sleep_timing') else n) for n in names}
def pred(tr,va,j,lam):
 a=X[tr];mu=a.mean(0);sd=a.std(0);sd[sd<1e-12]=1;a=(a-mu)/sd;v=(X[va]-mu)/sd;z=np.arctanh(Y[tr,j]);m=z.mean();b=np.linalg.lstsq(np.vstack((a,np.sqrt(lam)*np.eye(6))),np.r_[z-m,np.zeros(6)],rcond=None)[0];return np.tanh(m+v@b)
grid=[.1,1.,10.,100.];err=0
for row in rows:
 j=targets.index(row['target']);tr=np.array([i for i,n in enumerate(names) if groups[n]!=row['held_group']]);test=np.array([names.index(row['sleep_trait'])]); gs=list(dict.fromkeys(groups[names[i]] for i in tr));loss=[]
 for l in grid:
  e=[]
  for g in gs:
   a=np.array([i for i in tr if groups[names[i]]!=g]);b=np.array([i for i in tr if groups[names[i]]==g]);e.extend((pred(a,b,j,l)-Y[b,j])**2)
  loss.append(np.mean(e))
 k=min(range(4),key=lambda i:(loss[i],-grid[i]));assert grid[k]==float(row['lambda']);err=max(err,abs(float(pred(tr,test,j,grid[k])[0])-float(row['bounded_ridge'])))
assert len(rows)==36 and err<1e-10
(p/'verification.json').write_text(json.dumps({'rows':len(rows),'max_abs_difference':err,'method':'augmented least squares with independent nested search'},indent=2)+'\n');print('Verified',len(rows),'nested predictions; max difference',err)
