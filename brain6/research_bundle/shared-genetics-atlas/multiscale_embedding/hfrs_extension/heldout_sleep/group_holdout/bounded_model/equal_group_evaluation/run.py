import csv,hashlib,json,math
from pathlib import Path
from collections import defaultdict
p=Path(__file__).resolve().parent;src=p.parent/'predictions.tsv'
with src.open(newline='') as f:r=list(csv.DictReader(f,delimiter='\t'))
models=['training_mean','knn3','ridge','bounded_ridge']; targets=['FI','FFS','HFRS']; assert len(r)==36 and len({(x['sleep_trait'],x['target']) for x in r})==36 and len({x['held_group'] for x in r})==8 and len({x['sleep_trait'] for x in r})==12
loss=[]
for t in targets:
 for g in sorted({x['held_group'] for x in r}):
  z=[x for x in r if x['target']==t and x['held_group']==g];assert len(z) in (1,3)
  for m in models:
   e=[float(x[m])-float(x['observed']) for x in z];assert all(math.isfinite(x) for x in e);loss.append({'target':t,'held_group':g,'n_traits':len(z),'model':m,'mse':sum(x*x for x in e)/len(e),'mae':sum(abs(x) for x in e)/len(e)})
assert len(loss)==96
summary=[]
for t in ['ALL']+targets:
 for m in models:
  z=[x for x in loss if (t=='ALL' or x['target']==t) and x['model']==m];a=[x for x in r if t=='ALL' or x['target']==t];err=[float(x[m])-float(x['observed']) for x in a]
  summary.append({'target':t,'model':m,'equal_group_rmse':math.sqrt(sum(x['mse'] for x in z)/len(z)),'equal_group_mae':sum(x['mae'] for x in z)/len(z),'pooled_rmse':math.sqrt(sum(x*x for x in err)/len(err)),'pooled_mae':sum(abs(x) for x in err)/len(err),'groups':8,'rows':len(a)})
def write(name,rows):
 with (p/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0],delimiter='\t');w.writeheader();w.writerows(rows)
write('group_loss.tsv',loss);write('summary.tsv',summary)
(p/'receipt.json').write_text(json.dumps({'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'rows':len(r),'groups':8,'targets':3},indent=2)+'\n')
for x in summary:print(x)
