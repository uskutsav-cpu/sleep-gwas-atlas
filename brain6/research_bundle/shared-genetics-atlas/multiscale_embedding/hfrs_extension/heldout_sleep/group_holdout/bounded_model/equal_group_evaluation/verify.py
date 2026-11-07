import csv,math,json
from collections import Counter
from pathlib import Path
p=Path(__file__).resolve().parent
def read(f):
 with f.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
a=read(p.parent/'predictions.tsv');b=read(p/'summary.tsv');err=0
for row in b:
 r=[x for x in a if row['target']=='ALL' or row['target']==x['target']]; key=[(x['target'],x['held_group']) for x in r];counts=Counter(key); assert len(counts)==(24 if row['target']=='ALL' else 8)
 weights=[1/(len(counts)*counts[k]) for k in key]; mse=sum(w*(float(x[row['model']])-float(x['observed']))**2 for w,x in zip(weights,r));mae=sum(w*abs(float(x[row['model']])-float(x['observed'])) for w,x in zip(weights,r));err=max(err,abs(math.sqrt(mse)-float(row['equal_group_rmse'])),abs(mae-float(row['equal_group_mae'])))
assert err<1e-12 and len(b)==16
(p/'verification.json').write_text(json.dumps({'summary_rows':len(b),'max_absolute_error':err,'method':'independent row-level inverse group-size weights'},indent=2)+'\n');print('Verified 16 summaries, max error',err)
