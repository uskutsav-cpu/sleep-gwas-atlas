import json,hashlib,math
from pathlib import Path
p=Path(__file__).resolve().parent
s=json.loads((p/'spectrum.json').read_text());v=s['eigenvalues'];assert len(v)==s['n_snps']==2741
assert all(v[i]>=v[i+1] for i in range(len(v)-1))
assert sum(x < -1e-8 for x in v)==s['negative_below_1e_8']
assert sum(x < -1e-4 for x in v)==s['negative_below_1e_4']
cum=[];total=0
for x in v:total+=x;cum.append(total)
k=next(i+1 for i,x in enumerate(cum) if 100*x/max(cum)>=99)
assert k==s['threshold_99_K']==216 and all(x>0 for x in v[:k])
assert math.isclose(v[k-1],s['smallest_retained_eigenvalue'],rel_tol=1e-12)
assert 0<=s['offdiagonal_abs_gt_1']<=2741*2740//2
print(json.dumps({'n':len(v),'K':k,'negative_below_1e_8':s['negative_below_1e_8'],
 'offdiagonal_abs_gt_1':s['offdiagonal_abs_gt_1'],'sha256':hashlib.sha256((p/'spectrum.json').read_bytes()).hexdigest(),
 'scope':'saved-spectrum arithmetic, not an independent reconstruction or LAVA calibration'},indent=2))
