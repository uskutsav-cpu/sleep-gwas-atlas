import hashlib,json,math,statistics
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parent;r=json.loads((p/'results.json').read_text());rows=r['replicates']
assert r['seed']==46420261028 and len(rows)==30 and [x['replicate'] for x in rows]==list(range(1,31))
assert hashlib.sha256((p.parent/'eur.py').read_bytes()).hexdigest()==r['source_sha256']
for x in rows:
 assert 0<x['n_snps']<=2738 and 0<x['K']<503
 assert math.isclose(x['delta_sumsq']/x['K']-x['sigma'],x['fi_omega'],abs_tol=1e-15)
 assert math.isclose(x['delta_sumsq']-x['K']*x['sigma'],x['fi_h2_uncapped'],abs_tol=1e-15)
rng=np.random.default_rng(r['seed']);unique=[len(set(rng.integers(0,503,size=503).tolist())) for _ in rows]
def summary(xs):return {'min':min(xs),'median':statistics.median(xs),'max':max(xs)}
print(json.dumps({'replicates':len(rows),'K':summary([x['K'] for x in rows]),
 'omega':summary([x['fi_omega'] for x in rows]),'distinct_sampled_people':summary(unique),
 'negative':sum(x['fi_omega']<0 for x in rows),
 'scope':'saved arithmetic and regenerated draw counts; not independent source or calibrated bootstrap verification'},indent=2))
