import json,hashlib,math
from pathlib import Path
p=Path(__file__).resolve().parent
r=json.loads((p/'results.json').read_text());prior=json.loads((p.parent.parent/'insomnia.json').read_text())
assert [x['setting'] for x in r]==['unchanged','clipped','ridge']
assert all(x['n_snps']==2741 for x in r)
assert r[0]['K']==prior['K'] and math.isclose(r[0]['fi_omega'],prior['omega'][0][0],abs_tol=1e-15)
for x in r:
 assert math.isclose(x['fi_delta_sumsq']/x['K']-x['fi_sigma'],x['fi_omega'],abs_tol=1e-15)
 assert x['fi_failed']==(x['fi_omega']<0)
 assert math.isclose(x['fi_h2_capped'],max(0,x['K']*x['fi_omega']),abs_tol=1e-12)
 if x['setting']!='unchanged':assert x['minimum_eigenvalue']>-1e-7
print(json.dumps({'settings':len(r),'K':[x['K'] for x in r],'FI_failed':[x['fi_failed'] for x in r],
 'results_sha256':hashlib.sha256((p/'results.json').read_bytes()).hexdigest(),
 'scope':'saved numerical invariants only; does not validate adjusted LD or infer causation'},indent=2))
