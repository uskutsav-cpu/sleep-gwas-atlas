import json,math,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent
r=json.loads((p/'results.json').read_text());prev=json.loads((p.parent.parent/'insomnia.json').read_text())
assert {x['removed'] for x in r}=={'rs62262671','rs55754265'}
for x in r:
 assert x['n_snps']==prev['n_snps']-1 and x['K']==216
 assert math.isclose(x['fi_delta_sumsq']/x['K']-x['fi_sigma'],x['fi_omega'],abs_tol=1e-15)
 assert x['fi_failed'] and x['fi_omega']<0 and x['fi_h2_capped']==0
print(json.dumps({'runs':len(r),'both_failed':True,'baseline_fi_omega':prev['omega'][0][0],
 'results_sha256':hashlib.sha256((p/'results.json').read_bytes()).hexdigest(),
 'scope':'saved arithmetic and counts, not independent validation of reference or historical inputs'},indent=2))
