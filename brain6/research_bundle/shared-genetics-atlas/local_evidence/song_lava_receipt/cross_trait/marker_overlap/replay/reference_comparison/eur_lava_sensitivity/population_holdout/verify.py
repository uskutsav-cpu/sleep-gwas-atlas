import json,math,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent;r=json.loads((p/'results.json').read_text());h=r['holdouts']
assert [x['excluded_population'] for x in h]==['CEU','FIN','GBR','IBS','TSI']
assert sum(503-x['n_people'] for x in h)==503
assert hashlib.sha256((p.parent/'eur.py').read_bytes()).hexdigest()==r['source_sha256']
for x in h:
 assert x['n_snps']<=2738 and 0<x['K']<x['n_people']
 assert math.isclose(x['delta_sumsq']/x['K']-x['sigma'],x['fi_omega'],abs_tol=1e-15)
 assert math.isclose(x['delta_sumsq']-x['K']*x['sigma'],x['fi_h2_uncapped'],abs_tol=1e-15)
print(json.dumps({'holdouts':len(h),'n_people_range':[min(x['n_people'] for x in h),max(x['n_people'] for x in h)],
 'K_range':[min(x['K'] for x in h),max(x['K'] for x in h)],
 'omega_range':[min(x['fi_omega'] for x in h),max(x['fi_omega'] for x in h)],
 'negative':sum(x['fi_omega']<0 for x in h),
 'scope':'saved arithmetic and group-count checks, not independent reference validation'},indent=2))
