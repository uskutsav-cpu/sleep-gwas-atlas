import hashlib, json, math
from pathlib import Path
p=Path(__file__).resolve().parent
rows=json.loads((p/'results.json').read_text())
prior=json.loads((p.parent/'insomnia.json').read_text())
assert [x['threshold'] for x in rows]==[95,97,98,99,99.5]
assert all(x['n_snps']==prior['n_snps']==2741 for x in rows)
assert [x['K'] for x in rows]==sorted(x['K'] for x in rows)
for x in rows:
    assert math.isclose(x['fi_delta_sumsq']/x['K']-x['fi_sigma'],x['fi_omega'],abs_tol=1e-15)
    assert x['fi_failed']==(x['fi_omega']<0)
    assert math.isclose(x['fi_h2_capped'],max(0,x['K']*x['fi_omega']),abs_tol=1e-12)
x=rows[3]
assert x['K']==prior['K'] and math.isclose(x['fi_omega'],prior['omega'][0][0],abs_tol=1e-15)
print(json.dumps({'rows':len(rows),'negative':sum(x['fi_failed'] for x in rows),
 'prior_reproduced':True,'sha256':hashlib.sha256((p/'results.json').read_bytes()).hexdigest(),
 'limits':'arithmetic and saved replay consistency only, not historical input identity or statistical validation'},indent=2))
