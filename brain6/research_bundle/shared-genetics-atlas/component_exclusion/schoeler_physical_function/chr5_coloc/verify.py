import csv,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
with (p/'posteriors.tsv').open() as f: expected=list(csv.DictReader(f,delimiter='\t'))
maximum=0.
for trait in ('insomnia','adhd'):
 with (p/(trait+'_matched.tsv')).open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
 def abf(beta,se,sd):
  v=se*se;r=sd*sd/(sd*sd+v);return math.exp(.5*(math.log1p(-r)+r*(beta/se)**2))
 a=[abf(float(r['physical_beta']),float(r['physical_se']),.15) for r in rows]
 b=[abf(float(r['trait_beta']),float(r['trait_se']),.2) for r in rows]
 for row in (x for x in expected if x['trait']==trait):
  same=sum(x*y for x,y in zip(a,b));weights=[1.,1e-4*sum(a),1e-4*sum(b),1e-8*(sum(a)*sum(b)-same),float(row['p12'])*same];total=sum(weights)
  for j,w in enumerate(weights):maximum=max(maximum,abs(w/total-float(row['H'+str(j)])))
assert maximum<1e-12
(p/'verification.json').write_text(json.dumps({'maximum_absolute_difference':maximum,'scope':'Independent Python direct-ABF arithmetic; not statistical or input validation.'},indent=2)+'\n')
print('maximum difference',maximum)
