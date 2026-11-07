import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;parent=p.parent
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
r=json.loads((parent/'match_receipt.json').read_text())
rows={}
for trait in ('insomnia','adhd'):
 path=parent/(trait+'_matched.tsv');assert sha(path)==r['results'][trait]['matched_sha256']
 with path.open() as f:rows[trait]={int(x['pos38']):x for x in csv.DictReader(f,delimiter='\t')}
common=set(rows['insomnia']) & set(rows['adhd'])
with (p/'positions.tsv').open() as f:stored={int(x['pos38']) for x in csv.DictReader(f,delimiter='\t')}
assert common==stored and len(common)==2359
expected=list(csv.DictReader((p/'posteriors.tsv').open(),delimiter='\t'));maxdiff=0
for trait in ('insomnia','adhd'):
 d=[rows[trait][k] for k in sorted(common)]
 def abf(beta,se,sd):
  v=se*se;r=sd*sd/(sd*sd+v);return math.exp(.5*(math.log1p(-r)+r*(beta/se)**2))
 a=[abf(float(x['physical_beta']),float(x['physical_se']),.15) for x in d];b=[abf(float(x['trait_beta']),float(x['trait_se']),.2) for x in d];same=sum(x*y for x,y in zip(a,b))
 for r in (x for x in expected if x['trait']==trait):
  w=[1.,1e-4*sum(a),1e-4*sum(b),1e-8*(sum(a)*sum(b)-same),float(r['p12'])*same];tot=sum(w)
  for j,value in enumerate(w):maxdiff=max(maxdiff,abs(value/tot-float(r['H'+str(j)])))
assert maxdiff<1e-12
out={'common_n':len(common),'dropped_insomnia':len(rows['insomnia'])-len(common),'dropped_adhd':len(rows['adhd'])-len(common),'maximum_absolute_difference':maxdiff,'scope':'Checks input hashes, common positions and ABF arithmetic; not statistical validation.'}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
