from pathlib import Path
import csv, hashlib, json, math

root=Path('outputs/shared-genetics-atlas'); out=root/'frailty_coloc/adhd_same_support'
fi=root/'frailty_coloc/matched.tsv'; pilot=root/'variant_context/rs62264764_dense_qtl/multisignal_pilot/variants.tsv'
def read(path):
 with path.open(newline='') as h: return list(csv.DictReader(h,delimiter='\t'))
f=[r for r in read(fi) if r['chr']=='3']; p=read(pilot)
assert len({r['rsid'] for r in f})==len(f) and len({r['rsid'] for r in p})==len(p)
p={r['rsid']:r for r in p}; merged=[]
for a in f:
 b=p.get(a['rsid'])
 if not b or a['bp37']!=b['pos37'] or {a['effect_allele'],a['other_allele']}!={b['ref'],b['alt']}: continue
 if {b['ref'],b['alt']} in ({'A','T'},{'C','G'}): raise AssertionError('palindromic variant')
 orientation=1 if a['effect_allele']==b['alt'] else -1
 vals=[float(a[k]) for k in ('fi_beta','fi_se','r13_beta','r13_se')]+[float(b[k]) for k in ('adhd_beta','adhd_se','insomnia_beta','insomnia_se')]
 if not all(math.isfinite(v) for v in vals) or any(vals[i]<=0 for i in (1,3,5,7)): continue
 merged.append((vals,orientation))
assert len(merged)==1496
def logsumexp(xs):
 m=max(xs); return m+math.log(sum(math.exp(x-m) for x in xs))
def logabf(b,s,prior_sd):
 v=s*s; w=prior_sd**2; r=w/(v+w)
 return .5*(math.log1p(-r)+r*(b/s)**2)
a=[r[0] for r in merged]; l1=[logabf(v[0],v[1],.15) for v in a]; errors=[]
rows=read(out/'posteriors.tsv')
for row in rows:
 ix={'ADHD':(4,5),'FinnGen_insomnia':(2,3),'Brain6_insomnia':(6,7)}[row['trait']]
 l2=[logabf(v[ix[0]],v[ix[1]],.2) for v in a]
 s1=logsumexp(l1);s2=logsumexp(l2);shared=logsumexp([x+y for x,y in zip(l1,l2)])
 # Distinct variants: sum exp(l1_i + l2_j) for i != j.
 distinct=s1+s2+math.log1p(-math.exp(shared-s1-s2))
 weights=[0,math.log(1e-4)+s1,math.log(1e-4)+s2,math.log(1e-8)+distinct,math.log(float(row['p12']))+shared]
 post=[math.exp(w-logsumexp(weights)) for w in weights]; errors.extend(abs(post[i]-float(row[f'H{i}'])) for i in range(5))
assert max(errors)<1e-12
receipt={'fi_chr3_rows':len(f),'pilot_rows':len(p),'common_rows':len(merged),'max_absolute_posterior_difference':max(errors),'input_sha256':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in (fi,pilot)},'note':'Independent posterior arithmetic and intersection count; not model validation or external replication.'}
(out/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
