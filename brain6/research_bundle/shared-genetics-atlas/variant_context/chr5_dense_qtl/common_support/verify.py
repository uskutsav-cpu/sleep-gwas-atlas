import csv,hashlib,json,math,pathlib
p=pathlib.Path(__file__).resolve().parent
files=[p.parent/'Testis_ENSG00000251574_insomnia_matched.tsv',p.parent/'Testis_ENSG00000251574_adhd_matched.tsv']
def read(path):
 with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
a,b=({r['variant']:r for r in read(q)} for q in files); common=sorted(set(a)&set(b));assert len(common)==1633
for key in common:
 for field in ('ref','alt','qtl_beta','qtl_se'):assert a[key][field]==b[key][field]
def bf(x,s,sd):
 x=float(x);s=float(s);r=sd*sd/(sd*sd+s*s);return .5*(math.log1p(-r)+r*(x/s)**2)
def lsum(x):
 m=max(x);return m+math.log(math.fsum(math.exp(v-m) for v in x))
q=[bf(a[k]['qtl_beta'],a[k]['qtl_se'],.15) for k in common];i=[bf(a[k]['gwas_beta_alt'],a[k]['gwas_se'],.2) for k in common];j=[bf(b[k]['gwas_beta_alt'],b[k]['gwas_se'],.2) for k in common]; pairs={'qtl_insomnia':(q,i),'qtl_adhd':(q,j),'insomnia_adhd':(i,j)}
results=read(p/'pair_posteriors.tsv');snps=read(p/'conditional_snps.tsv');maxerr=0
for row in results:
 x,y=pairs[row['pair']];sx=lsum(x);sy=lsum(y);shared=lsum([u+v for u,v in zip(x,y)]);distinct=sx+sy+math.log1p(-math.exp(shared-sx-sy));w=[0,math.log(1e-4)+sx,math.log(1e-4)+sy,math.log(1e-8)+distinct,math.log(float(row['p12']))+shared];z=lsum(w);prob=[math.exp(v-z) for v in w]
 maxerr=max(maxerr,max(abs(prob[k]-float(row[f'H{k}'])) for k in range(5)))
assert maxerr<1e-10
sets={}
for name,(x,y) in pairs.items():
 z=lsum([u+v for u,v in zip(x,y)]); expected={k:math.exp(u+v-z) for k,u,v in zip(common,x,y)};r=[s for s in snps if s['pair']==name];assert len(r)==len(common)
 assert max(abs(expected[s['variant']]-float(s['conditional_h4'])) for s in r)<1e-10
 assert abs(sum(float(s['conditional_h4']) for s in r)-1)<1e-10
 cum=0;cs=set()
 for s in r:
  keep=cum<.95;assert keep==(s['credible95']=='TRUE')
  if keep:cs.add(s['variant'])
  cum+=float(s['conditional_h4'])
 sets[name]=cs
receipt={'common_count':len(common),'max_posterior_difference':maxerr,'conditional_credible_set_sizes':{k:len(v) for k,v in sets.items()},'intersection_count':len(set.intersection(*sets.values())),'sha256':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in files},'note':'Independent arithmetic and support verification, not statistical validation.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
