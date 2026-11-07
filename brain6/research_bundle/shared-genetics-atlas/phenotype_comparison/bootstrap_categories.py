#!/usr/bin/env python3
import csv,json,hashlib,pathlib,numpy as np
p=pathlib.Path('outputs/shared-genetics-atlas/phenotype_comparison');source=p/'ukb_cross_tab.csv'
with source.open() as f:rows=list(csv.DictReader(f))
fp=['robust','prefrail','frail'];fi=['robust','mild','moderate_severe'];counts=np.zeros((3,3),dtype=np.int64);seen=set()
for x in rows:
 i=fp.index(x['fried_category']);j=fi.index(x['fi_category']);assert (i,j) not in seen;seen.add((i,j));counts[i,j]=int(x['n'])
assert len(seen)==9 and counts.sum()==383845 and np.all(counts>=0)
def rho(n):
 n=np.asarray(n,dtype=float);p=n/n.sum();r=p.sum(axis=1);c=p.sum(axis=0);rx=np.cumsum(r)-r/2;ry=np.cumsum(c)-c/2;mx=(r*rx).sum();my=(c*ry).sum();vx=(r*(rx-mx)**2).sum();vy=(c*(ry-my)**2).sum();return ((p*(rx-mx)[:,None]*(ry-my)[None,:]).sum())/np.sqrt(vx*vy)
obs=float(rho(counts));assert abs(obs-.28447194444596086)<1e-12
rng=np.random.default_rng(20261027);samples=rng.multinomial(int(counts.sum()),(counts/counts.sum()).ravel(),size=10000);boot=np.array([rho(z.reshape(3,3)) for z in samples]);assert np.isfinite(boot).all();lo,hi=np.quantile(boot,[.025,.975]);
with (p/'bootstrap_rho.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['replicate','rho']);w.writerows(enumerate(boot,1))
meta={'n':int(counts.sum()),'observed_category_spearman':obs,'n_replicates':len(boot),'seed':20261027,'rng':'NumPy default_rng PCG64','numpy':np.__version__,'bootstrap_percentile_95CI':[float(lo),float(hi)],'bootstrap_se':float(np.std(boot,ddof=1)),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
(p/'bootstrap_summary.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
