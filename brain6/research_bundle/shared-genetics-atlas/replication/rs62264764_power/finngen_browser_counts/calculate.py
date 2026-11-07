import csv,json,math,hashlib
from pathlib import Path
from statistics import NormalDist
out=Path(__file__).resolve().parent
src=out.parents[2]/'variant_context/rs62264764_dense_qtl/multisignal_pilot/variants.tsv'
with src.open() as f: matches=[r for r in csv.DictReader(f,delimiter='\t') if r['rsid']=='rs62264764']
assert len(matches)==1
beta=abs(float(matches[0]['adhd_beta']));freq=float(matches[0]['reference_alt_frequency']);cases=5559;controls=489493;neff=4*cases*controls/(cases+controls);N=NormalDist();var=2*freq*(1-freq)
def power(scale,alpha):
 z=N.inv_cdf(1-alpha/2);mu=scale*beta*math.sqrt(neff*var/4);return N.cdf(mu-z)+N.cdf(-z-mu)
rows=[{'scale':scale,'alpha':alpha,'power':power(scale,alpha)} for scale in [0,.5,1,1.5] for alpha in [.05,.05/6]]
with (out/'scenarios.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
receipt={'phenotype':'F5_ADHD','release':'FinnGen R13','gw_browser_case_count':cases,'gw_browser_control_count':controls,'browser_url':'https://r13.finngen.fi/','browser_accessed':'2026-10-27','n_eff':neff,'discovery_abs_beta':beta,'reference_alt_frequency':freq,'variant_source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scenario_sha256':hashlib.sha256((out/'scenarios.tsv').read_bytes()).hexdigest(),'warning':'count-based illustration, not variant-level validation, actual Finnish allele frequency and association information unknown'}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));print(rows)
