from pathlib import Path
import csv, hashlib, json, math
from statistics import NormalDist
root=Path('outputs/shared-genetics-atlas'); out=root/'replication/rs62264764_power'
path=root/'variant_context/rs62264764_dense_qtl/multisignal_pilot/variants.tsv'
with path.open(newline='') as f: matches=[x for x in csv.DictReader(f,delimiter='\t') if x['rsid']=='rs62264764']
assert len(matches)==1
v=matches[0]; f=float(v['reference_alt_frequency']); beta=abs(float(v['adhd_beta'])); norm=NormalDist(); genotype_var=2*f*(1-f)
def power(n_eff,b,alpha):
 z=norm.inv_cdf(1-alpha/2); mu=abs(b)*math.sqrt(n_eff*genotype_var/4)
 return norm.cdf(mu-z)+norm.cdf(-z-mu)
def solve(alpha):
 lo,hi=0.,1e7
 for _ in range(100):
  mid=(lo+hi)/2
  if power(mid,beta,alpha)<.8: lo=mid
  else: hi=mid
 return hi
rows=[]
for controls in (100000,None):
 cases=5559; n_eff=4*cases if controls is None else 4*cases*controls/(cases+controls)
 for alpha in (.05,.05/6):
  for scale in (0,.5,1,1.5):
   rows.append(dict(cases=cases,controls='infinite' if controls is None else controls,n_eff=n_eff,alpha=alpha,beta_scale=scale,beta=beta*scale,power=power(n_eff,beta*scale,alpha)))
with (out/'power_scenarios.tsv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
target=[{'alpha':alpha,'n_eff_80pct':solve(alpha),'cases_if_infinite_controls_80pct':solve(alpha)/4} for alpha in (.05,.05/6)]
with (out/'target_samples.tsv').open('w',newline='') as h:
 w=csv.DictWriter(h,fieldnames=target[0].keys(),delimiter='\t');w.writeheader();w.writerows(target)
receipt={'variant':'rs62264764','effect':'absolute ADHD ALT log odds','discovery_beta':beta,'reference_alt_frequency':f,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'model':'approximate additive logistic Wald; hypothetical samples','cases_count_status':'endpoint descriptive figure used as hypothetical case count, not verified GWAS cases','two_sided_alpha_scenarios':[.05,.05/6]}
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'rows':rows,'targets':target},indent=2))
