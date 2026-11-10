import csv,math,json,hashlib
import numpy as np
from pathlib import Path
src=Path('outputs/reference_22_score_variance_decomposition.csv')
rows=list(csv.DictReader(src.open()))
assert len(rows)==22 and len({r['variant'] for r in rows})==22
# Verify the covariance identity against direct deletion on synthetic data.
rng=np.random.default_rng(9302); x=rng.normal(size=(500,5)); weights=rng.normal(size=5); s=x@weights
for i in range(5):
 t=x[:,i]*weights[i]; var=np.var(s,ddof=1); a=np.cov(t,s,ddof=1)[0,1]/var; d=np.var(t,ddof=1)/var
 assert np.isclose((1-a)/math.sqrt(1-2*a+d),np.corrcoef(s,s-t)[0,1],atol=1e-12)
out=[]
for name in ('FI','FFS','insomnia'):
 for r in rows:
  a=float(r[name+'_allocated_covariance_fraction']);d=float(r[name+'_standalone_variance_fraction'])
  reduced_variance=1-2*a+d
  assert reduced_variance>0
  correlation=(1-a)/math.sqrt(reduced_variance)
  assert -1.000000001<=correlation<=1.000000001
  out.append({'score':name,'omitted_variant':r['variant'],'correlation_reduced_full':correlation,'reduced_to_full_sd_ratio':math.sqrt(reduced_variance)})
target=Path('outputs/stage6_weighted22_missing_marker_sensitivity.csv')
with target.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
receipt={'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'rows':len(out),'formula':'a=Cov(w_i X_i,S)/Var(S), d=Var(w_i X_i)/Var(S); corr(S-w_i X_i,S)=(1-a)/sqrt(1-2a+d)','limitation':'Reference-only descriptive deletion, no imputation, no phenotype or prediction validation.'}
Path('outputs/stage6_weighted22_missing_marker_sensitivity.json').write_text(json.dumps(receipt,indent=2)+'\n')
for name in ('FI','FFS','insomnia'):
 z=sorted((r for r in out if r['score']==name),key=lambda r:r['correlation_reduced_full'])
 print(name, [(r['omitted_variant'],round(r['correlation_reduced_full'],4),round(r['reduced_to_full_sd_ratio'],4)) for r in z[:3]])
