import csv,math,json,hashlib
import numpy as np
from scipy.special import logsumexp
from pathlib import Path
p=Path(__file__).resolve().parent;base=Path('/Volumes/Brain6Work/xqtl-chr5/methylation_candidate')
with (p/'coloc_posteriors.tsv').open() as f:results=list(csv.DictReader(f,delimiter='\t'))
assert len(results)==16
worst=0
for r in results:
 with (base/(r['trait']+'.tsv')).open() as f:d=list(csv.DictReader(f,delimiter='\t'))
 assert len(d)==int(r['n'])
 b1=np.array([float(x['qtl_beta']) for x in d]);s1=np.array([float(x['qtl_se']) for x in d]);b2=np.array([float(x['gwas_beta_alt']) for x in d]);s2=np.array([float(x['gwas_se']) for x in d]);sd=float(r['qtl_prior_sd'])
 w1=sd**2/(sd**2+s1*s1);w2=.2**2/(.2**2+s2*s2);l1=.5*(np.log1p(-w1)+w1*(b1/s1)**2);l2=.5*(np.log1p(-w2)+w2*(b2/s2)**2)
 a=logsumexp(l1);b=logsumexp(l2);c=logsumexp(l1+l2);h3=a+b+math.log1p(-math.exp(c-a-b));ls=np.array([0,math.log(1e-4)+a,math.log(1e-4)+b,math.log(1e-8)+h3,math.log(float(r['p12']))+c]);post=np.exp(ls-logsumexp(ls));worst=max(worst,max(abs(post-np.array([float(r[f'H{i}']) for i in range(5)]))))
assert worst<1e-12
out={'posterior_rows':len(results),'maximum_absolute_difference':worst,'independent_check':'Python direct Wakefield ABF and posterior recomputation; arithmetic only, no model validation','posterior_sha256':hashlib.sha256((p/'coloc_posteriors.tsv').read_bytes()).hexdigest()}
(p/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
