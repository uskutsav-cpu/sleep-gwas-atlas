import csv,math,json,hashlib
import numpy as np
from scipy.stats import norm
from scipy.special import logsumexp
from pathlib import Path
p=Path(__file__).resolve().parent;parent=p.parent;base=Path('/Volumes/Brain6Work/xqtl-chr5/methylation_candidate')
results=[];integrity={}
for trait in ['insomnia','adhd']:
 path=base/(trait+'.tsv');integrity[trait]=hashlib.sha256(path.read_bytes()).hexdigest()
 with path.open() as f:d=list(csv.DictReader(f,delimiter='\t'))
 b=np.array([float(x['qtl_beta']) for x in d]);se=np.array([float(x['qtl_se']) for x in d]);q=np.array([float(x['qtl_p']) for x in d]);bg=np.array([float(x['gwas_beta_alt']) for x in d]);sg=np.array([float(x['gwas_se']) for x in d]);assert np.all((q>0)&(q<=1))
 zraw=b/se;zequiv=np.sign(zraw)*norm.isf(q/2);ratio=.2**2/(.2**2+sg**2);lg=.5*(np.log1p(-ratio)+ratio*(bg/sg)**2)
 for sd in [.001,.005,.01,.02]:
  r=sd**2/(sd**2+se**2)
  for mode,z in [('original',zraw),('p_value_equivalent_normal',zequiv)]:
   lq=.5*(np.log1p(-r)+r*z*z);a=logsumexp(lq);c=logsumexp(lg);both=logsumexp(lq+lg);logdistinct=a+c+math.log1p(-math.exp(both-a-c))
   for prior in [1e-5,1e-6]:
    logs=np.array([0,math.log(1e-4)+a,math.log(1e-4)+c,math.log(1e-8)+logdistinct,math.log(prior)+both]);post=np.exp(logs-logsumexp(logs))
    results.append({'trait':trait,'n':len(d),'qtl_prior_sd':sd,'p12':prior,'mode':mode,**{f'H{i}':post[i] for i in range(5)}})
with (p/'comparison.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
receipt={'source_sha256':integrity,'rows':len(results),'output_sha256':hashlib.sha256((p/'comparison.tsv').read_bytes()).hexdigest(),'warning':'heuristic p-matched z replacement with unchanged SE and effect prior; not a coherent t-likelihood or calibrated Bayesian analysis'}
(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
for trait in ['insomnia','adhd']:
 for mode in ['original','p_value_equivalent_normal']:
  a=[r['H4'] for r in results if r['trait']==trait and r['mode']==mode];print(trait,mode,min(a),max(a))
