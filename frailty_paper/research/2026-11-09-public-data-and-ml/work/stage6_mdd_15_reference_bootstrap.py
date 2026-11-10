import csv,json,math
from pathlib import Path
import numpy as np
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
hits=json.loads(Path('work/mdd_symptoms/candidate_hits_private.json').read_text())
lookup={(h['chr'],h['pos'],frozenset(h['pair'])):h for h in hits}
data=np.load('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True)
assert list(data['variants'])==[r['variant'] for r in rows]
x=data['genotypes']; assert x.shape==(503,15)
chrom=np.array([r['chr'] for r in rows]); same=chrom[:,None]==chrom[None,:]
fields=['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']
w=np.array([[float(r[f]) for r in rows] for f in fields]); zs=[]
for trait in ['insomnia','fatigue']:
 z=[]
 for r in rows:
  h=lookup[(r['chr'],int(r['position_GRCh37_inferred']),frozenset((r['counted_allele'],r['other_allele'])))]['hits'][trait][0]
  z.append(float(h['beta'])/float(h['se'])*(1 if h['a1']==r['counted_allele'] else -1))
 zs.append(z)
zs=np.array(zs);rng=np.random.default_rng(20261028); values=[]; rejected=0
for i in range(2000):
 sample=x[rng.integers(0,503,503)]; sd=sample.std(axis=0,ddof=1)
 if np.any(sd==0): rejected+=1;continue
 R=np.where(same,np.corrcoef(sample,rowvar=False),0); a=w*sd; denominator=np.sqrt(np.einsum('ij,jk,ik->i',a,R,a)); values.append((zs@a.T)/denominator)
values=np.array(values); out=[]
for ti,t in enumerate(['insomnia','fatigue']):
 for wi,f in enumerate(fields):
  v=values[:,ti,wi]; out.append({'trait':t,'weight':f,'median':float(np.median(v)),'reference_resampling_2_5_percentile':float(np.quantile(v,.025)),'reference_resampling_97_5_percentile':float(np.quantile(v,.975))})
Path('outputs/stage6_mdd_15_reference_bootstrap.json').write_text(json.dumps({'seed':20261028,'requested':2000,'used':len(values),'zero_variance_resamples_excluded':rejected,'results':out,'interpretation':'Descriptive sensitivity to resampling reference individuals, conditional on fixed GWAS statistics and other assumptions; not a confidence interval for an effect or association.'},indent=2)+'\n');print(json.dumps(out,indent=2));print('excluded',rejected)
