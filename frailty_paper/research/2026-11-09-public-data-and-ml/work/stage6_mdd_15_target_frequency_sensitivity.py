import csv,json,math
from pathlib import Path
import numpy as np
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
hits=json.loads(Path('work/mdd_symptoms/candidate_hits_private.json').read_text())
lookup={(x['chr'],x['pos'],frozenset(x['pair'])):x for x in hits}
d=np.load('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True)
assert list(d['variants'])==[r['variant'] for r in rows]
x=d['genotypes']; R=np.corrcoef(x,rowvar=False); chrom=np.array([r['chr'] for r in rows]); R=np.where(chrom[:,None]==chrom[None,:],R,0)
old=json.loads(Path('outputs/stage6_mdd_15_weighted_sensitivity.json').read_text())['results']
results=[]
for trait in ('insomnia','fatigue'):
 z=[]; sd=[]
 for r in rows:
  k=(r['chr'],int(r['position_GRCh37_inferred']),frozenset((r['counted_allele'],r['other_allele'])))
  h=lookup[k]['hits'][trait][0]; n1=float(h['ncase']); n0=float(h['ncontrol']); p=(n1*float(h['fcas'])+n0*float(h['fcon']))/(n1+n0)
  assert 0<p<1 and float(h['se'])>0
  z.append(float(h['beta'])/float(h['se'])*(1 if h['a1']==r['counted_allele'] else -1)); sd.append(math.sqrt(2*p*(1-p)))
 for field in ('FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'):
  a=np.array([float(r[field]) for r in rows])*np.array(sd); value=float(a@z/np.sqrt(a@R@a)); prior=next(o['z'] for o in old if o['trait']==trait and o['weight']==field)
  results.append({'trait':trait,'weight':field,'target_frequency_approx_z':value,'reference_sd_approx_z':prior,'nominal_p':math.erfc(abs(value)/math.sqrt(2))})
Path('outputs/stage6_mdd_15_target_frequency_sensitivity.json').write_text(json.dumps({'results':results,'method':'Use case-count weighted reported case/control A1 frequencies to estimate genotype SD under HWE; retain EUR503 reference within-chromosome correlation and zero cross-chromosome correlations. Compare to empirical EUR503 genotype SD.','limitations':'Case/control frequencies in a meta-analysis may not be exactly pooled, HWE and common target covariance are approximations, estimates have different marker sample sizes and possible cohort overlap; no individual-level score validation.'},indent=2)+'\n')
print(json.dumps(results,indent=2))
