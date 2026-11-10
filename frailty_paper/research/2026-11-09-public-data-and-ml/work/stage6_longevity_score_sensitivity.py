import csv,json,math,hashlib
from pathlib import Path
import numpy as np
hits=json.loads(Path('outputs/stage6_longevity_candidate_screen.json').read_text())['matched']
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')))
p=np.load('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True);assert list(p['variants'])==[r['variant'] for r in rows];x=p['genotypes'];sd=x.std(axis=0,ddof=1);R=np.corrcoef(x,rowvar=False);chrom=np.array([r['chr'] for r in rows]);R=np.where(chrom[:,None]==chrom[None,:],R,0)
results=[]
for f in ('FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight'):
 h=[next(y for y in hits if y['score']=='15' and y['variant']==r['variant']) for r in rows];a=np.array([float(r[f]) for r in rows])*sd;z=np.array([y['beta_counted']/y['se'] for y in h]);v=float(a@z/np.sqrt(a@R@a));q=np.array([y['se'] for y in h])*sd
 results.append({'prototype':'15','weight':f,'z_approx':v,'nominal_p_approx':math.erfc(abs(v)/math.sqrt(2)),'se_sd_min':float(q.min()),'se_sd_max':float(q.max())})
score=list(csv.DictReader(open('outputs/unvalidated_22_shared_direction_score.csv')));d=json.loads(Path('outputs/shared22_reference_diagnostic.json').read_text());assert hashlib.sha256(Path('outputs/unvalidated_22_shared_direction_score.csv').read_bytes()).hexdigest()==d['shared_score_sha256'];freq={r['variant']:r for r in csv.DictReader(open('outputs/22_score_reference_group_allele_frequencies.csv'))};a=[];z=[];q=[]
for r in score:
 h=next(y for y in hits if y['score']=='22' and y['variant']==r['variant']);v=freq[r['variant']];p0=float(v['reference_EUR_frequency']);p0=p0 if v['counted_allele']==r['counted_allele'] else 1-p0;sd0=math.sqrt(2*p0*(1-p0));a.append(sd0);z.append(h['beta_counted']/h['se']);q.append(h['se']*sd0)
val=float(np.dot(a,z)/d['shared_allele_count_sd']);results.append({'prototype':'22','weight':'equal_count','z_approx':val,'nominal_p_approx':math.erfc(abs(val)/math.sqrt(2)),'se_sd_min':min(q),'se_sd_max':max(q)})
Path('outputs/stage6_longevity_score_sensitivity.json').write_text(json.dumps({'results':results,'assumptions':'Approximate weighted marginal Z using 503 1000G EUR reference. 15 uses empirical genotype SD and within-chromosome correlation; 22 uses HWE SD and empirical whole-count SD including cross-chromosome sample correlations.','limitations':'postselected, target LD and sample overlap unknown, differing variant effective Ns, meta-analysis SEs and reference variance; longevity rather than frailty; not participant-level score prediction or direct comparison'},indent=2)+'\n');print(json.dumps(results,indent=2))
