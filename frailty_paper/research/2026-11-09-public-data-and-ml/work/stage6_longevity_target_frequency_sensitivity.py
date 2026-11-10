import csv,json,math
from pathlib import Path
import numpy as np
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')));hits=json.loads(Path('outputs/stage6_longevity_candidate_screen.json').read_text())['matched'];d=np.load('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True);assert list(d['variants'])==[r['variant'] for r in rows];x=d['genotypes'];chrom=np.array([r['chr'] for r in rows]);R=np.where(chrom[:,None]==chrom[None,:],np.corrcoef(x,rowvar=False),0)
h=[next(v for v in hits if v['score']=='15' and v['variant']==r['variant']) for r in rows];sd=np.array([math.sqrt(2*v['eaf_counted']*(1-v['eaf_counted'])) for v in h]);z=np.array([v['beta_counted']/v['se'] for v in h]);prior=json.loads(Path('outputs/stage6_longevity_score_sensitivity.json').read_text())['results'];out=[]
for field in ['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']:
 a=np.array([float(r[field]) for r in rows])*sd;Z=float(a@z/np.sqrt(a@R@a));out.append({'weight':field,'original_z':next(v['z_approx'] for v in prior if v['weight']==field),'target_eaf_z':Z})
Path('outputs/stage6_longevity_target_frequency_sensitivity.json').write_text(json.dumps({'results':out,'limits':'HWE, reported meta-analysis EAF, external EUR LD and heterogeneous per-marker effective samples are approximations; not individual score validation.'},indent=2)+'\n');print(json.dumps(out,indent=2))
