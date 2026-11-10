import csv,json,math
import numpy as np
from pathlib import Path
rows=list(csv.DictReader(open('outputs/unvalidated_15_variant_score_prototype.csv')));d=json.loads(Path('work/mdd_symptoms/candidate_hits_private.json').read_text());matches={(v['chr'],v['pos'],frozenset(v['pair'])):v for v in d};p=np.load('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz',allow_pickle=True);assert list(p['variants'])==[r['variant'] for r in rows];X=p['genotypes'];assert X.shape==(503,15) and np.all(np.isfinite(X));sd=X.std(axis=0,ddof=1);assert np.all(sd>0);R=np.corrcoef(X,rowvar=False);c=np.array([r['chr'] for r in rows]);R=np.where(c[:,None]==c[None,:],R,0)
out=[]
for trait in ['insomnia','fatigue']:
 z=[];se=[]
 for r in rows:
  k=(r['chr'],int(r['position_GRCh37_inferred']),frozenset((r['counted_allele'],r['other_allele'])));h=matches[k]['hits'][trait][0];v=float(h['se']);assert v>0;se.append(v);z.append(float(h['beta'])/v*(1 if h['a1']==r['counted_allele'] else -1))
 se=np.array(se);z=np.array(z);q=se*sd
 for f in ['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']:
  a=np.array([float(r[f]) for r in rows])*sd;Z=float(a@z/np.sqrt(a@R@a));out.append({'trait':trait,'weight':f,'n':15,'z':Z,'nominal_p':math.erfc(abs(Z)/math.sqrt(2)),'se_times_reference_sd_min':float(q.min()),'se_times_reference_sd_max':float(q.max()),'se_times_reference_sd_median':float(np.median(q))})
Path('outputs/stage6_mdd_15_weighted_sensitivity.json').write_text(json.dumps({'results':out,'reference':'previously cached 1000 Genomes EUR 503 scored-allele genotypes, per prior receipt','formula':'a = score weights * reference genotype SD; Z=a.T z / sqrt(a.T R a)','caveats':'Exploratory summary approximation, heterogeneous case-enriched phenotype and sample sizes; reference rather than target genotype variance and LD; not score prediction or independent replication'},indent=2)+'\n');print(json.dumps(out,indent=2))
