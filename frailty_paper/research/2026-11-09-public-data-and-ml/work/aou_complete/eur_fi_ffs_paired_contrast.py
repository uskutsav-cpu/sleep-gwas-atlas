import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm
O=Path('outputs');src=O/'stage6_aou_eur_weighted_z_components.csv';d=pd.read_csv(src);out=[]
for score,ldfile in [('15_weighted',O/'stage6_hfrs_sparse_score_EUR503_ld.csv'),('22_count',O/'stage6_aou_count22_EUR503_ld.csv')]:
 for ph in ['NS_333.2','SS_807.21']:
  left=d[(d.score==score)&(d.weights=='FI')&(d.outcome==ph)].set_index('variant');right=d[(d.score==score)&(d.weights=='FFS')&(d.outcome==ph)].set_index('variant').loc[left.index]
  assert len(left)==(15 if score=='15_weighted' else 22)
  assert np.allclose(left.conditional_counted_z,right.conditional_counted_z) and np.allclose(left.hwe_genotype_sd,right.hwe_genotype_sd)
  R=pd.read_csv(ldfile,index_col=0).loc[left.index,left.index].to_numpy();assert np.allclose(R,R.T) and np.linalg.eigvalsh(R).min()>-1e-8
  a=left.weight_on_counted_allele.to_numpy()*left.hwe_genotype_sd.to_numpy();b=right.weight_on_counted_allele.to_numpy()*right.hwe_genotype_sd.to_numpy();z=left.conditional_counted_z.to_numpy();a/=np.sqrt(a@R@a);b/=np.sqrt(b@R@b)
  r=float(a@R@b);v=float((a-b)@R@(a-b));assert v>0 and np.isclose(v,2-2*r)
  za=float(a@z);zb=float(b@z);zd=(za-zb)/np.sqrt(v)
  out.append({'score_set':'22_weighted' if score=='22_count' else score,'outcome':ph,'n_markers':len(a),'FI_weighted_z':za,'FFS_weighted_z':zb,'assumed_score_correlation':r,'contrast_z':zd,'contrast_nominal_p':float(2*norm.sf(abs(zd)))})
p=O/'stage6_aou_eur_fi_ffs_paired_contrast.csv';pd.DataFrame(out).to_csv(p,index=False)
receipt={'component_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'ld_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [O/'stage6_hfrs_sparse_score_EUR503_ld.csv',O/'stage6_aou_count22_EUR503_ld.csv']},'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'method':'Normalize each score weight vector a = weight*sqrt(2f(1-f)) by sqrt(aRa); compare their projected marginal Z statistics as (Za-Zb)/sqrt(2-2aRb).','limitations':'Exploratory same-outcome Z contrast under common effective N and matched LD; is not a patient-level association or comparison of prediction, FI vs FFS phenotypes, or across-outcome causality. The score weights and candidate/outcome selection are post hoc. FFS weight direction and source alleles require verification; outcome and discovery sample overlap unknown.'}
(O/'stage6_aou_eur_fi_ffs_paired_contrast_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(pd.DataFrame(out).to_string(index=False))
