import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm
O=Path('outputs');p=O/'stage6_aou_eur_weighted_z_components.csv';d=pd.read_csv(p)
left=d[(d.score=='15_weighted')&(d.weights=='FI')&(d.outcome=='NS_333.2')].set_index('variant'); right=d[(d.score=='15_weighted')&(d.weights=='FFS')&(d.outcome=='NS_333.2')].set_index('variant').loc[left.index]
assert len(left)==15 and np.allclose(left.conditional_counted_z,right.conditional_counted_z)
R=pd.read_csv(O/'stage6_hfrs_sparse_score_EUR503_ld.csv',index_col=0).loc[left.index,left.index].to_numpy();weights=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv').set_index('variant').loc[left.index]
a=left.weight_on_counted_allele.to_numpy()*left.hwe_genotype_sd.to_numpy();b=right.weight_on_counted_allele.to_numpy()*right.hwe_genotype_sd.to_numpy();z=left.conditional_counted_z.to_numpy();out=[]
masks=[('none',np.ones(15,dtype=bool))]+[(i,np.asarray(left.index!=i)) for i in left.index]
for chrom in sorted(weights.chr.unique()):
 keep=np.asarray(weights.chr!=chrom)
 if keep.sum()<15:masks.append(('chromosome_'+str(chrom),keep))
for dropped,keep in masks:
 x=a[keep];y=b[keep];r=R[np.ix_(keep,keep)];v=z[keep];x=x/np.sqrt(x@r@x);y=y/np.sqrt(y@r@y);rho=float(x@r@y);var=float((x-y)@r@(x-y));assert var>0 and np.isclose(var,2-2*rho)
 za=float(x@v);zb=float(y@v);contrast=(za-zb)/np.sqrt(var)
 out.append({'omitted':dropped,'retained':int(keep.sum()),'FI_z':za,'FFS_z':zb,'assumed_correlation':rho,'contrast_z':contrast,'contrast_nominal_p':float(2*norm.sf(abs(contrast)))})
pd.DataFrame(out).to_csv(O/'stage6_aou_eur_paired_contrast_omission.csv',index=False)
receipt={'component_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'ld_sha256':hashlib.sha256((O/'stage6_hfrs_sparse_score_EUR503_ld.csv').read_bytes()).hexdigest(),'output_sha256':hashlib.sha256((O/'stage6_aou_eur_paired_contrast_omission.csv').read_bytes()).hexdigest(),'limitations':'Exploratory post hoc reference-LD-based sensitivity; changing marker set creates a different score and does not validate a prediction model.'}
(O/'stage6_aou_eur_paired_contrast_omission_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(pd.DataFrame(out).to_string(index=False))
