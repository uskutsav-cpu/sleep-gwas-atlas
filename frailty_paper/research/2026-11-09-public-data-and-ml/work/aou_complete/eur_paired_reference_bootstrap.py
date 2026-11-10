import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
O=Path('outputs');inp=O/'stage6_aou_eur_weighted_z_components.csv';gen=Path('work/stage6_hfrs_sparse_score_EUR503_genotypes.npz');sc=O/'unvalidated_15_variant_score_prototype.csv'
d=pd.read_csv(sc);comp=pd.read_csv(inp);left=comp[(comp.score=='15_weighted')&(comp.weights=='FI')&(comp.outcome=='NS_333.2')].set_index('variant').loc[d.variant];right=comp[(comp.score=='15_weighted')&(comp.weights=='FFS')&(comp.outcome=='NS_333.2')].set_index('variant').loc[d.variant];assert np.allclose(left.conditional_counted_z,right.conditional_counted_z)
with np.load(gen,allow_pickle=True) as f:G=np.asarray(f['genotypes'],float);ids=f['variants'].astype(str)
assert G.shape==(503,15) and list(ids)==list(d.variant) and np.isfinite(G).all()
chrom=d.chr.astype(str).to_numpy();block=chrom[:,None]==chrom[None,:];A=left.weight_on_counted_allele.to_numpy()*left.hwe_genotype_sd.to_numpy();B=right.weight_on_counted_allele.to_numpy()*right.hwe_genotype_sd.to_numpy();z=left.conditional_counted_z.to_numpy()
def calc(R,keep):
 r=R[np.ix_(keep,keep)];a=A[keep];b=B[keep];zz=z[keep];a=a/np.sqrt(a@r@a);b=b/np.sqrt(b@r@b);den=(a-b)@r@(a-b)
 if den<=0:raise ValueError('invalid denominator')
 return float((a@zz-b@zz)/np.sqrt(den))
R=np.where(block,np.corrcoef(G,rowvar=False),0);saved=pd.read_csv(O/'stage6_hfrs_sparse_score_EUR503_ld.csv',index_col=0).loc[d.variant,d.variant].to_numpy();assert np.allclose(R,saved,atol=1e-10)
keeps={'all':np.ones(15,dtype=bool),'without_rs113851554':np.asarray(d.variant!='rs113851554')};point={name:calc(R,k) for name,k in keeps.items()};prev=pd.read_csv(O/'stage6_aou_eur_paired_contrast_omission.csv').set_index('omitted');assert np.isclose(point['all'],prev.loc['none','contrast_z']) and np.isclose(point['without_rs113851554'],prev.loc['rs113851554','contrast_z'])
rng=np.random.default_rng(20261029);vals={name:[] for name in keeps};discard=0
for i in range(5000):
 B_=G[rng.integers(0,len(G),len(G))]
 if np.any(B_.std(axis=0)==0):discard+=1;continue
 r=np.where(block,np.corrcoef(B_,rowvar=False),0)
 for name,k in keeps.items():vals[name].append(calc(r,k))
out=[]
for name,v in vals.items():
 q=np.quantile(v,[.025,.5,.975]);out.append({'scenario':name,'fixed_reference_z':point[name],'bootstrap_2_5':q[0],'bootstrap_median':q[1],'bootstrap_97_5':q[2],'valid_draws':len(v),'discarded':discard})
output=O/'stage6_aou_paired_reference_bootstrap.csv';pd.DataFrame(out).to_csv(output,index=False)
receipt={'seed':20261029,'draws_requested':5000,'genotypes_sha256':hashlib.sha256(gen.read_bytes()).hexdigest(),'components_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'output_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),'limits':'Reference individuals resampled with target GWAS and weights fixed; percentiles are not confidence intervals for target associations, effects, prediction, or sample overlap. Reference population may be structured or related; bootstrap does not address reference mismatch.'};(O/'stage6_aou_paired_reference_bootstrap_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(pd.DataFrame(out).to_string(index=False))
