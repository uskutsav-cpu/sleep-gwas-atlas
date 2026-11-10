import csv, hashlib,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import norm
O=Path('outputs');hits=['rs113851554','rs12536335'];out=[]
for pop in ('eur','afr'):
 comp=pd.read_csv(O/f'stage6_aou_{pop}_weighted_z_components.csv');old=pd.read_csv(O/f'stage6_aou_{pop}_weighted_z_sensitivity.csv')
 for score in ('15_weighted','22_count'):
  ldname={'eur':{'15_weighted':'stage6_hfrs_sparse_score_EUR503_ld.csv','22_count':'stage6_aou_count22_EUR503_ld.csv'},'afr':{'15_weighted':'stage6_aou_weighted15_AFR661_ld.csv','22_count':'stage6_aou_count22_AFR661_ld.csv'}}[pop][score]
  ld=pd.read_csv(O/ldname,index_col=0)
  for wt in (['FI','FFS','insomnia'] if score=='15_weighted' else ['count','FI','FFS','insomnia']):
   c=comp[(comp.score==score)&(comp.weights==wt)&(comp.outcome=='NS_333.2')].copy();assert len(c)==(15 if score=='15_weighted' else 22)
   R=ld.loc[c.variant,c.variant].to_numpy();a=(c.hwe_genotype_sd*c.weight_on_counted_allele).to_numpy();z=c.conditional_counted_z.to_numpy();orig=float(a@z/np.sqrt(a@R@a))
   ref=old[(old.score==score)&(old.weights==wt)&(old.outcome=='NS_333.2')];assert len(ref)==1 and np.isclose(orig,float(ref.iloc[0].weighted_z_conditional),rtol=1e-12,atol=1e-12)
   targets=[[]]+[[h] for h in hits if h in set(c.variant)]
   if set(hits).issubset(set(c.variant)):targets+=[hits]
   for omit in targets:
    keep=~c.variant.isin(omit).to_numpy();d=float(a[keep]@R[np.ix_(keep,keep)]@a[keep]);assert d>0
    v=float(a[keep]@z[keep]/np.sqrt(d))
    out.append({'population':pop.upper(),'score':score,'weights':wt,'omitted':'|'.join(omit) if omit else 'none','n_retained':int(keep.sum()),'conditional_z':v,'nominal_two_sided_p':float(2*norm.sf(abs(v)))})
p=O/'stage6_aou_insomnia_two_locus_removal.csv';pd.DataFrame(out).to_csv(p,index=False)
receipt={'rows':len(out),'output_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'method':'Post hoc complete versus selected one or two loci omitted; variant weights, marginal z and population-matched reference LD inherited unchanged','limitations':'Changing the marker set changes the score and its standardization; these are dependent exploratory quantities, not predictive comparison, score improvement or validation. Effect-allele, ancestry LD and discovery-overlap assumptions persist.'}
(O/'stage6_aou_insomnia_two_locus_removal_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(pd.DataFrame(out).query("population=='EUR'").to_string(index=False))
