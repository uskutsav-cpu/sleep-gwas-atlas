#!/usr/bin/env python3
"""Postselected marginal-Z comparison, conditional on complete verified inputs."""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy.stats import norm
O=Path('outputs')
receipt=json.loads((O/'stage6_grip_gwas_extract_receipt.json').read_text())
assert receipt['source_checksum_matches'] is True
score=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv')
g=pd.read_csv(O/'stage6_grip_gwas_candidate_aligned.csv')
g=g[g.score.eq('prototype15')]
assert len(score)==len(g)==15 and score.variant.is_unique and g.variant.is_unique
x=score.merge(g,on='variant',suffixes=('_score','_target'),validate='one_to_one')
assert len(x)==15 and (x.exact_rows==1).all()
assert (x.counted_allele_score==x.counted_allele_target).all() and (x.other_allele_score==x.other_allele_target).all()
for c in ['beta_counted','se','p','counted_allele_frequency','FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']:
 assert np.isfinite(x[c]).all(),c
assert (x.se>0).all() and x.p.between(0,1).all() and x.counted_allele_frequency.between(0,1,inclusive='neither').all()
ld=pd.read_csv(O/'stage6_hfrs_sparse_score_EUR503_ld.csv',index_col=0)
assert set(ld.index)==set(x.variant)==set(ld.columns)
R=ld.loc[x.variant,x.variant].to_numpy(dtype=float)
assert np.allclose(R,R.T) and np.allclose(np.diag(R),1) and np.linalg.eigvalsh(R).min()>-1e-8
z=x.beta_counted.to_numpy()/x.se.to_numpy()
sd=np.sqrt(2*x.counted_allele_frequency.to_numpy()*(1-x.counted_allele_frequency.to_numpy()))
results=[]
for weight in ['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']:
 w=x[weight].to_numpy()*sd
 for omitted in [None,*x.variant]:
  keep=np.ones(len(x),bool) if omitted is None else x.variant.ne(omitted).to_numpy()
  denom=float(np.sqrt(w[keep]@R[np.ix_(keep,keep)]@w[keep]));assert denom>0
  Z=float(w[keep]@z[keep]/denom)
  results.append({'weight':weight,'omitted_variant':omitted or 'none','n':int(keep.sum()),'z_approx':Z,'two_sided_p_approx':float(2*norm.sf(abs(Z)))})
out=O/'stage6_grip_sparse_score_summary.csv';pd.DataFrame(results).to_csv(out,index=False)
meta={'source':receipt['source'],'source_md5':receipt['observed_md5'],'reference_ld':'1000 Genomes phase 3 EUR503; zero cross-chromosome LD as previously constructed','score_sha256':hashlib.sha256((O/'unvalidated_15_variant_score_prototype.csv').read_bytes()).hexdigest(),'caveats':['Exploratory postselected weighted marginal Z, not individual-level prediction or calibration.','UK Biobank overlaps score discovery; low grip strength is a component of physical frailty, not FI or FFS.','Reference LD approximates null Z correlation, requiring compatible study samples, effects and errors; per-variant sample and heterogeneity may differ across cohorts.','Hardy-Weinberg genotype standard deviation uses target effect-allele frequency; score weights and FFS orientation remain research prototypes.']}
(O/'stage6_grip_sparse_score_summary_receipt.json').write_text(json.dumps(meta,indent=2)+'\n')
print(pd.DataFrame(results).query("omitted_variant == 'none'").to_string(index=False))
