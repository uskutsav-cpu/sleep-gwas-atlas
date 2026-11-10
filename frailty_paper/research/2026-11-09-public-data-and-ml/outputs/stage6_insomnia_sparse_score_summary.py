#!/usr/bin/env python3
"""Post hoc weighted-Z summary association, not individual score validation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.stats import norm
O=Path('outputs'); score=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv');lookup=pd.read_csv(O/'prototype_15_variants_full_external_insomnia_lookup.csv');out=[];details=[]
for source in ['MVP','R13']:
 d=score.merge(lookup[lookup.source.eq(source)],left_on='variant',right_on='variant',validate='one_to_one');assert len(d)==15
 assert (d.counted_allele==d.candidate_effect_allele).all() and (d.other_allele==d.candidate_other_allele).all();assert d.allele_relation.isin(['direct','swap']).all()
 f=np.where(d.source_effect_allele==d.counted_allele,d.source_EAF,1-d.source_EAF);assert np.isfinite(f).all() and ((f>0)&(f<1)).all()
 if source=='MVP':
  p=d.source_p.to_numpy();beta=d.aligned_beta_to_candidate_allele.to_numpy();assert (p>0).all() and (p<1).all() and np.isfinite(beta).all() and (beta!=0).all();z=np.sign(beta)*norm.isf(p/2)
 else:
  se=pd.to_numeric(d.source_SE,errors='raise').to_numpy();assert (se>0).all();z=d.aligned_beta_to_candidate_allele.to_numpy()/se
 for pop in ['EUR503','FIN99']:
  ld=pd.read_csv(O/f'stage6_hfrs_sparse_score_{pop}_ld.csv',index_col=0);R=ld.loc[d.variant,d.variant].to_numpy();assert np.allclose(R,R.T) and np.linalg.eigvalsh(R).min()>-1e-8
  for weight in ['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']:
   a=d[weight].to_numpy()*np.sqrt(2*f*(1-f));
   for omit in [None,*d.variant]:
    keep=np.ones(len(d),dtype=bool) if omit is None else (d.variant!=omit).to_numpy();denom=float(np.sqrt(a[keep]@R[np.ix_(keep,keep)]@a[keep]));stat=float(a[keep]@z[keep]/denom);out.append({'source':source,'reference':pop,'weight':weight,'omitted':omit or 'none','n':int(sum(keep)),'weighted_z':stat,'approx_two_sided_p':2*norm.sf(abs(stat))})
   if pop==('EUR503' if source=='MVP' else 'FIN99'):
    for i,v in enumerate(d.variant):details.append({'source':source,'weight':weight,'variant':v,'z_counted_allele':z[i],'frequency_counted_allele':f[i],'weight_counted_allele':d[weight].iloc[i],'weighted_numerator':a[i]*z[i]})
pd.DataFrame(out).to_csv(O/'stage6_insomnia_sparse_score_summary.csv',index=False);pd.DataFrame(details).to_csv(O/'stage6_insomnia_sparse_score_components.csv',index=False);receipt={'weight_sha256':hashlib.sha256((O/'unvalidated_15_variant_score_prototype.csv').read_bytes()).hexdigest(),'lookup_sha256':hashlib.sha256((O/'prototype_15_variants_full_external_insomnia_lookup.csv').read_bytes()).hexdigest(),'MVP_z':'sign of aligned log odds ratio times two-sided normal quantile from reported p; odds ratios rounded; no usable SE','R13_z':'aligned beta / source SE','frequency':'source effect-allele frequency oriented to prototype counted allele','reference':'same previously verified ordered 1000 Genomes dosage LD matrices; cross-chromosomal correlations set to zero','limitations':'Post hoc, weights and candidates previously inspected against target data, target overlap and variant interpretation uncertain. Approximate marginal-Z covariance and common effective sample; not individual-level score prediction or independent replication.'};(O/'stage6_insomnia_sparse_score_summary_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(pd.DataFrame(out).query("omitted == 'none'").to_string(index=False))
