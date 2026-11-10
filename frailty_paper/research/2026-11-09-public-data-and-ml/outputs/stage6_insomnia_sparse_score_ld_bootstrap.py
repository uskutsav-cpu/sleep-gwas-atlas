#!/usr/bin/env python3
"""Reference-only bootstrap of existing post hoc insomnia summary statistics."""
from pathlib import Path
import numpy as np,pandas as pd,json
O=Path('outputs');W=Path('work');d=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv');comp=pd.read_csv(O/'stage6_insomnia_sparse_score_components.csv');observed=pd.read_csv(O/'stage6_insomnia_sparse_score_summary.csv');rng=np.random.default_rng(61028);B=5000;records=[];summary=[];chrom=d['chr'].astype(str).to_numpy();mask=chrom[:,None]==chrom[None,:];weights=['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']
for src,pop in [('MVP','EUR503'),('R13','FIN99')]:
 with np.load(W/f'stage6_hfrs_sparse_score_{pop}_genotypes.npz',allow_pickle=True) as f:
  G=np.asarray(f['genotypes'],dtype=float);ids=f['variants'].astype(str)
 assert list(ids)==list(d.variant) and G.shape==(503 if pop=='EUR503' else 99,15) and np.isfinite(G).all();a=[];zs=[]
 for weight in weights:
  g=comp[(comp.source==src)&(comp.weight==weight)].set_index('variant').loc[d.variant];freq=g.frequency_counted_allele.to_numpy();a.append(g.weight_counted_allele.to_numpy()*np.sqrt(2*freq*(1-freq)));zs.append(g.z_counted_allele.to_numpy())
 a=np.asarray(a);zs=np.asarray(zs);R=np.where(mask,np.corrcoef(G,rowvar=False),0.0);point=np.sum(a*zs,axis=1)/np.sqrt(np.einsum('ij,jk,ik->i',a,R,a));target=observed[(observed.source==src)&(observed.reference==pop)&(observed.omitted=='none')].set_index('weight').loc[weights,'weighted_z'].to_numpy();assert np.allclose(point,target,atol=1e-10)
 valid=[];dropped=0
 for k in range(B):
  boot=G[rng.integers(0,len(G),len(G))];v=boot.std(axis=0)
  if np.any(v==0):dropped+=1;continue
  r=np.where(mask,np.corrcoef(boot,rowvar=False),0.0);den=np.einsum('ij,jk,ik->i',a,r,a)
  if np.any(den<=0):dropped+=1;continue
  valid.append(np.sum(a*zs,axis=1)/np.sqrt(den))
 valid=np.asarray(valid)
 for j,weight in enumerate(weights):
  q=np.quantile(valid[:,j],[.025,.5,.975]);summary.append({'source':src,'reference':pop,'weight':weight,'observed_z':float(point[j]),'reference_bootstrap_2_5':float(q[0]),'reference_bootstrap_median':float(q[1]),'reference_bootstrap_97_5':float(q[2]),'valid_draws':len(valid),'discarded_draws':dropped})
  for k,value in enumerate(valid[:,j]):records.append({'source':src,'reference':pop,'weight':weight,'valid_draw_index':k,'weighted_z':value})
pd.DataFrame(summary).to_csv(O/'stage6_insomnia_sparse_score_ld_bootstrap_summary.csv',index=False);pd.DataFrame(records).to_csv(W/'stage6_insomnia_sparse_score_ld_bootstrap_draws.csv.gz',index=False,compression='gzip');(O/'stage6_insomnia_sparse_score_ld_bootstrap_receipt.json').write_text(json.dumps({'seed':61028,'draws_requested_per_reference':B,'sampling':'reference individuals resampled with replacement; target GWAS data and weights fixed','limitations':'Reference-only Monte Carlo sensitivity, not a confidence interval for target association, effect, replication or performance. Does not incorporate target GWAS uncertainty, missingness, ascertainment, selection, relatedness or reference mismatch.'},indent=2)+'\n');print(pd.DataFrame(summary).to_string(index=False))
