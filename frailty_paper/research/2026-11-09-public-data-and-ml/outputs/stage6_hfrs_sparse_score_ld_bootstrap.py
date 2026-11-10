#!/usr/bin/env python3
"""Resample FIN reference individuals only; study summary statistics and weights remain fixed."""
import numpy as np,pandas as pd,json
from pathlib import Path
O=Path('outputs');W=Path('work');v=np.load(W/'stage6_hfrs_sparse_score_FIN99_genotypes.npz',allow_pickle=True);G=v['genotypes'];variants=v['variants'];g=pd.read_csv(O/'candidate_finngen_hfrs_full_gwas_lookup.csv').set_index('variant').loc[variants];s=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv').set_index('variant').loc[variants];f=np.where(g.hfrs_source_effect_allele==g.candidate_counted_allele,g.hfrs_source_effect_allele_frequency,1-g.hfrs_source_effect_allele_frequency);scale=np.sqrt(2*f*(1-f));z=(g.hfrs_beta_on_counted_allele/g.hfrs_se).to_numpy();chr=s.chr.to_numpy();mask=chr[:,None]==chr[None,:];rng=np.random.default_rng(20261028);all_results=[];invalid=0
for j in range(5000):
 x=G[rng.integers(0,len(G),size=len(G))]; sd=x.std(axis=0)
 if (sd==0).any():invalid+=1;continue
 r=np.where(mask,np.corrcoef(x,rowvar=False),0)
 for wname in ['FI_beta_weight','FFS_beta_weight_inferred','insomnia_beta_weight']:
  a=s[wname].to_numpy()*scale;value=float(a@z/np.sqrt(a@r@a));all_results.append((j,wname,value))
d=pd.DataFrame(all_results,columns=['replicate','weight','weighted_z']);d.to_csv(O/'stage6_hfrs_sparse_score_FIN99_ld_bootstrap.csv.gz',index=False);summary=d.groupby('weight').weighted_z.agg(median='median',q025=lambda a:a.quantile(.025),q975=lambda a:a.quantile(.975));summary.to_csv(O/'stage6_hfrs_sparse_score_FIN99_ld_bootstrap_summary.csv');se_scale=g.hfrs_se.to_numpy()*scale;receipt={'seed':20261028,'attempted':5000,'discarded_zero_variance_reference_draws':invalid,'bootstrap_unit':'reference individuals','quantile_interpretation':'Monte Carlo sensitivity to reference-sampling variation only; not confidence intervals for effect or prediction','min_SE_times_HWE_SD':float(se_scale.min()),'median_SE_times_HWE_SD':float(np.median(se_scale)),'max_SE_times_HWE_SD':float(se_scale.max())};(O/'stage6_hfrs_sparse_score_FIN99_ld_bootstrap_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(summary);print(receipt)
