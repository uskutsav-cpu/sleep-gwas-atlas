#!/usr/bin/env python3
import hashlib,json,math
from pathlib import Path
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parent.parent; out=root/'outputs'
base=root/'work/niagads_smad5_susie';base.mkdir(exist_ok=True)
ref=pd.read_csv(out/'stage6_tenk10k_expanded_eur_ld_order.csv').reset_index()
data=np.load(out/'stage6_tenk10k_expanded_eur_ld.npz')
assert np.array_equal(ref.variant_GRCh38.to_numpy(),data['variant'])
d=pd.read_csv(out/'stage6_niagads_smad5_exc_Mak_FI_aligned.csv')
d=d[(d.qtl_alt_af-d.gwas_af_alt).abs()<=.1].copy()
d['key']=d.variant_GRCh38.str.removeprefix('chr')
r=ref.merge(d,left_on='variant_GRCh38',right_on='key',how='inner',validate='one_to_one')
r=r[(r.alt_frequency-r.qtl_alt_af).abs()<=.1].sort_values('index')
assert len(r)>100 and r['index'].is_unique
def bf(beta,se,sd):
 v=se*se;w=sd*sd; z=beta/se; return .5*(np.log1p(-w/(w+v))+w/(w+v)*z*z)
def frac(x,mask):
 x=np.asarray(x); m=float(x.max()); a=np.exp(x-m);return float(a[mask].sum()/a.sum())
mask=d.key.isin(r.key).to_numpy()
stats={'source_n':len(d),'matched_ld_n':len(r),'reference_n':503,
       'qtl_support_fraction':frac(bf(d.qtl_beta_alt,d.qtl_se,.15),mask),
       'fi_support_fraction':frac(bf(d.gwas_beta_alt,d.gwas_se,1.1004195728364333),mask),
       'qtl_lead_included':'5:136206558:C:T' in set(r.key),
       'candidate_included':'5:136202653:A:C' in set(r.key),
       'assumptions':'Conditional on already matched nonpalindromic SNPs, frequency filter and cis window; external EUR LD, not cohort LD.'}
r.to_csv(base/'aligned.csv',index=False)
idx=r['index'].to_numpy();matrix=np.asarray(data['r'][np.ix_(idx,idx)],dtype='<f8',order='F')
matrix.ravel(order='F').tofile(base/'ld_f64.bin')
stats['ld_sha256']=hashlib.sha256((base/'ld_f64.bin').read_bytes()).hexdigest()
(out/'stage6_niagads_smad5_susie_coverage.json').write_text(json.dumps(stats,indent=2)+'\n')
print(json.dumps(stats,indent=2))
