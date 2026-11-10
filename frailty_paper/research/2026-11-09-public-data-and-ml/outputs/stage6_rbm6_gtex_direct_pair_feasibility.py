#!/usr/bin/env python3
"""Count exact RBM6 cortex–disease matches before reference LD matching."""
from pathlib import Path
import json
import pandas as pd
root=Path(__file__).resolve().parents[1];o=root/'outputs'; w=root/'work'
q=pd.read_csv(o/'stage6_rbm6_gtex_cortex_region.tsv',sep='\t');q=q[(q.type=='SNP') & (q.ref.str.len()==1) & (q.alt.str.len()==1)].copy();q['pair']=q.apply(lambda r:''.join(sorted([r.ref,r.alt])),axis=1);q=q[~q.pair.isin(['AT','CG']) & q.se.gt(0)].copy()
for key,g in q.groupby('variant'):
 for field in ['position','ref','alt','beta','se','ac','an']:
  assert g[field].nunique(dropna=False)==1,(key,field)
q=q.drop_duplicates('variant'); assert not q.duplicated(['position','pair']).any()
cs=pd.read_csv(o/'stage6_rbm6_gtex_cortex_published_cs_coverage.tsv',sep='\t').drop_duplicates(['cs_id','variant']); result={}
for trait in ['Mak_FI','Ye_FFS']:
 d=pd.read_csv(w/f'stage3_inka1_{trait}_aligned.csv');d['pair']=d.apply(lambda r:''.join(sorted([r.ref,r.alt])),axis=1); assert not d.duplicated(['position_GRCh38','pair']).any()
 x=q.merge(d,left_on=['position','pair'],right_on=['position_GRCh38','pair'],validate='one_to_one');ids=set(x['variant_x']);match=cs[cs.variant.isin(ids)];top=cs.loc[cs.pip.idxmax()]
 result[trait]={'qtl_valid_nonpalindromic_snps':len(q),'direct_exact_matches_before_reference_ld':len(x),'published_cs_unique':len(cs),'published_cs_matched':len(match),'published_cs_pip_sum_matched':float(match.pip.sum()),'published_cs_pip_sum_total':float(cs.pip.sum()),'top_variant_matched':bool(top.variant in ids),'top_rsid':top.rsid,'scope':'Position and allele-pair match before LD and frequency QC; NOT colocalization or independent validation.'}
 x[['variant_x','position','ref_x','alt_x','beta','se','gwas_beta_alt','gwas_se','gwas_p']].to_csv(o/f'stage6_rbm6_gtex_direct_{trait}_candidates.csv',index=False)
(o/'stage6_rbm6_gtex_direct_pair_feasibility.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
