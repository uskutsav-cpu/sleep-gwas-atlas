#!/usr/bin/env python3
from pathlib import Path
import json,pandas as pd,numpy as np
root=Path(__file__).resolve().parents[1];o=root/'outputs'
a=pd.read_csv(o/'stage5_singlebrain_ext_smad5_full_regional.tsv.gz',sep='\t',dtype={'ref':str,'alt':str,'Allele':str}); a['variant_key']=a.chr+':'+a.pos.astype(str)+':'+a.ref+'>'+a.alt
source=a.copy()
a=a[a.fixed_sd>0].copy()
assert a.variant_key.is_unique
p=pd.read_csv(o/'stage5_singlebrain_priority_variant_matches.csv');p=p[(p.gene=='SMAD5')&(p.cell_type=='Ext')];m=p.merge(a,on='variant_key',how='left',validate='one_to_one',indicator=True); assert (m._merge=='both').all()
valid=m.fixed_sd.gt(0)&np.isfinite(m.fixed_beta); assert valid.all();m['beta_source_alt']=np.where(m.Allele==m.alt,m.fixed_beta,np.where(m.Allele==m.ref,-m.fixed_beta,np.nan));m['portal_beta_source_diff']=m.portal_beta-m.fixed_beta
assert m.beta_source_alt.notna().all(); assert m.portal_beta_source_diff.abs().max()<1e-12
m['p_round_relative_difference']=(m.portal_random_p-m.Random_P).abs()/m.Random_P
assert m.p_round_relative_difference.max()<.006
matrix=pd.read_csv(o/'stage5_variant_multitrait_gwas_matrix.csv'); assert matrix.variant_GRCh38_ref_alt.is_unique
m=m.merge(matrix,on='variant_GRCh38_ref_alt',how='left') if 'variant_GRCh38_ref_alt' in m else m.merge(matrix,left_on='variant_key',right_on='variant_GRCh38_ref_alt',how='left',validate='one_to_one')
m.to_csv(o/'stage5_singlebrain_ext_smad5_source_validated_matches.csv',index=False)
receipt={'source_rows_in_gene_region':len(source),'source_rows_positive_fixed_se':len(a),'source_allele_alt':int((a.Allele==a.alt).sum()),'source_allele_ref':int((a.Allele==a.ref).sum()),'portal_queue_match_rows':len(m),'max_portal_beta_absolute_difference':float(m.portal_beta_source_diff.abs().max()),'max_random_p_round_relative_difference':float(m.p_round_relative_difference.max()),'original_marker_in_source':bool((a.variant_key=='chr5:136202653:A>C').any()),'alternative_in_source':bool((a.variant_key=='chr5:136207208:G>A').any())}
(o/'stage5_singlebrain_ext_smad5_source_audit.json').write_text(json.dumps(receipt,indent=2)+'\n'); print(json.dumps(receipt,indent=2))
