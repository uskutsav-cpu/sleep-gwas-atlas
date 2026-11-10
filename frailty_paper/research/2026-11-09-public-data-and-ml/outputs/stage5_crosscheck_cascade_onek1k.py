#!/usr/bin/env python3
"""Cross-check direction at a selected, correlated SMAD5 molecular credible set."""
import csv,collections
with open('outputs/stage5_variant_gene_tissue_evidence_atlas.csv') as f:q=[x for x in csv.DictReader(f) if x['qtl_source']=='outputs/stage4_cd4_smad5_credible_set.csv' and x['gene_id']=='ENSG00000113658']
with open('outputs/stage5_cascade_top_qtl_by_cell.csv') as f:c={(x['variant_GRCh38'],x['feature_id'],x['cell_type']):x for x in csv.DictReader(f) if x['qtl_type']=='eQTL' and x['cell_type']=='l2.CD4_Naive'}
rows=[]
for x in q:
 y=c.get((x['variant_GRCh38_ref_alt'],x['gene_id'],'l2.CD4_Naive'))
 if not y:continue
 try:
  a=float(x['qtl_alt_effect']);b=float(y['beta_as_reported']);sgn='same' if a*b>0 else ('opposite' if a*b<0 else 'zero')
 except (ValueError,TypeError):sgn='unavailable'
 rows.append({'variant_GRCh38':x['variant_GRCh38_ref_alt'],'gene_id':x['gene_id'],'oneK_cell_type':x['tissue_or_cell'],'oneK_alt_beta':x['qtl_alt_effect'],'oneK_p':x['qtl_p_value'],'oneK_pip':x['molecular_pip'],'cascade_cell_type':y['cell_type'],'cascade_beta_as_reported':y['beta_as_reported'],'cascade_p':y['p_value'],'cascade_af_as_reported':y['allele_frequency_as_reported'],'sign_screen':sgn,'caveat':'Variants are correlated within one molecular credible set; samples, methods, scales and selection differ; comparison does not prove portal effect orientation or disease colocalization.'})
if not rows:raise RuntimeError('No matched records; check inputs')
with open('outputs/stage5_cascade_onek1k_smad5_crosscheck.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('oneK set',len(q),'matched records',len(rows),'signs',dict(collections.Counter(x['sign_screen'] for x in rows)))
