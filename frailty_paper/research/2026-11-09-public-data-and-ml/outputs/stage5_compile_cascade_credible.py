#!/usr/bin/env python3
"""Extract portal-returned credible-set molecular evidence, retaining missingness."""
import csv,gzip,json,pathlib,hashlib
root=pathlib.Path('work/stage5_cascade_raw')
with open('outputs/stage5_score_queue.csv') as f:q=list(csv.DictReader(f))
fields=['variant_GRCh38','input_classes','input_labels','cell_type','qtl_type','feature_id','gene_symbol','pip','nlog10p','beta_as_reported','se','credible_set','peak_overlaps_variant','mechanism','portal_partial','portal_partial_reason','caveat'];out=[];available=0
for x in q:
 key=x['variant_key'];ident=f"{x['chr']}_{x['position_GRCh38']}_{x['reference_allele']}_{x['alternate_allele']}";f=root/(hashlib.sha256(key.encode()).hexdigest()[:20]+'.json.gz')
 if not f.exists():continue
 with gzip.open(f,'rt') as h:d=json.load(h)
 if d.get('empty'):continue
 assert d['variant']['variant_id']==ident;available+=1
 for group in ('per_ct','per_ct_l2'):
  for _,ct in (d.get(group) or {}).items():
   for typ,field in [('eQTL','eqtls'),('caQTL','caqtls')]:
    for item in ct.get(field,[]):
     out.append({'variant_GRCh38':key,'input_classes':x['input_classes'],'input_labels':x['input_labels'],'cell_type':ct.get('cell_type'),'qtl_type':typ,'feature_id':item.get('ensembl') or item.get('peak_id'),'gene_symbol':item.get('gene',''),'pip':item.get('pip'),'nlog10p':item.get('nlog10p'),'beta_as_reported':item.get('beta'),'se':item.get('se'),'credible_set':item.get('cs'),'peak_overlaps_variant':item.get('overlap',''),'mechanism':ct.get('mechanism'),'portal_partial':d.get('partial'),'portal_partial_reason':d.get('partial_reason'),'caveat':'Only portal-returned credible hits; allele orientation and source selection need verification; no causal or independent disease claim'})
with open('outputs/stage5_cascade_credible_qtl_by_cell.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
print('available variants',available,'credible records',len(out))
