#!/usr/bin/env python3
"""Preserve selected CASCADE top-per-cell and credible QTL rows; not exhaustive."""
import argparse,csv,gzip,json,pathlib,hashlib,math
p=argparse.ArgumentParser();p.add_argument('--rawdir',default='work/stage5_cascade_raw');a=p.parse_args();root=pathlib.Path(a.rawdir)
with open('outputs/stage5_score_queue.csv') as f:q=list(csv.DictReader(f))
fields=['variant_GRCh38','input_classes','input_labels','cell_type','qtl_type','record_scope','feature_id','gene_symbol','p_value','q_value','beta_as_reported','se','allele_frequency_as_reported','pip','credible_set','partial','partial_reason','portal_url','caveat'];rows=[];coverage=[]
for x in q:
 key=x['variant_key'];variant=f"{x['chr']}_{x['position_GRCh38']}_{x['reference_allele']}_{x['alternate_allele']}";f=root/(hashlib.sha256(key.encode()).hexdigest()[:20]+'.json.gz')
 if not f.exists():coverage.append({'variant_GRCh38':key,'status':'not available','partial':'','n_cell_contexts':'','note':'No retrieved API record; not evidence of no QTL'});continue
 try:
  with gzip.open(f,'rt') as z:d=json.load(z)
  if d.get('empty') is True:
   assert f"{x['chr']}:{x['position_GRCh38']}" in d.get('reason','')
   coverage.append({'variant_GRCh38':key,'status':'portal reports no position results','partial':'','n_cell_contexts':0,'note':d.get('reason','')});continue
  assert d['variant']['variant_id']==variant
 except Exception:coverage.append({'variant_GRCh38':key,'status':'invalid','partial':'','n_cell_contexts':'','note':'Response could not be verified'});continue
 marginal=d.get('marginal_per_ct') or {};url='https://cascade.finngen.fi/variant/'+variant
 for ct,v in marginal.items():
  for typ,field in [('eQTL','eqtl_top'),('caQTL','caqtl_top')]:
   item=v.get(field)
   if not item:continue
   rows.append({'variant_GRCh38':key,'input_classes':x['input_classes'],'input_labels':x['input_labels'],'cell_type':ct,'qtl_type':typ,'record_scope':'browser top hit per variant/cell/type; not all features','feature_id':item.get('phenotype_id',''),'gene_symbol':item.get('gene_symbol',''),'p_value':item.get('p',''),'q_value':item.get('qval',''),'beta_as_reported':item.get('beta',''),'se':item.get('se',''),'allele_frequency_as_reported':item.get('af',''),'pip':'','credible_set':'','partial':d.get('partial'),'partial_reason':d.get('partial_reason'),'portal_url':url,'caveat':'Finnish PBMC multiome; beta allele convention must be checked for this portal; top-hit selection, sample overlap and colocalization limitations apply'})
 coverage.append({'variant_GRCh38':key,'status':'retrieved exact variant ID','partial':d.get('partial'),'n_cell_contexts':len(marginal),'note':d.get('partial_reason') or ''})
with open('outputs/stage5_cascade_top_qtl_by_cell.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
with open('outputs/stage5_cascade_query_coverage.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(coverage[0]));w.writeheader();w.writerows(coverage)
print('queue',len(q),'retrieved',sum(x['status']=='retrieved exact variant ID' for x in coverage),'top rows',len(rows))
