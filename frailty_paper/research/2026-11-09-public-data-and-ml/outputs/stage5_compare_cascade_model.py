#!/usr/bin/env python3
"""Exact allele and gene, explicit cell-label links; no unsupported directional claims."""
import argparse,csv,gzip,hashlib,pathlib,collections
p=argparse.ArgumentParser();p.add_argument('--store',required=True);a=p.parse_args();store=pathlib.Path(a.store)
with open('outputs/stage5_cascade_top_qtl_by_cell.csv') as f:q=[r for r in csv.DictReader(f) if r['qtl_type']=='eQTL' and r['feature_id'].startswith('ENSG')]
labels={'l2.CD4_Naive':'naive thymus-derived CD4-positive, alpha-beta T cell','l1.CD4_T':'CD4-positive, alpha-beta T cell','l2.CD8_Naive':'naive thymus-derived CD8-positive, alpha-beta T cell','l2.B_naive':'naive B cell','l1.B':'B cell','l1.NK':'natural killer cell','l1.Mono':'monocyte'}
index=collections.defaultdict(list)
for r in q:
 if r['cell_type'] in labels:index[(r['variant_GRCh38'],r['feature_id'].split('.')[0],labels[r['cell_type']])].append(r)
bykey=collections.defaultdict(list)
for k in index:bykey[k[0]].append(k)
fields=['variant_GRCh38','gene_id','gene_name','cascade_cell_type','model_biosample','cascade_p','cascade_q','cascade_beta_as_reported','cascade_se','cascade_af_as_reported','model_assay','model_track','model_raw_alt_minus_ref','model_quantile','cascade_partial','comparison_status','caveat'];n=0
with open('outputs/stage5_cascade_model_eqtl_comparison.csv','w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
 for variant in bykey:
  f=store/(hashlib.sha256(variant.encode()).hexdigest()[:20]+'.csv.gz')
  if not f.is_file():continue
  with gzip.open(f,'rt',newline='') as source:
   for s in csv.DictReader(source):
    if s['output_type']!='RNA_SEQ' or not s['gene_id']:continue
    k=(variant,s['gene_id'].split('.')[0],s['biosample_name'])
    for qtl in index.get(k,[]):
     w.writerow({'variant_GRCh38':variant,'gene_id':s['gene_id'],'gene_name':s['gene_name'],'cascade_cell_type':qtl['cell_type'],'model_biosample':s['biosample_name'],'cascade_p':qtl['p_value'],'cascade_q':qtl['q_value'],'cascade_beta_as_reported':qtl['beta_as_reported'],'cascade_se':qtl['se'],'cascade_af_as_reported':qtl['allele_frequency_as_reported'],'model_assay':s['Assay title'],'model_track':s['track_name'],'model_raw_alt_minus_ref':s['raw_score'],'model_quantile':s['quantile_score'],'cascade_partial':qtl['partial'],'comparison_status':'allele and gene match; cell labels mapped, beta orientation not yet confirmed','caveat':'Portal top feature only; samples and assays differ; beta allele convention, model training independence and cohort overlap unverified; not colocalization'});n+=1
print('comparison rows',n,'variants with candidate label matches',len(bykey))
