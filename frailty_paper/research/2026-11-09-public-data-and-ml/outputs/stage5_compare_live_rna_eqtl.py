#!/usr/bin/env python3
"""Conservative allele/gene/context matches; direction screen, not validation."""
import argparse,csv,gzip,hashlib,json,math,re
from pathlib import Path
from collections import defaultdict,Counter
p=argparse.ArgumentParser();p.add_argument('--store',required=True);p.add_argument('--output',default='outputs/stage5_live_rna_eqtl_direction_screen.csv');a=p.parse_args()
root=Path('outputs');store=Path(a.store)
def norm(s):return re.sub('[^a-z0-9]+','',s.lower())
with (root/'stage5_variant_gene_tissue_evidence_atlas.csv').open(newline='') as f:qtl=list(csv.DictReader(f))
lookup=defaultdict(list)
for q in qtl:
 if not q['gene_id'] or not q['qtl_alt_effect']:continue
 try:beta=float(q['qtl_alt_effect'])
 except ValueError:continue
 if not math.isfinite(beta) or beta==0:continue
 if q['tissue_or_cell']=='CD4_Naive' and ('eqtl' in q['molecular_measure'].lower() or 'gene expression' in q['molecular_measure'].lower()):context=('cell','naive thymus-derived CD4-positive, alpha-beta T cell')
 elif q['molecular_measure'].lower() in ('eqtl','gene expression') and q['tissue_or_cell']:context=('gtex',norm(q['tissue_or_cell']))
 else:continue
 lookup[(q['variant_GRCh38_ref_alt'],q['gene_id'].split('.')[0],*context)].append(q)
fields=['variant_GRCh38_ref_alt','gene_id','gene_name','tissue_match_method','biosample_name','gtex_tissue','assay','track_name','ontology_curie','model_raw_score','model_quantile_score','qtl_source_collection','qtl_source','qtl_dataset','qtl_tissue','qtl_measure','qtl_event','qtl_alt_effect','qtl_p_value','possible_exact_duplicate_across_collections','sign_screen','caveat']
count=Counter()
with open(a.output,'w',newline='') as out:
 w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
 for receiptfile in store.glob('*.json'):
  if receiptfile.name=='run_config.json' or receiptfile.name.startswith('._'):continue
  r=json.loads(receiptfile.read_text());key=r.get('variant_key');score=store/r.get('file','')
  if r.get('status')!='complete' or not score.is_file() or hashlib.sha256(score.read_bytes()).hexdigest()!=r.get('sha256'):continue
  count['valid_variants']+=1
  with gzip.open(score,'rt',newline='') as f:
   for s in csv.DictReader(f):
    if s['output_type']!='RNA_SEQ' or not s['gene_id']:continue
    contexts=[]
    if s.get('gtex_tissue'):contexts.append(('gtex',norm(s['gtex_tissue'])))
    if s.get('biosample_name')=='naive thymus-derived CD4-positive, alpha-beta T cell':contexts.append(('cell',s['biosample_name']))
    for mode,context in contexts:
     for q in lookup.get((key,s['gene_id'].split('.')[0],mode,context),[]):
      try:pred=float(s['raw_score']);beta=float(q['qtl_alt_effect'])
      except (ValueError,TypeError):continue
      if not math.isfinite(pred):continue
      direction='zero_or_unresolved' if pred==0 else ('same_sign' if pred*beta>0 else 'opposite_sign')
      w.writerow(dict(variant_GRCh38_ref_alt=key,gene_id=s['gene_id'],gene_name=s['gene_name'],tissue_match_method='tissue labels normalized to model GTEx tissue' if mode=='gtex' else 'naive CD4 label mapping; samples differ',biosample_name=s.get('biosample_name',''),gtex_tissue=s.get('gtex_tissue',''),assay=s.get('Assay title',''),track_name=s.get('track_name',''),ontology_curie=s.get('ontology_curie',''),model_raw_score=s['raw_score'],model_quantile_score=s.get('quantile_score',''),qtl_source_collection=q['qtl_source_collection'],qtl_source=q['qtl_source'],qtl_dataset=q['qtl_dataset'],qtl_tissue=q['tissue_or_cell'],qtl_measure=q['molecular_measure'],qtl_event=q['molecular_event'],qtl_alt_effect=q['qtl_alt_effect'],qtl_p_value=q['qtl_p_value'],possible_exact_duplicate_across_collections=q['possible_exact_duplicate_across_collections'],sign_screen=direction,caveat='Exploratory direction only; tracks and studies differ, small scores may be unstable; GTEx may overlap model training; no causal or independent validation claim'))
      count['rows']+=1;count[direction]+=1
print(dict(count))
