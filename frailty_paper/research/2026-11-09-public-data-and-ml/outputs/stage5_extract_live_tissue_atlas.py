#!/usr/bin/env python3
"""Stream verified AlphaGenome outputs into a labeled, non-exhaustive tissue view."""
import argparse,csv,gzip,hashlib,json,re
from pathlib import Path
from collections import Counter
p=argparse.ArgumentParser();p.add_argument('--store',required=True);p.add_argument('--output',default='outputs/stage5_live_tissue_evidence_atlas.csv.gz');a=p.parse_args();store=Path(a.store)
with open('outputs/stage5_variant_multitrait_gwas_matrix.csv',newline='') as f:gw={r['variant_GRCh38_ref_alt']:r for r in csv.DictReader(f)}
qtl_gene_count=Counter()
with open('outputs/stage5_variant_gene_tissue_evidence_atlas.csv',newline='') as f:
 for q in csv.DictReader(f):
  if q['gene_id'] and q['qtl_source']:
   qtl_gene_count[(q['variant_GRCh38_ref_alt'],q['gene_id'].split('.')[0])]+=1
traits=['Atkins_FI','Jansen_insomnia','Mak_FI','Mak_FP','Ye_FFS','FinnGen_R13_insomnia']
base=['variant_GRCh38_ref_alt','candidate_gene_context','ccre_overlap_ids','same_variant_gene_qtl_source_rows','tissue_category_screen','gene_id','gene_name','gene_type','gene_strand','output_type','variant_scorer','track_name','track_strand','Assay title','ontology_curie','biosample_name','biosample_type','biosample_life_stage','gtex_tissue','data_source','endedness','genetically_modified','histone_mark','transcription_factor','junction_Start','junction_End','raw_score','quantile_score','scored_interval']
extra=[t+'_'+x for t in traits for x in ['match_status','beta_alt','se','p']]
patterns={'brain_neural':r'brain|cerebr|cerebell|cortex|cortical|hippocamp|hypothalam|thalam|amygdal|caudate|putamen|substantia nigra|neuron|neuronal|astrocyte|oligodendro|microglia|glial|glia','pituitary':r'pituitary','muscle':r'muscle|myoblast|myocyte|myotube','immune':r'immune|lymphocyt|leukocyt|monocyt|macrophage|dendritic cell|granulocyt|neutrophil|eosinophil|basophil|natural killer|\bt.cell\b|\bb.cell\b|\bcd4\b|\bcd8\b|thymus|spleen|blood|hematopoietic'}
rx={k:re.compile(v,re.I) for k,v in patterns.items()};count=Counter()
with gzip.open(a.output,'wt',newline='') as out:
 w=csv.DictWriter(out,fieldnames=base+extra,extrasaction='ignore');w.writeheader()
 for receiptfile in store.glob('*.json'):
  if receiptfile.name=='run_config.json' or receiptfile.name.startswith('._'):continue
  r=json.loads(receiptfile.read_text());key=r.get('variant_key');file=store/r.get('file','')
  if r.get('status')!='complete' or key not in gw or not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest()!=r.get('sha256'):continue
  count['verified_variants']+=1
  with gzip.open(file,'rt',newline='') as f:
   for s in csv.DictReader(f):
    text=' '.join(s.get(k,'') for k in ['biosample_name','gtex_tissue','track_name']);labels=[label for label,pattern in rx.items() if pattern.search(text)]
    if s['output_type']=='SPLICE_SITES' and not labels:labels=['splice_site_no_tissue_assignment']
    if not labels:continue
    row={k:s.get(k,'') for k in base};row['variant_GRCh38_ref_alt']=key;row['candidate_gene_context']=gw[key]['candidate_gene_context'];row['ccre_overlap_ids']=gw[key]['ccre_ids'];row['same_variant_gene_qtl_source_rows']=qtl_gene_count[(key,s['gene_id'].split('.')[0])] if s['gene_id'] else '';row['tissue_category_screen']='|'.join(labels);row.update({col:gw[key][col] for col in extra});w.writerow(row);count['rows']+=1
    for label in labels:count[label]+=1
print(dict(count))
