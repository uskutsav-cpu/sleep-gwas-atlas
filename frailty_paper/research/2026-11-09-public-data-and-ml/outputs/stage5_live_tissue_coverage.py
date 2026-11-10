#!/usr/bin/env python3
"""Summarize observed rows by screened category and output type for all queue entries."""
import argparse,csv,gzip,hashlib,json
from pathlib import Path
from collections import Counter
p=argparse.ArgumentParser();p.add_argument('--store',required=True);p.add_argument('--atlas',default='outputs/stage5_live_tissue_evidence_atlas.csv.gz');p.add_argument('--output',default='outputs/stage5_live_tissue_coverage.csv');a=p.parse_args();store=Path(a.store)
with open('outputs/stage5_score_queue.csv',newline='') as f:queue=list(csv.DictReader(f))
counts=Counter()
with gzip.open(a.atlas,'rt',newline='') as f:
 for r in csv.DictReader(f):
  for label in r['tissue_category_screen'].split('|'):counts[(r['variant_GRCh38_ref_alt'],label,r['output_type'])]+=1
labels=['brain_neural','pituitary','muscle','immune','splice_site_no_tissue_assignment'];types=['RNA_SEQ','SPLICE_SITES','SPLICE_SITE_USAGE','SPLICE_JUNCTIONS','ATAC','DNASE','CHIP_HISTONE','CHIP_TF','CAGE','PROCAP']
fields=['variant_GRCh38_ref_alt','input_classes','input_labels','result_status','category','output_type','matching_score_rows','interpretation']
with open(a.output,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
 for q in queue:
  key=q['variant_key'];name=hashlib.sha256(key.encode()).hexdigest()[:20];rp=store/(name+'.json');sp=store/(name+'.csv.gz');status='not scored'
  if rp.exists():
   try:
    receipt=json.loads(rp.read_text());valid=receipt.get('status')=='complete' and receipt.get('variant_key')==key and sp.is_file() and hashlib.sha256(sp.read_bytes()).hexdigest()==receipt.get('sha256');status='file verified' if valid else 'failed or invalid'
   except (OSError,ValueError):status='failed or invalid'
  for label in labels:
   for t in types:
    n=counts[(key,label,t)];w.writerow({'variant_GRCh38_ref_alt':key,'input_classes':q['input_classes'],'input_labels':q['input_labels'],'result_status':status,'category':label,'output_type':t,'matching_score_rows':n if status=='file verified' else '','interpretation':'Observed screen-matched rows; not an exhaustive ontology or unique biological samples' if n else ('No rows matched this label screen; not evidence of a zero effect or universal absence' if status=='file verified' else 'Not assessed')})
print('queued variants',len(queue),'coverage rows',len(queue)*len(labels)*len(types))
