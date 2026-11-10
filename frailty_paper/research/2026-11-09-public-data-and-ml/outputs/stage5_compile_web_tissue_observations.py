#!/usr/bin/env python3
"""Normalize manually observed, rounded Atlas rows without inventing track IDs."""
import csv
from pathlib import Path
root=Path(__file__).resolve().parent

def rows(name):
 with (root/name).open(newline='') as f:return list(csv.DictReader(f))

queue={r['variant_key'] for r in rows('stage5_score_queue.csv')}
out=[]
for name in ['stage5_initial_atlas_tissue_observations.csv','stage5_rs60271_atlas_tissue_observations.csv','stage5_rs323509_atlas_initial_observations.csv']:
 for i,r in enumerate(rows(name),1):
  variant=r.get('variant_id') or r.get('variant_GRCh38');scorer=r.get('scorer') or r.get('output_type')
  if scorer=='AVI':continue
  assert variant in queue
  out.append({'variant_GRCh38_ref_alt':variant,'source_file':name,'source_data_row':i,'scorer':scorer,'assay':r.get('assay',''),'biosample':r.get('biosample',''),'reported_gene':r.get('reported_gene') or r.get('gene_id',''),'track_strand':r.get('track_strand') or r.get('strand',''),'displayed_raw_score':r.get('displayed_raw_score',''),'displayed_other_score':r.get('displayed_website_quantile') or r.get('displayed_other_score',''),'displayed_other_score_label':r.get('displayed_other_score_label') or 'Website quantile','source_url':r.get('source_url') or 'https://deepmind.google.com/science/alphagenome/atlas','source_limitations':r.get('caveat') or r.get('precision_and_limits',''),'comparability':'No full precision; event and/or tissue matching not established; not a validation result'})
assert len(out)==12
with (root/'stage5_rounded_tissue_observations_normalized.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=out[0].keys());w.writeheader();w.writerows(out)
print(len(out),'tissue observations;',len({r['variant_GRCh38_ref_alt'] for r in out}),'variants')
