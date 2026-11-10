#!/usr/bin/env python3
"""Inventory existing GRCh37 reference flags; do not infer NCDS exclusions."""
import csv,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent.parent;out=root/'outputs'
spec=[('15','unvalidated_15_variant_score_prototype.csv','GRCh37_ensembl_reference_alleles'),
      ('22','unvalidated_22_shared_direction_score.csv','reference_alleles')]
result={}
for label,name,col in spec:
 p=out/name
 with p.open() as f:rows=list(csv.DictReader(f))
 assert len(rows)==int(label) and len({r['variant'] for r in rows})==len(rows)
 hits=[{'variant':r['variant'],'reference_alleles':r[col], 'counted':r['counted_allele'],'comparator':r['other_allele']} for r in rows if len(r[col].split('/'))>2]
 result[label]={'n':len(rows),'reference_multiple_allele_labels':hits,'count':len(hits),'input_sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
receipt={'results':result,'limitations':'Source-file GRCh37 reference annotation, not NCDS variant file, site frequencies or proof of removal; splits, genomic build and alternate representations require review.'}
(out/'stage6_ncds_multiallelic_reference_inventory.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:{'n':v['n'],'multiple_label_count':v['count'],'variants':[x['variant'] for x in v['reference_multiple_allele_labels']]} for k,v in result.items()},indent=2))
