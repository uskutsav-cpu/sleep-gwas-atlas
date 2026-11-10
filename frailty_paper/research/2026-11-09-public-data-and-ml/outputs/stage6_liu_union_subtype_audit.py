#!/usr/bin/env python3
"""Compare rs9536347-overlapping subtype intervals with the public union peaks."""
import csv,hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent.parent;out=root/'outputs';union=root/'work/stage6_liu_union_peaks.bed'
receipt=json.loads((out/'stage6_liu_brain_ccre_candidate_overlaps.json').read_text());meta=next(x for x in receipt['files'] if x['file']=='Union_PeakSet.bed');digest=hashlib.sha256(union.read_bytes()).hexdigest();assert digest==meta['sha256']
with (out/'stage6_liu_brain_ccre_candidate_overlaps.csv').open(newline='') as f:sub=[x for x in csv.DictReader(f) if x['variant']=='rs9536347' and x['source_file']!='Union_PeakSet.bed']
assert len(sub)==3
unions=[]
with union.open() as f:
 for line in f:
  c,a,b=line.rstrip('\n').split('\t')[:3]
  if c=='chr13':
   a,b=int(a),int(b)
   if a<53071000 and b>53068000:unions.append((a,b))
rows=[]
for s in sub:
 a,b=int(s['peak_start']),int(s['peak_end']);pos=int(s['position_grch38']);hits=[(u,v,min(b,v)-max(a,u)) for u,v in unions if max(a,u)<min(b,v)]
 for u,v,overlap in hits:rows.append({'variant':s['variant'],'position':pos,'subtype_file':s['source_file'],'subtype_start':a,'subtype_end':b,'union_start':u,'union_end':v,'interval_overlap_length_assuming_half_open':overlap,'variant_inside_union_strict':u<pos<v})
assert len(rows)==3 and all(not x['variant_inside_union_strict'] for x in rows)
with (out/'stage6_liu_union_subtype_audit.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
result={'source_union_sha256':digest,'subtype_overlap_table_sha256':hashlib.sha256((out/'stage6_liu_brain_ccre_candidate_overlaps.csv').read_bytes()).hexdigest(),'rows':rows,'interpretation':'Subtype intervals overlap union intervals even though the variant lies outside them. Consistent with iterative peak merging but does not prove actual file provenance or which peak won; source-file assembly and coordinate conventions remain conditional.'}
(out/'stage6_liu_union_subtype_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(rows,indent=2))
