#!/usr/bin/env python3
import csv,gzip,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parents[2]
prior=json.loads((root/'replication/mvp_six_leads/receipt.json').read_text());src=Path(prior['mvp_path'])
h=hashlib.sha256()
with src.open('rb') as f:
 while b:=f.read(1024*1024):h.update(b)
assert h.hexdigest()==prior['mvp_sha256']
ids={'rs6765037','rs7134291'};rows=[];n=0
with gzip.open(src,'rt') as f:
 reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames
 for n,x in enumerate(reader,1):
  if x['rsid'] in ids:rows.append(x)
assert n==prior['mvp_rows']
with (p/'mvp_source_rows.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
(p/'receipt.json').write_text(json.dumps({'input_path':str(src),'sha256':h.hexdigest(),'rows':n,'target_ids':sorted(ids),'counts':{x:sum(r['rsid']==x for r in rows) for x in sorted(ids)},'limitation':'rsID lookup is not full build, allele or cohort-overlap validation'},indent=2)+'\n')
print('rows',n,'target matches',len(rows))
for x in rows:print(x)
