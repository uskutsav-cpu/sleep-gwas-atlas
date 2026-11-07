#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;d=json.loads((p/'receipt.json').read_text())
with (p/'ld.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==15 and len(d['variants'])==6
assert all(x['matching_records']==1 and x['eligible']==503 for x in d['variants'].values())
for x in rows:
 r=float(x['r']);r2=float(x['r2']);assert abs(r)<=1 and 0<=r2<=1 and math.isclose(r*r,r2,rel_tol=1e-12)
 assert int(x['n_individuals'])==503 and int(x['n_haplotypes'])==1006
assert hashlib.sha256((p/'ld.tsv').read_bytes()).hexdigest()==d['ld_table_sha256']
result={'n_variants':6,'pairs':len(rows),'n_samples':503,'table_sha256':d['ld_table_sha256'],
 'source_receipt_sha256':hashlib.sha256((p/'receipt.json').read_bytes()).hexdigest(),
 'note':'Arithmetic and saved receipt check only; remote VCF not independently hashed.'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
