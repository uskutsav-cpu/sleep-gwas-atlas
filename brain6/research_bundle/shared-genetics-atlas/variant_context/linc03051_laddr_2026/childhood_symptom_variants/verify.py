#!/usr/bin/env python3
import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;d=json.loads((p/'receipt.json').read_text())
with (p/'results.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==6 and all(x['status']=='MATCHED' for x in rows)
for x in rows:
 calc=math.erfc(abs(float(x['symptom_beta_a1'])/float(x['symptom_se']))/math.sqrt(2))
 assert abs(calc-float(x['symptom_p']))<1e-12
 assert (float(x['symptom_beta_a1'])*float(x['brain_adhd_beta_a1'])>0)==(x['direction_agreement']=='YES')
assert hashlib.sha256((p/'results.tsv').read_bytes()).hexdigest()==d['table_sha256']
o={'rows':6,'minimum_symptom_p':min(float(x['symptom_p']) for x in rows),
 'maximum_symptom_p':max(float(x['symptom_p']) for x in rows),'all_nominal_above_0_05':all(float(x['symptom_p'])>.05 for x in rows),
 'receipt_sha256':hashlib.sha256((p/'receipt.json').read_bytes()).hexdigest(),
 'note':'Arithmetic checks, not independent source authentication or overlap verification.'}
(p/'verification.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))
