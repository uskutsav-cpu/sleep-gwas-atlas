#!/usr/bin/env python3
"""Read-only coordinate overlap audit of the archived insomnia LAVA receipts."""
import hashlib, json
from pathlib import Path

base=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1')
out=Path(__file__).resolve().parent
files=sorted((base/'pairs/insomnia/loci').glob('locus_*.json'))
assert len(files)==2495
matches=[]
for f in files:
 raw=f.read_bytes();d=json.loads(raw)
 assert d['trait']=='insomnia'
 if d['CHR']==3 and d['START']<=49316972 and d['STOP']>=47727212:
  log=base/'pairs/insomnia/logs'/f.with_suffix('.log').name
  matches.append({'receipt_path':str(f),'receipt_sha256':hashlib.sha256(raw).hexdigest(),
   'log_path':str(log),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'receipt':d})
assert len(matches)==1
record={'source_interval':{'chromosome':3,'start':47727212,'stop':49316972,'build':'not verified'},
 'receipts_read':len(files),'matching_rule':'inclusive numerical overlap; coordinate-build equivalence unverified',
 'matches':matches}
(out/'receipt.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
