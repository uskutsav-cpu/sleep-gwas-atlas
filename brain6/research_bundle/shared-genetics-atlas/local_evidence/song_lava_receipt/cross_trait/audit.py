#!/usr/bin/env python3
"""Read saved locus 464 receipts and chromosome-3 input descriptors."""
import csv,hashlib,json
from pathlib import Path

base=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1/pairs')
out=Path(__file__).resolve().parent
items=[]
for file in sorted(base.glob('*/loci/locus_0464.json')):
 b=file.read_bytes();d=json.loads(b);trait=file.parents[1].name
 assert d['trait']==trait and d['LOC']=='464' and d['CHR']==3
 info=file.parents[1]/'input_info_chr03.tsv'
 with info.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 frailty=[r for r in rows if r['phenotype']=='frailty'];assert len(frailty)==1
 log=file.parents[1]/'logs/locus_0464.log';lb=log.read_bytes()
 items.append({'trait':trait,'receipt_sha256':hashlib.sha256(b).hexdigest(),
  'input_info_sha256':hashlib.sha256(info.read_bytes()).hexdigest(),
  'frailty_input_path':frailty[0]['filename'],'log_sha256':hashlib.sha256(lb).hexdigest(),
  'fi_negative_variance_warning':"Negative variance estimate for phenotype(s) 'frailty'" in lb.decode(errors='replace'),
  'receipt':d})
assert len(items)==12 and len({x['frailty_input_path'] for x in items})==1
assert sum(x['fi_negative_variance_warning'] for x in items)==11
frailty_file=Path(items[0]['frailty_input_path']);h=hashlib.sha256()
with frailty_file.open('rb') as f:
 for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
result={'frailty_input_sha256':h.hexdigest(),'n_pairs':len(items),'records':items}
(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps([{'trait':x['trait'],'n_snps':x['receipt']['n_snps'],'K':x['receipt']['K'],
 'fi_status':x['receipt']['fi_status'],'bivar_status':x['receipt']['bivar_status']} for x in items],indent=2))
