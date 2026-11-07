#!/usr/bin/env python3
import csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent;receipt=json.loads((p/'receipt.json').read_text());selected=[]
for name,info in receipt['files'].items():
 path=p/(name.replace('.tsv.gz','')+'_selected.tsv');assert hashlib.sha256(path.read_bytes()).hexdigest()==info['selected_sha256']
 with path.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 assert len(rows)==info['rows_selected']
 if name.startswith('twas_hits'):
  for r in rows:
   if r['TRAIT']=='UKB_1200_Sleeplessness_or_insomnia':
    assert r['ID'].startswith('ENSG00000239268__PC') and r['BEST.GWAS.ID']=='chr3_117918850_T_C_b38'
    selected.append((name,r))
assert sum(x[0]=='twas_hits.gtex-ddp.tsv.gz' for x in selected)==3
assert sum(x[0]=='twas_hits.gtex-rddp.tsv.gz' for x in selected)==5
assert sum(x[0]=='twas_hits.gtex-kdp.tsv.gz' for x in selected)==0
out=p/'sleep_hits.tsv'
fields=['model_set','TISSUE','TRAIT','ID','BEST.GWAS.ID','EQTL.ID','TWAS.Z','TWAS.P','COLOC.PP3','COLOC.PP4']
with out.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
 for name,r in selected:w.writerow({'model_set':name.split('.')[1],**{k:r[k] for k in fields[1:]}})
result={'sleep_hits':len(selected),'ddp':3,'rddp':5,'kdp':0,'all_have_lead':'chr3_117918850_T_C_b38',
 'sleep_hits_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),
 'receipt_sha256':hashlib.sha256((p/'receipt.json').read_bytes()).hexdigest(),
 'note':'Checks selected files and recorded full-stream receipts, not external source identity or statistical calibration.'}
(p/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
