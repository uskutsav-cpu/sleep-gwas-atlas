#!/usr/bin/env python3
import csv,gzip,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
src=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/track_b_finngen_r13_F5_ADHD.sumstats.gz')
h=hashlib.sha256()
with src.open('rb') as f:
 while b:=f.read(1024*1024):h.update(b)
assert h.hexdigest()=='24d4d72ebdf1f83e7af4dd5e46aa5208266e121f14bd3d5a26217a87f861f59b'
selected=[];n=0
with gzip.open(src,'rt') as f:
 for n,r in enumerate(csv.DictReader(f,delimiter='\t'),1):
  if n%1000==500:selected.append({'row':n,**r})
assert n==1156825 and len(selected)==1157 and len({x['SNP'] for x in selected})==len(selected)
with (p/'sample.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=selected[0],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(selected)
(p/'selection_receipt.json').write_text(json.dumps({'source':str(src),'source_sha256':h.hexdigest(),'rows':n,'selected':len(selected),'rule':'data row modulo 1000 equals 500'},indent=2)+'\n')
print('Selected',len(selected),'from',n,'rows')
