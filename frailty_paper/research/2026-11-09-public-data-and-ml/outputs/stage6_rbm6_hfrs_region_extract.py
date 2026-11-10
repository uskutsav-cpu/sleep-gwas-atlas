#!/usr/bin/env python3
import csv,hashlib,gzip,json
from pathlib import Path
root=Path(__file__).resolve().parents[1];p=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90707291.tsv');h=hashlib.md5();n=kept=0
with p.open('rb') as f,gzip.open(root/'work/stage6_rbm6_hfrs_2mb.tsv.gz','wt') as o:
 for l in f:
  h.update(l);n+=1
  if n==1:o.write(l.decode());continue
  a=l.split(b'\t',3)
  if a[0]==b'3' and 49008613<=int(a[1])<=51008613:o.write(l.decode());kept+=1
assert h.hexdigest()=='38d9434f50a4493def9d2c9d1cc1b4b3'
(root/'outputs/stage6_rbm6_hfrs_region_receipt.json').write_text(json.dumps({'source':str(p),'md5':h.hexdigest(),'source_rows':n-1,'selected':kept,'GRCh38_region':'chr3:49008613-51008613'},indent=2)+'\n');print('selected',kept,'total',n-1)
