#!/usr/bin/env python3
import sqlite3,csv,hashlib,gzip,json
from pathlib import Path
root=Path(__file__).resolve().parents[1];p=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/raw/GCST90707291.tsv');expected='38d9434f50a4493def9d2c9d1cc1b4b3';db=sqlite3.connect(root/'work/stage6_psychad_rbm6_model.db');ids={r[0].encode() for r in db.execute('select rsid from weights')};db.close();h=hashlib.md5();results=[];n=0
with p.open('rb') as f:
 for l in f:
  h.update(l);n+=1
  if n==1:header=l.decode().strip().split('\t');continue
  if l.rstrip(b'\r\n').rsplit(b'\t',1)[-1] in ids:results.append(l.decode().rstrip('\r\n').split('\t'))
if h.hexdigest()!=expected:raise ValueError('source checksum mismatch')
with gzip.open(root/'work/stage3_region_cache/FinnGen_HFRS_R12.regional.tsv.gz','wt') as out:
 w=csv.writer(out,delimiter='\t');w.writerow(header);w.writerows(results)
(root/'outputs/stage6_psychad_hfrs_extract_receipt.json').write_text(json.dumps({'source':str(p),'source_md5_verified':h.hexdigest(),'source_rows':n-1,'selected_rows':len(results)},indent=2)+'\n');print('selected',len(results),'records',n-1,flush=True)
