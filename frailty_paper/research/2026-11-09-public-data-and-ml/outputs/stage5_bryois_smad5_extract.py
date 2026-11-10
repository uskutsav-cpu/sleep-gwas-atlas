#!/usr/bin/env python3
"""Download checksum-verified public Bryois summary files and retain SMAD5 records."""
import csv,gzip,hashlib,json,os,tempfile
from pathlib import Path
import requests
out=Path(__file__).resolve().parent
files=[('Excitatory.neurons.5.gz','86dc91ce73d8ea13b6e263a0b667305d'),('snp_pos.txt.gz','309df692b9a7e23553445395ae1c7eab')]; receipt=[]
def download(name,expected):
 h=hashlib.md5();url=f'https://zenodo.org/api/records/7276971/files/{name}/content';r=requests.get(url,stream=True,timeout=90);r.raise_for_status();n=0
 fd,path=tempfile.mkstemp(prefix='bryois_',dir=out.parent/'work')
 try:
  with os.fdopen(fd,'wb') as f:
   for b in r.iter_content(1024*1024):
    if b:f.write(b);h.update(b);n+=len(b)
  if h.hexdigest()!=expected:raise ValueError('checksum mismatch')
 except Exception:
  os.unlink(path);raise
 receipt.append({'url':url,'bytes':n,'md5':h.hexdigest(),'checksum_verified':True});return path
p=download(*files[0]);rows=[]
try:
 with gzip.open(p,'rt') as f:
  for line in f:
   v=line.strip().split()
   if len(v)==5 and v[0]=='SMAD5_ENSG00000113658':rows.append(v)
finally:os.unlink(p)
if not rows:raise ValueError('No SMAD5 results')
ids={r[1] for r in rows};print('SMAD5 rows',len(rows),'distinct SNP IDs',len(ids),flush=True)
p=download(*files[1]);positions=[]
try:
 with gzip.open(p,'rt') as inner:
  for line in inner:
   v=line.strip().split('\t')
   if len(v)>=6 and v[0] in ids:positions.append(v[:6])
finally:os.unlink(p)
with open(out/'stage5_bryois_exc_smad5_associations.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['gene_id','snp_id','distance_to_tss','nominal_p','beta']);w.writerows(rows)
with open(out/'stage5_bryois_exc_smad5_variant_positions.csv','w',newline='') as f:
 w=csv.writer(f);w.writerow(['snp_id','snp_id_hg38','snp_id_hg19','effect_allele','other_allele','MAF']);w.writerows(positions)
(out/'stage5_bryois_exc_smad5_extract_receipt.json').write_text(json.dumps({'files':receipt,'smad5_rows':len(rows),'distinct_snp_ids':len(ids),'matching_position_rows':len(positions)},indent=2)+'\n')
print('matching position rows',len(positions),flush=True)
