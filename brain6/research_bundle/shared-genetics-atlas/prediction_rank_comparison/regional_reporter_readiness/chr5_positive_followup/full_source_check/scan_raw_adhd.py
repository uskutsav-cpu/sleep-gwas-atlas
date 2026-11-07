import gzip,hashlib,json
from pathlib import Path
p=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/raw/adhd.txt.gz')
h=hashlib.sha256()
with p.open('rb') as f:
 for block in iter(lambda:f.read(1048576),b''):h.update(block)
with gzip.open(p,'rt') as f:
 header=next(f).split(); count=0; malformed=0; matches=[]
 for line in f:
  parts=line.split(); count+=1
  if len(parts)!=len(header):malformed+=1;continue
  row=dict(zip(header,parts))
  if row['SNP']=='rs12659965' or (row['CHR'] in ('5','chr5') and row['BP']=='103904399'):matches.append(row)
r={'path':str(p),'sha256':h.hexdigest(),'rows':count,'malformed_rows':malformed,'matches':matches}
Path(__file__).with_name('raw_adhd_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
