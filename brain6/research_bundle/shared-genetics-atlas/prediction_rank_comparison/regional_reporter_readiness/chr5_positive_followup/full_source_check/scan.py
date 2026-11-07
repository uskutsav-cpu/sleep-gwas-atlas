import csv, gzip, hashlib, json
from pathlib import Path
out=Path(__file__).resolve().parent
files={
 'frailty_raw':Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/data/raw/frailty.txt.gz'),
 'adhd_ldsc_munged':Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/adhd.sumstats.gz')}
results={}
for name,path in files.items():
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
 matches=[]; count=0; malformed=0
 with gzip.open(path,'rt') as f:
  reader=csv.DictReader(f,delimiter='\t')
  for row in reader:
   count+=1
   if name=='frailty_raw':
    match=(row['variant_id']=='rs12659965' or (row['chromosome'] in ('5','chr5') and row['base_pair_location']=='103904399'))
   else: match=row['SNP']=='rs12659965'
   if match: matches.append(row)
   if None in row.values(): malformed+=1
 results[name]={'path':str(path),'sha256':h.hexdigest(),'rows':count,'malformed_rows':malformed,'matches':matches}
(out/'receipt.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2))
