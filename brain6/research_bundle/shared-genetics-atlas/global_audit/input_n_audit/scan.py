import csv,gzip,hashlib,json,collections,math
from pathlib import Path
out=Path(__file__).resolve().parent; atlas=out.parent.parent
m=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged');f=Path('/Volumes/Extreme SSD/sleep-frailty-full-tp-v3-20261009/munged_full_native')
files={'insomnia':m/'insomnia.sumstats.gz','adhd':m/'adhd.sumstats.gz','frailty_index':f/'frailty.sumstats.gz','fried_frailty_score':f/'fried_frailty_score.sumstats.gz'}
old=[json.loads((atlas/'global_audit'/x/'receipt.json').read_text())['source_sha256'] for x in ('insomnia_adhd_rerun','fi_ffs_rerun')]
records={}
for name,path in files.items():
 h=hashlib.sha256()
 with path.open('rb') as fd:
  for b in iter(lambda:fd.read(1048576),b''):h.update(b)
 digest=h.hexdigest(); previous=[d[str(path)] for d in old if str(path) in d];assert previous and all(x==digest for x in previous)
 c=collections.Counter();invalid=0;missing=0;n=0
 with gzip.open(path,'rt') as fd:
  r=csv.DictReader(fd,delimiter='\t');assert 'N' in r.fieldnames
  for row in r:
   n+=1
   if row['N'] is None or row['N'].strip().upper() in ('','NA','NAN'):
    missing+=1;continue
   try:
    value=float(row['N'])
    if not math.isfinite(value) or value<=0:raise ValueError()
    c[value]+=1
   except (ValueError,TypeError):invalid+=1
 records[name]={'path':str(path),'sha256':digest,'matches_rerun_receipt':True,'rows':n,'missing_n':missing,'invalid_n':invalid,'distinct_n':len(c),'min_n':min(c),'max_n':max(c),'most_common_n':c.most_common(5)}
(out/'receipt.json').write_text(json.dumps(records,indent=2)+'\n');print(json.dumps(records,indent=2))
