#!/usr/bin/env python3
"""Query public CASCADE variant browser; preserve top-per-cell and raw evidence."""
import argparse,csv,hashlib,json,pathlib,time,gzip,requests,concurrent.futures
p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=0);p.add_argument('--rawdir',default='work/stage5_cascade_raw');a=p.parse_args();root=pathlib.Path(a.rawdir);root.mkdir(parents=True,exist_ok=True)
with open('outputs/stage5_score_queue.csv') as f:q=list(csv.DictReader(f))
if a.limit:q=q[:a.limit]
def task(x):
 key=x['variant_key'];variant=f"{x['chr']}_{x['position_GRCh38']}_{x['reference_allele']}_{x['alternate_allele']}";f=root/(hashlib.sha256(key.encode()).hexdigest()[:20]+'.json.gz')
 if f.exists():
  try:
   with gzip.open(f,'rt') as z:d=json.load(z)
   if d.get('variant',{}).get('variant_id')==variant or (d.get('empty') is True and f"{x['chr']}:{x['position_GRCh38']}" in d.get('reason','')):return (key,'cached',len(d.get('marginal_per_ct') or {}))
  except Exception:pass
 url='https://cascade.finngen.fi/api/cascade/variant/'+variant
 for i in range(4):
  try:
   res=requests.get(url,timeout=45);res.raise_for_status();d=res.json()
   if d.get('variant',{}).get('variant_id')!=variant and not (d.get('empty') is True and f"{x['chr']}:{x['position_GRCh38']}" in d.get('reason','')):raise ValueError('response variant mismatch')
   tmp=f.with_suffix('.tmp');
   with gzip.open(tmp,'wt') as z:json.dump(d,z,separators=(',',':'))
   tmp.replace(f);time.sleep(0.2);return(key,'downloaded',len(d.get('marginal_per_ct') or {}))
  except Exception as e:
   if i==3:return(key,'failed:'+type(e).__name__,0)
   time.sleep(2**i)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 for key,status,n in pool.map(task,q):print(status,key,'contexts',n,flush=True)
