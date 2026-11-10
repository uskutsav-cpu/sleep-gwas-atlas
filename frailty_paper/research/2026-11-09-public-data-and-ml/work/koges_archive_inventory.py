import csv, json, struct, subprocess
from pathlib import Path
records=[('SAIGE_summaries.zip',26072103381),('POLMM_summaries.zip',3429270714)]
rows=[]
for name,size in records:
 url=f'https://zenodo.org/api/records/7042518/files/{name}/content'
 start=max(0,size-100000)
 cmd=['curl','-fsSL','-r',f'{start}-{size-1}','--max-filesize','200000','-w','\nSTATUS:%{http_code}',url]
 raw=subprocess.run(cmd,check=True,stdout=subprocess.PIPE).stdout
 body,status=raw.rsplit(b'\nSTATUS:',1)
 if status!=b'206' or len(body)!=size-start:raise RuntimeError('Range response incorrect; abort rather than download full archive')
 # Require end-of-central-directory plus consistent record counts; these archives use ZIP64.
 if b'PK\x05\x06' not in body:raise RuntimeError('Missing archive end record')
 names=[];p=0
 while True:
  i=body.find(b'PK\x01\x02',p)
  if i<0:break
  if i+46>len(body):raise RuntimeError('Truncated entry')
  n,e,c=struct.unpack_from('<HHH',body,i+28);end=i+46+n+e+c
  if end>len(body):raise RuntimeError('Truncated entry name')
  names.append(body[i+46:i+46+n].decode('utf-8'));p=end
 expected={'SAIGE_summaries.zip':43,'POLMM_summaries.zip':10}[name]
 if len(names)!=expected:raise RuntimeError(f'Expected {expected}, found {len(names)}')
 for entry in names:rows.append({'archive':name,'file':entry,'source':url})
target=Path('outputs/stage6_koges_archive_trait_inventory.csv')
with target.open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['archive','file','source']);w.writeheader();w.writerows(rows)
print(json.dumps({'files':len(rows),'archive_counts':{n:sum(x['archive']==n for x in rows) for n,_ in records},'keyword_hits':[x for x in rows if any(t in x['file'].lower() for t in ['sleep','insom','frail','grip'])]},indent=2))
