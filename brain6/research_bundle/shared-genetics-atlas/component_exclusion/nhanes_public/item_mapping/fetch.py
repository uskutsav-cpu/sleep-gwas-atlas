#!/usr/bin/env python3
import concurrent.futures,hashlib,json,subprocess
from pathlib import Path
out=Path(__file__).resolve().parent
root=Path('/Volumes/Brain6NHANES/source');assert root.is_dir()
files=[(y,f'{n}_{s}') for y,s in [('2011','G'),('2013','H')] for n in
       ['MCQ','BPQ','DIQ','KIQ_U','HUQ','RXQ_RX','BMX','GHB','CBC','PAQ','WHQ','SLQ']]
def fetch(x):
 y,n=x;url=f'https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{y}/DataFiles/{n}.XPT';path=root/(n+'.XPT')
 head=subprocess.run(['curl','-fsSIL','--max-time','30',url],capture_output=True,text=True,check=True)
 sizes=[int(x.split(':',1)[1].strip()) for x in head.stdout.splitlines() if x.lower().startswith('content-length:')]
 assert sizes;expected=sizes[-1]
 if not path.exists():
  subprocess.run(['curl','-fsSL','--max-time','90',url,'-o',str(path)],check=True)
 assert path.is_file() and path.stat().st_size==expected
 h=hashlib.sha256()
 with path.open('rb') as f:
  assert f.read(80).startswith(b'HEADER RECORD*******LIBRARY HEADER RECORD');f.seek(0)
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return {'file':path.name,'url':url,'bytes':expected,'sha256':h.hexdigest()}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex: rows=list(ex.map(fetch,files))
assert len(rows)==24
(out/'download_receipt.json').write_text(json.dumps(rows,indent=2)+'\n')
print('Verified',len(rows),'files,',sum(r['bytes'] for r in rows),'bytes; existing files not overwritten.')
