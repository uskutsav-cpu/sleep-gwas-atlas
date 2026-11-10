import concurrent.futures,hashlib,os,requests
from pathlib import Path
u='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90854001-GCST90855000/GCST90854255/GCST90854255.tsv.gz';size=141860844;chunk=(size+7)//8;root=Path('work/hansen2026/parts');root.mkdir(exist_ok=True)
def get(i):
 lo=i*chunk;hi=min(size-1,(i+1)*chunk-1);p=root/str(i)
 with requests.get(u,headers={'Range':f'bytes={lo}-{hi}'},stream=True,timeout=(20,180)) as r:
  r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {lo}-{hi}/{size}'
  with p.open('wb') as f:
   for b in r.iter_content(1024*256):f.write(b)
 assert p.stat().st_size==hi-lo+1;print('completed',i,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as e:list(e.map(get,range(8)))
p=Path('work/hansen2026/parallel_verified.tsv.gz');h=hashlib.md5()
with p.open('wb') as w:
 for i in range(8):
  with (root/str(i)).open('rb') as f:
   for b in iter(lambda:f.read(1024*1024),b''):h.update(b);w.write(b)
assert h.hexdigest()=='f634c1ceadf8bd32e98b311aab8a973f';print('verified',h.hexdigest(),flush=True)
