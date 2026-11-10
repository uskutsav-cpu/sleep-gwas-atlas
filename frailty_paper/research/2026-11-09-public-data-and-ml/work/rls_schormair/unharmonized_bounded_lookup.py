'''Bounded lookups in original Schormair files; finding rows does not prove the file globally sorted.'''
import urllib.request,concurrent.futures,json,hashlib
from pathlib import Path
studies=[('GCST90399568',506157142),('GCST90399571',902988634)]
targets=[('rs113851554',66750564,{'G','T'}),('rs78206187',67022234,{'G','A'})]
def run(job):
 acc,size,rs,pos,pair=job;url=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90399001-GCST90400000/{acc}/{acc}.tsv'
 def get(a,b):
  with urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}'}),timeout=30) as r:
   if r.status!=206 or r.headers.get('Content-Range')!=f'bytes {a}-{b}/{size}':raise RuntimeError('bad range')
   x=r.read();assert len(x)==b-a+1;return x
 def key(line):
  f=line.decode().split('\t');return (int(f[0]),int(f[1]))
 lo=0;hi=size-1;trace=[]
 while hi-lo>60000:
  mid=(lo+hi)//2;ls=get(mid,min(mid+2048,size-1)).split(b'\n');assert len(ls)>=3;k=key(ls[1]);trace.append({'offset':mid,'key':k})
  if k<(2,pos):lo=mid
  else:hi=mid
 a=max(0,lo-2000);b=min(size-1,hi+2000);chunk=get(a,b);hits=[]
 for ln in chunk.split(b'\n')[1:-1]:
  f=ln.decode().split('\t')
  if len(f)>=9 and f[0]=='2' and f[1]==str(pos) and set(f[2:4])==pair and f[8]==rs:hits.append(f[:9])
 if len(hits)!=1:raise RuntimeError(f'{acc} {rs} expected 1 got {len(hits)}; no assumption of absence')
 return {'study':acc,'url':url,'variant':rs,'input_build':'GRCh37','row':hits[0],'byte_range':[a,b],'chunk_sha256':hashlib.sha256(chunk).hexdigest(),'probes':trace,'whole_file_checksum_verified':False}
jobs=[(acc,size,rs,pos,pair) for acc,size in studies for rs,pos,pair in targets]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:out=list(ex.map(run,jobs))
Path('outputs/stage6_rls_schormair_original_allele_lookup.json').write_text(json.dumps(out,indent=2)+'\n')
for v in out:print(v['study'],v['variant'],v['row'])
