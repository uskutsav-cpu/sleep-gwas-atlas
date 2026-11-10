'''Search unharmonized text by observed lexical chromosome:position order; fail closed if not found.'''
import urllib.request,concurrent.futures,json,hashlib,csv
from pathlib import Path
u='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90435001-GCST90436000/GCST90435387/GCST90435387.tsv';total=347498277

def get(a,b):
 with urllib.request.urlopen(urllib.request.Request(u,headers={'Range':f'bytes={a}-{b}'}),timeout=30) as r:
  if r.status!=206 or r.headers.get('Content-Range')!=f'bytes {a}-{b}/{total}':raise RuntimeError('bad range')
  x=r.read()
  if len(x)!=b-a+1:raise RuntimeError('incomplete range')
  return x

def rowkey(line):
 a=line.decode().split('\t')
 if len(a)<8 or not a[1].isdigit():raise RuntimeError('invalid source row')
 return a[0]+':'+a[1]

def find(variant,pos,pair):
 key='2:'+str(pos);lo=0;hi=total-1;probes=[]
 while hi-lo>60000:
  mid=(lo+hi)//2;buf=get(mid,min(mid+1024,total-1));lines=buf.split(b'\n')
  if len(lines)<3:raise RuntimeError('short record')
  got=rowkey(lines[1]);probes.append({'offset':mid,'observed_key':got})
  if got<key:lo=mid
  else:hi=mid
 a=max(0,lo-2000);b=min(total-1,hi+2000);chunk=get(a,b);lines=chunk.split(b'\n')[1:-1];hits=[]
 for ln in lines:
  f=ln.decode().split('\t')
  if len(f)>=8 and f[0]=='2' and f[1]==str(pos) and set(f[2:4])==pair:hits.append(f[:8])
 if len(hits)!=1:raise RuntimeError(f'Expected one exact row for {variant}, found {len(hits)}; sortedness assumption not established')
 return {'variant':variant,'query_key':key,'line':hits[0],'byte_range':[a,b],'chunk_sha256':hashlib.sha256(chunk).hexdigest(),'probes':probes}
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:
 out=list(e.map(lambda t:find(*t),[('rs113851554',66523432,{'G','T'}),('rs78206187',66795102,{'G','A'})]))
Path('outputs/stage6_rls2024_unharmonized_bounded_lookup.json').write_text(json.dumps({'url':u,'file_size_bytes':total,'full_file_md5_not_verified':True,'ordering_note':'Binary lookup based on observed lexical ordering; exact rows found are evidence, inability to find a row would not establish absence','results':out},indent=2)+'\n')
for x in out:print(x['variant'],x['line'])
