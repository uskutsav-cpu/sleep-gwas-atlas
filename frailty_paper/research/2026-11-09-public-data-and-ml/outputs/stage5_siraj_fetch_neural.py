import requests,io,zipfile
class RangeFile(io.RawIOBase):
 def __init__(self,u):
  self.u=u;self.s=requests.Session();r=self.s.get(u,headers={'Range':'bytes=0-0'},timeout=30);r.raise_for_status();self.size=int(r.headers['Content-Range'].split('/')[-1]);self.pos=0
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,x,whence=0):
  self.pos=({0:0,1:self.pos,2:self.size}[whence]+x);return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  if n==0 or self.pos>=self.size:return b''
  end=min(self.pos+n,self.size)-1;r=self.s.get(self.u,headers={'Range':f'bytes={self.pos}-{end}'},timeout=120);r.raise_for_status()
  if r.status_code!=206 or len(r.content)!=end-self.pos+1:raise IOError('range mismatch')
  self.pos=end+1;return r.content

import sys,zipfile,pathlib,hashlib,json,concurrent.futures

u='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content';root=pathlib.Path(sys.argv[1]);root.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(RangeFile(u)) as z:members=[x.filename for x in z.infolist() if '/processed_results/DE/' in x.filename and '_SKNSH_' in x.filename and x.filename.endswith('.gz')]
def fetch(name):
 target=root/name.split('/')[-1]
 if not target.exists():
  with zipfile.ZipFile(RangeFile(u)) as z:data=z.read(name)
  target.write_bytes(data)
 data=target.read_bytes();print('downloaded',target.name,len(data),flush=True);return {'member':name,'file':str(target),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:r=list(pool.map(fetch,members))
pathlib.Path('outputs/stage5_siraj2026_neural_source_receipt.json').write_text(json.dumps({'url':u,'selection':'all SKNSH processed allele-differential result files','files':r},indent=2)+'\n');print('completed',len(r))
