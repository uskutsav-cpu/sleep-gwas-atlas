#!/usr/bin/env python3
import io,json,subprocess,zipfile
from pathlib import Path
p=Path(__file__).resolve().parent;meta=json.loads((p/'archive_inventory.json').read_text());size=meta['archive_bytes'];url=meta['url']
class Remote(io.RawIOBase):
 def __init__(self):self.pos=0;self.offset=-1;self.cache=b''
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):
  self.pos=n if whence==0 else self.pos+n if whence==1 else size+n
  if self.pos<0:raise ValueError('negative offset')
  return self.pos
 def read(self,n=-1):
  n=min(size-self.pos,size if n<0 else n)
  if n<=0:return b''
  if n>16*1024*1024:raise ValueError('request too large')
  if self.offset<=self.pos and self.pos+n<=self.offset+len(self.cache):
   i=self.pos-self.offset;out=self.cache[i:i+n]
  else:
   start=self.pos;end=min(size,start+max(n,2*1024*1024))-1
   out=subprocess.run(['curl','-fsSL','--max-time','40','--range',f'{start}-{end}',url],capture_output=True,check=True).stdout
   if len(out)!=end-start+1:raise IOError('bad range')
   self.offset=start;self.cache=out;out=out[:n]
  self.pos+=n;return out
with zipfile.ZipFile(Remote()) as z:
 content=z.read('ldm13M/ldm.info');(p/'ldm.info').write_bytes(content)
 lines=content.decode().splitlines();print('block info header',lines[0]);print('block341',*[x for x in lines if x.split() and x.split()[0]=='341'])
 with z.open('ldm13M/snp.info') as f:
  for i in range(3):print('snp info',f.readline().decode().rstrip())
