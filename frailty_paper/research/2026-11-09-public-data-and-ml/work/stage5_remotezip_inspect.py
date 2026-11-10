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
u='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content'
with zipfile.ZipFile(RangeFile(u)) as z:
 for x in z.infolist():print(x.filename,x.file_size,x.compress_size)
