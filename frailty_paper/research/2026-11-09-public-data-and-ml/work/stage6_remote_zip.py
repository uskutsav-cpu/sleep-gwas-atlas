import io,requests
class HTTPFile(io.RawIOBase):
 def __init__(self,url,size,chunk=4194304):self.url=url;self.size=size;self.chunk=chunk;self.pos=0;self.cache_start=-1;self.cache=b'';self.session=requests.Session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset
  if self.pos<0:raise ValueError('negative offset')
  return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(n,self.size-self.pos);out=bytearray()
  while n>0:
   if not(self.cache_start<=self.pos<self.cache_start+len(self.cache)):
    start=(self.pos//self.chunk)*self.chunk;end=min(start+self.chunk,self.size)-1;r=self.session.get(self.url,headers={'Range':f'bytes={start}-{end}'},timeout=90);r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range','').startswith(f'bytes {start}-');self.cache_start=start;self.cache=r.content
   i=self.pos-self.cache_start;k=min(n,len(self.cache)-i);assert k>0;out+=self.cache[i:i+k];self.pos+=k;n-=k
  return bytes(out)
