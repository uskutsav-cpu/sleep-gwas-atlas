import io,subprocess
class Remote(io.RawIOBase):
 def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.offset=-1;self.cache=b'';self.transferred=0
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):
  self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n
  if self.pos<0:raise ValueError('negative offset')
  return self.pos
 def read(self,n=-1):
  n=min(self.size-self.pos,self.size if n<0 else n)
  if n<=0:return b''
  if n>16*1024*1024:raise ValueError('request too large')
  if not (self.offset<=self.pos and self.pos+n<=self.offset+len(self.cache)):
   start=self.pos;end=min(self.size,start+max(n,4*1024*1024))-1
   b=subprocess.run(['curl','-fsSL','--max-time','60','--retry','2','--range',f'{start}-{end}',self.url],capture_output=True,check=True).stdout
   if len(b)!=end-start+1:raise IOError('bad range')
   self.offset=start;self.cache=b;self.transferred+=len(b)
  i=self.pos-self.offset;out=self.cache[i:i+n];self.pos+=n;return out
