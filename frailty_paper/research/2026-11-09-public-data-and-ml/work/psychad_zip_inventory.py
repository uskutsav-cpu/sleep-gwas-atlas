import io,requests,zipfile
class HTTPFile(io.RawIOBase):
 def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.s=requests.Session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):self.pos=offset if whence==0 else self.pos+offset if whence==1 else self.size+offset;return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  if n<=0:return b''
  a=self.pos;b=min(a+n-1,self.size-1);r=self.s.get(self.url,headers={'Range':f'bytes={a}-{b}','Accept-Encoding':'identity'},timeout=45);r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range','').startswith(f'bytes {a}-{b}/');self.pos+=len(r.content);return r.content
j=requests.get('https://zenodo.org/api/records/17435917',timeout=30).json()
for f in j['files']:
 z=zipfile.ZipFile(HTTPFile(f['links']['self'],f['size']));names=z.infolist();print('\n',f['key'],'members',len(names));
 for x in names[:35]:print(x.filename,'uncompressed',x.file_size,'compressed',x.compress_size)
