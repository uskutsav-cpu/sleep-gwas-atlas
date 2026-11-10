import io,subprocess,zipfile,json,collections
url='https://zenodo.org/api/records/6883292/files/GSTRIDE_database.zip/content'; size=1371421673
class Range(io.RawIOBase):
 def __init__(self):self.pos=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,offset,whence=0):
  self.pos=(offset if whence==0 else self.pos+offset if whence==1 else size+offset);return self.pos
 def read(self,n=-1):
  if n<0:n=size-self.pos
  n=min(n,size-self.pos)
  if not n:return b''
  r=subprocess.run(['curl','-fsSL','--max-time','20','-r',str(self.pos)+'-'+str(self.pos+n-1),url],capture_output=True,check=True).stdout
  if len(r)!=n:raise ValueError('server did not provide exact range')
  self.pos+=n;return r
with zipfile.ZipFile(Range()) as z:
 entries=[x for x in z.infolist() if not x.is_dir()]
 print('files',len(entries),'types',dict(collections.Counter(x.filename.rsplit('.',1)[-1].lower() for x in entries)))
 for x in entries:
  if any(y in x.filename.lower() for y in ['register','description','readme','metadata']):print('candidate',x.filename,'uncompressed_bytes',x.file_size)
