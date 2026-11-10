#!/usr/bin/env python3
"""Extract selected public ZIP members by range and verify their ZIP CRC."""
import io,requests,zipfile,hashlib,json,sqlite3
from pathlib import Path
root=Path(__file__).resolve().parents[1];work=root/'work';out=root/'outputs'
class Remote(io.RawIOBase):
 def __init__(self,url,size):self.url=url;self.size=size;self.pos=0;self.cache=b'';self.cache_start=0;self.s=requests.Session();self.transferred=0
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else self.size+n;return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  if n<=0 or self.pos>=self.size:return b''
  if not(self.cache_start<=self.pos and self.pos+n<=self.cache_start+len(self.cache)):
   a=self.pos;b=min(a+max(n,4*1024*1024)-1,self.size-1);r=self.s.get(self.url,headers={'Range':f'bytes={a}-{b}','Accept-Encoding':'identity'},timeout=60);r.raise_for_status();assert r.status_code==206 and r.headers.get('Content-Range','').startswith(f'bytes {a}-{b}/');self.cache=r.content;self.cache_start=a;self.transferred+=len(r.content)
  start=self.pos-self.cache_start;data=self.cache[start:start+n];self.pos+=len(data);return data
record=requests.get('https://zenodo.org/api/records/17435917',timeout=30).json();f=next(x for x in record['files'] if x['key']=='snTIMs.zip');remote=Remote(f['links']['self'],f['size']);z=zipfile.ZipFile(remote);receipt=[]
for c in ['Endo','Mural']:
 member=f'snTIMs/PsychAD_snTWAS_EUR_Class_{c}.db'; info=z.getinfo(member);path=work/Path(member).name;h=hashlib.sha256()
 with z.open(info) as src,open(path,'wb') as dst:
  while True:
   chunk=src.read(1024*1024)
   if not chunk:break
   h.update(chunk);dst.write(chunk)
 # Reading to EOF verified the ZIP member's CRC.
 db=sqlite3.connect(f'file:{path}?mode=ro',uri=True);print(c,'tables',db.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall(),flush=True);db.close()
 receipt.append({'member':member,'crc32_hex':f'{info.CRC:08x}','member_crc_verified':True,'bytes':info.file_size,'sha256':h.hexdigest()})
(out/'stage6_psychad_small_models_receipt.json').write_text(json.dumps({'source':f['links']['self'],'whole_zip_checksum_verified':False,'bytes_downloaded':remote.transferred,'members':receipt},indent=2)+'\n');print('downloaded',remote.transferred,flush=True)
