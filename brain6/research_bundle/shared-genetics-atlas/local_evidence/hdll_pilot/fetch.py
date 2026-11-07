import hashlib,io,json,zipfile
from pathlib import Path
import requests

p=Path(__file__).resolve().parent;work=Path('/Volumes/Brain6Work/hdll-chr3-44-20261028')
assert not work.exists();work.mkdir();(work/'bim').mkdir();(work/'ld').mkdir()
base='https://zenodo.org/api/records/14825987/files/'
bim=requests.get(base+'bimfile.zip/content',timeout=90);bim.raise_for_status();data=bim.content
assert len(data)==9673272 and hashlib.md5(data).hexdigest()=='608f037e4fc4d69f0ee9231c7fc005bb'
bimname='bimfile/ukb_chr3.44_n336000.imputed_clean.bim'
with zipfile.ZipFile(io.BytesIO(data)) as z:
 item=z.getinfo(bimname);bb=z.read(item);assert len(bb)==item.file_size
 (work/'bim'/Path(bimname).name).write_bytes(bb)
positions=[int(x.split()[3]) for x in bb.decode().splitlines()];assert len(positions)==681 and min(positions)==47598663 and max(positions)==50380661
url=base+'LDfile.zip/content';size=1830679554
class Remote(io.RawIOBase):
 def __init__(self):self.pos=0
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else size+n;return self.pos
 def read(self,n=-1):
  if n<0:n=size-self.pos
  n=min(n,size-self.pos)
  if n<=0:return b''
  r=requests.get(url,headers={'Range':f'bytes={self.pos}-{self.pos+n-1}'},params={'start':self.pos},timeout=90);r.raise_for_status()
  assert r.status_code==206 and r.headers.get('Content-Range')==f'bytes {self.pos}-{self.pos+n-1}/{size}' and len(r.content)==n
  self.pos+=n;return r.content
with zipfile.ZipFile(Remote()) as z:
 names=[n for n in z.namelist() if n.startswith('LD.path/ukb_chr3.44_') and n.endswith('_LDSVD.rda')];assert len(names)==1
 item=z.getinfo(names[0]);ld=z.read(item);assert len(ld)==item.file_size
 (work/'ld'/Path(names[0]).name).write_bytes(ld)
commit='cadb27c9aaa0e64a06a480cf2b0390a38f4f1d85'
codeurl=f'https://raw.githubusercontent.com/zhenin/HDL/{commit}/HDL/R/HDL.L.R';res=requests.get(codeurl,timeout=30);res.raise_for_status();code=res.content
assert b'HDL.L <- function' in code
(work/'HDL.L.R').write_bytes(code)
receipt={'zenodo_record':'https://doi.org/10.5281/zenodo.14825987','license':'CC BY 4.0','bim_archive_bytes':len(data),
 'bim_archive_md5':hashlib.md5(data).hexdigest(),'bim_archive_sha256':hashlib.sha256(data).hexdigest(),
 'bim_member':bimname,'bim_member_sha256':hashlib.sha256(bb).hexdigest(),'bim_count':len(positions),'min_pos':min(positions),'max_pos':max(positions),
 'ld_archive_bytes':size,'ld_archive_md5_publisher':'1d2aa5d3afee8b3bef7836f943476696','ld_archive_fully_verified':False,
 'ld_member':names[0],'ld_member_crc32':f'{item.CRC:08x}','ld_member_sha256':hashlib.sha256(ld).hexdigest(),'ld_member_bytes':len(ld),
 'source_commit':commit,'source_url':codeurl,'source_sha256':hashlib.sha256(code).hexdigest(),'work_directory':str(work)}
(p/'fetch_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
