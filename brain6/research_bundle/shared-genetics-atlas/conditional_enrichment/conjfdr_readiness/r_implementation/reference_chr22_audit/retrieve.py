import hashlib,io,json,os,requests,zipfile
from pathlib import Path
url='https://zenodo.org/records/5750318/files/genref.zip?download=1';total=3012662640
out=Path('/Volumes/Extreme SSD/brain6-work/cfdr-reference-check');out.mkdir(exist_ok=False)
class Remote(io.RawIOBase):
 def __init__(self):self.p=0;self.s=requests.Session()
 def readable(self):return True
 def seekable(self):return True
 def tell(self):return self.p
 def seek(self,o,w=0):self.p={0:0,1:self.p,2:total}[w]+o;return self.p
 def read(self,n=-1):
  if n<0:n=total-self.p
  if not n:return b''
  end=min(total,self.p+n)-1;res=self.s.get(url,headers={'Range':f'bytes={self.p}-{end}'},timeout=90);res.raise_for_status();assert res.status_code==206 and len(res.content)==end-self.p+1;self.p=end+1;return res.content
z=zipfile.ZipFile(Remote());names=['all_chr_perVariant.rds','chr22_LDpairs.rds'];records=[]
for name in names:
 info=z.getinfo(name);dst=out/name;h=hashlib.sha256();count=0
 with z.open(info) as source,dst.open('xb') as target:
  while True:
   b=source.read(4*1024*1024)
   if not b:break
   target.write(b);h.update(b);count+=len(b)
 assert count==info.file_size
 records.append({'name':name,'bytes':count,'sha256':h.hexdigest(),'zip_crc32':format(info.CRC,'08x'),'crc_verified_on_full_read':True,'local_path':str(dst)})
 print(name,count,flush=True)
p=Path(__file__).resolve().parent;(p/'retrieval_receipt.json').write_text(json.dumps({'url':url,'archive_bytes':total,'members':records,'limitation':'HTTP ranges only; member CRC verified, full archive checksum not verified.'},indent=2)+'\n')
