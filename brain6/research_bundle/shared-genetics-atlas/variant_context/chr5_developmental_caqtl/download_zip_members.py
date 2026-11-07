import io,subprocess,zipfile,json,hashlib,struct,zlib,sys
from pathlib import Path
p=Path(__file__).resolve().parent;url='https://zenodo.org/api/records/7820257/files/data.zip/content';size=6322256088;tail_start=size-131072;tail=Path('work/zenodo782_zip_tail.bin').read_bytes();assert len(tail)==131072
class Remote(io.RawIOBase):
 def __init__(self):self.pos=0
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.pos
 def seek(self,n,whence=0):self.pos=n if whence==0 else self.pos+n if whence==1 else size+n;return self.pos
 def read(self,n=-1):
  if n<0:n=size-self.pos
  n=min(n,size-self.pos);assert n<=4_000_000
  if n==0:return b''
  if self.pos>=tail_start:b=tail[self.pos-tail_start:self.pos-tail_start+n]
  else:b=subprocess.run(['curl','-fsSL','--max-time','30','-r',f'{self.pos}-{self.pos+n-1}',url],capture_output=True,check=True).stdout
  assert len(b)==n;self.pos+=n;return b
out=Path('/Volumes/Brain6Work/xqtl-chr5/developmental_caqtl');out.mkdir(exist_ok=True);receipt=[]
with zipfile.ZipFile(Remote()) as z:
 for cell in ('neuron','progenitor'):
  name=f'data/Model1A/caQTL/{cell}/chr5.txt.gz';item=z.getinfo(name);assert item.compress_type==8 and item.file_size<100_000_000
  h=Remote();h.seek(item.header_offset);data=h.read(30);sig,ver,flags,method,tm,dt,crc,cs,us,nlen,elen=struct.unpack('<IHHHHHIIIHH',data);assert sig==0x04034b50 and method==8 and not flags&1
  filename=h.read(nlen).decode();assert filename==name;start=item.header_offset+30+nlen+elen;end=start+item.compress_size-1
  dest=out/(cell+'_chr5.txt.gz');assert not dest.exists();hdr=out/(cell+'_range_headers.txt')
  proc=subprocess.Popen(['curl','-fsSL','--max-time','300','-r',f'{start}-{end}','-D',str(hdr),url],stdout=subprocess.PIPE)
  dec=zlib.decompressobj(-15); checksum=0;count=0;hsh=hashlib.sha256()
  with dest.open('wb') as o:
   while True:
    chunk=proc.stdout.read(1024*1024)
    if not chunk:break
    plain=dec.decompress(chunk);o.write(plain);hsh.update(plain);checksum=zlib.crc32(plain,checksum);count+=len(plain)
   plain=dec.flush();o.write(plain);hsh.update(plain);checksum=zlib.crc32(plain,checksum);count+=len(plain)
  assert proc.wait()==0 and dec.eof and count==item.file_size and (checksum & 0xffffffff)==item.CRC
  hh=hdr.read_text().lower();assert f'content-range: bytes {start}-{end}/{size}' in hh
  print(cell,'inner_bytes',count,'CRC matches',flush=True)
  receipt.append({'cell':cell,'member':name,'zip_local_header_offset':item.header_offset,'zip_data_start':start,'zip_data_end':end,'compressed_bytes':item.compress_size,'inner_bytes':count,'zip_crc32':format(item.CRC,'08x'),'crc_verified':True,'inner_gzip_sha256':hsh.hexdigest(),'file':str(dest)})
(p/'retrieval_receipt.json').write_text(json.dumps({'source':'https://doi.org/10.5281/zenodo.7820257','archive_bytes':size,'archive_md5_publisher':'60afce365b0c76a795c8863a99496bb6','archive_fully_hashed':False,'members':receipt,'limitation':'CRC verifies member transport consistency, not full archive authenticity or phenotype interpretation.'},indent=2)+'\n')
