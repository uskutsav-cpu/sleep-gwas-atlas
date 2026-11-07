#!/usr/bin/env python3
import hashlib,io,json,subprocess,zipfile
from pathlib import Path
p=Path(__file__).resolve().parent
url='https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/ukbEUR_13M_FullLDM.zip';size=206566549726
class Remote(io.RawIOBase):
 def __init__(self):self.offset=0;self.ranges=[]
 def seekable(self):return True
 def readable(self):return True
 def tell(self):return self.offset
 def seek(self,offset,whence=0):
  self.offset=offset if whence==0 else self.offset+offset if whence==1 else size+offset
  if self.offset<0:raise ValueError('negative offset')
  return self.offset
 def read(self,n=-1):
  n=min(size-self.offset,size if n<0 else n)
  if n<=0:return b''
  if n>32*1024*1024:raise ValueError('inspection read too large')
  first=self.offset;last=first+n-1
  data=subprocess.run(['curl','-fsSL','--max-time','30','--range',f'{first}-{last}',url],capture_output=True,check=True).stdout
  if len(data)!=n:raise IOError(f'wrong range length {len(data)} != {n}')
  self.offset+=n;self.ranges.append({'start':first,'end':last,'sha256':hashlib.sha256(data).hexdigest()});return data
r=Remote()
with zipfile.ZipFile(r) as z:
 all=z.infolist(); matches=[x for x in all if any(y in x.filename.lower() for y in ('block341.','block_341.','block341/','ldm.info','snp.info'))]
 result={'url':url,'archive_bytes':size,'member_count':len(all),'matched_members':[{'filename':x.filename,'uncompressed_bytes':x.file_size,'compressed_bytes':x.compress_size,'compression_method':x.compress_type,'crc32':f'{x.CRC:08x}','header_offset':x.header_offset} for x in matches],'ranges':r.ranges,'note':'Central-directory metadata only; neither selected members nor the full archive have been verified.'}
(p/'archive_inventory.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
