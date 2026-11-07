#!/usr/bin/env python3
"""Stream a public tar.bz2, recording a digest and selected gene rows."""
import hashlib,json,tarfile,urllib.request
from pathlib import Path

out=Path(__file__).resolve().parent
url='https://pantry.pejlab.org/data/twas_associations/PGC_ADHD_EUR_2017.tar.bz2'
class Reader:
 def __init__(self,s):self.s=s;self.h=hashlib.sha256();self.n=0
 def read(self,n=-1):
  b=self.s.read(n);self.h.update(b);self.n+=len(b);return b
with urllib.request.urlopen(url,timeout=120) as response:
 size=int(response.headers['Content-Length']);headers={k:response.headers.get(k) for k in ('Content-Length','Last-Modified','ETag','Content-Type')};raw=Reader(response);members=[];matches=[]
 with tarfile.open(fileobj=raw,mode='r|bz2') as archive:
  for member in archive:
   if not member.isfile():continue
   members.append({'name':member.name,'size':member.size})
   if len(members)<=3:print('member',member.name,member.size,flush=True)
   handle=archive.extractfile(member);header=None;number=0
   for line in handle:
    number+=1
    if number==1:header=line.decode('utf-8',errors='replace').rstrip('\n')
    if b'ENSG00000239268' in line or b'LINC03051' in line:
     matches.append({'member':member.name,'line':number,'header':header,'row':line.decode('utf-8',errors='replace').rstrip('\n')})
   if len(members)%100==0:print('members read',len(members),'matches',len(matches),flush=True)
 while raw.read(1024*1024):pass
 assert raw.n==size
 receipt={'url':url,'headers':headers,'bytes_read':raw.n,'sha256':raw.h.hexdigest(),
  'publisher_checksum':'not identified','regular_members':members,'matches':matches,
  'note':'Full transport stream hashed; archive not retained; source authenticity and cohort overlap not independently validated.'}
 (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('completed',len(members),'members',len(matches),'matching rows',flush=True)
