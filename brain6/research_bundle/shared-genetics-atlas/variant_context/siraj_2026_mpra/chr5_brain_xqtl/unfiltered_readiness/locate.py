import urllib.request,json,hashlib,time
from pathlib import Path
p=Path(__file__).resolve().parent
url='https://st1.niagads.org/portal/v1/download-public/NG00184/fileset/57073';total=1471550000000
request=urllib.request.Request(url,headers={'Range':'bytes=0-2047','Accept-Encoding':'identity'});resp=urllib.request.urlopen(request,timeout=60);signed=resp.geturl();first=resp.read();assert resp.status==206 and len(first)==2048

def get(a,n=2048):
 req=urllib.request.Request(signed,headers={'Range':f'bytes={a}-{a+n-1}','Accept-Encoding':'identity'});r=urllib.request.urlopen(req,timeout=60);b=r.read();assert r.status==206 and len(b)==n,(a,r.status,len(b));return b
def parse_size(b):
 if b[0]&128:
  x=bytearray(b);x[0]&=127;return int.from_bytes(x,'big')
 return int(b.strip(b'\0 ') or b'0',8)
wanted='FunGenADv1.Knight-ADRC.mQTL.PC.NGFGXQTL0002_v1_5_snp_all.20260330.bed.gz';offset=0;found=[];n=0
while offset<total and n<500:
 b=first if offset==0 else get(offset);h=b[:512]
 if h==bytes(512):break
 typ=h[156:157];sz=parse_size(h[124:136])
 if typ==b'L':
  name=b[512:512+sz].split(b'\0')[0].decode();hpos=offset+512+((sz+511)//512)*512
  if hpos+512>offset+len(b):h=get(hpos,512)
  else:h=b[hpos-offset:hpos-offset+512]
  typ2=h[156:157];assert typ2 in (b'0',b'\0'),(offset,name,typ2)
  size=parse_size(h[124:136]);start=hpos+512
 else:
  name=h[:100].split(b'\0')[0].decode();size=sz;start=offset+512
 n+=1
 if wanted in name:
  found.append({'member':name,'header_offset':hpos if typ==b'L' else offset,'data_start':start,'size':size})
  print('found',name,size,flush=True)
 if n%20==0:print('examined',n,'current',name[:65],flush=True)
 if len(found)==2:break
 offset=start+((size+511)//512)*512
assert len(found)==2,(n,found)
index=next(x for x in found if x['member'].endswith('.tbi'));data=next(x for x in found if x['member'].endswith('.bed.gz'))
b=get(index['data_start'],index['size']);md5=hashlib.md5(b).hexdigest();assert md5=='57fe8e4233f4a4ecea73bb321a20bf2e',md5
(p/'chr5_index.tbi').write_bytes(b)
receipt={'catalog_public_url':url,'archive_size':total,'headers_examined':n,'members':found,'index_md5':md5,'full_archive_downloaded':False,'full_member_md5_verified':False,'signed_url_not_retained':True}
(p/'locate_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
