import csv,gzip,hashlib,io,json,ssl,struct,urllib.request,zlib
from pathlib import Path
p=Path(__file__).resolve().parent;url='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content';size=929286104;ctx=ssl.create_default_context(cafile='/etc/ssl/cert.pem')
def req(a,b):
 r=urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}','Accept-Encoding':'identity'}),context=ctx,timeout=90)
 with r: data=r.read();assert r.status==206 and len(data)==b-a+1;return data
tail=req(size-65536,size-1);idx=tail.rfind(b'PK\x05\x06');assert idx>=0;e=struct.unpack_from('<4s4H2LH',tail,idx);off=e[6]-(size-65536);entries=[]
for _ in range(e[4]):
 x=struct.unpack_from('<4s6H3I5H2I',tail,off);assert x[0]==b'PK\x01\x02';nl,xl,cl=x[10:13];name=tail[off+46:off+46+nl].decode();off+=46+nl+xl+cl
 if '/DE/OL33_' in name and any('_'+c+'_' in name for c in ('A549','HepG2','K562','GM12878')):entries.append(dict(name=name,method=x[4],crc=x[7],compressed_size=x[8],size=x[9],offset=x[16]))
assert len(entries)==4
rows=[];receipt=[]
for item in entries:
 off=item['offset'];h=req(off,off+29);lh=struct.unpack('<4s5H3I2H',h);assert lh[0]==b'PK\x03\x04';start=off+30+lh[9]+lh[10];raw=req(start,start+item['compressed_size']-1)
 data=zlib.decompress(raw,-15) if item['method']==8 else raw;assert len(data)==item['size'] and zlib.crc32(data)&0xffffffff==item['crc']
 count=0;matches=[]
 with gzip.open(io.BytesIO(data),'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   count+=1
   if row['chr'].removeprefix('chr')=='3' and row['pos']=='117637697' and row['ref_allele']=='T' and row['alt_allele']=='C':matches.append(row)
 assert len(matches)<=1,(item['name'],len(matches))
 for row in matches:rows.append(dict(source=item['name'],**row))
 receipt.append(dict(source=item['name'],compressed_size=item['compressed_size'],gzip_sha256=hashlib.sha256(data).hexdigest(),zip_crc_verified=True,rows=count,matches=len(matches)))
 print(item['name'],'rows',count,'matches',len(matches),flush=True)
if rows:
 with (p/'source_rows.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
(p/'receipt.json').write_text(json.dumps(dict(url=url,archive_size=size,whole_archive_hash_verified=False,query='GRCh37 chr3:117637697 T>C',files=receipt),indent=2)+'\n')
