import csv,struct,urllib.request,zipfile,zlib,gzip,io,hashlib,json,re
from pathlib import Path
p=Path(__file__).resolve().parent;ssd=Path('/Volumes/Extreme SSD/brain6-large-mpra');src=p.parents[1]/'variant_context/rs40465_mpra/ranked_screen.tsv';url='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content';size=929286104
with src.open() as f:rank=list(csv.DictReader(f,delimiter='\t'))
wanted={r['grch37']+':'+r['grch38'].split(':')[-1].replace('>',':'):r for r in rank if r.get('grch37')}
def req(start,end):
 r=urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'}),timeout=90)
 b=r.read();assert r.status==206 and len(b)==end-start+1,(r.status,len(b),start,end);return b
buf=req(size-65536,size-1);idx=buf.rfind(b'PK\x05\x06');e=struct.unpack_from('<4s4H2LH',buf,idx);offset=e[6]-(size-65536);entries=[]
for _ in range(e[4]):
 x=struct.unpack_from('<4s6H3I5H2I',buf,offset);assert x[0]==b'PK\x01\x02';nl,xl,cl=x[10:13];name=buf[offset+46:offset+46+nl].decode();offset+=46+nl+xl+cl
 if '/DE/' in name and '_SKNSH_' in name:entries.append({'name':name,'method':x[4],'crc':x[7],'compressed_size':x[8],'size':x[9],'offset':x[16]})
assert len(entries)==9,len(entries)
records=[];receipt=[]
for item in entries:
 off=item['offset'];h=req(off,off+29);lh=struct.unpack('<4s5H3I2H',h);assert lh[0]==b'PK\x03\x04';start=off+30+lh[9]+lh[10];raw=req(start,start+item['compressed_size']-1)
 data=zlib.decompress(raw,-15) if item['method']==8 else raw
 assert len(data)==item['size'] and zlib.crc32(data)&0xffffffff==item['crc']
 with gzip.open(io.BytesIO(data),'rt') as f:
  reader=csv.DictReader(f,delimiter='\t');header=reader.fieldnames;n=0
  for row in reader:
   n+=1
   chrom=row.get('chr','');chrom=chrom if chrom.startswith('chr') else 'chr'+chrom
   v=f"{chrom}:{row.get('pos','')}:{row.get('ref_allele','')}:{row.get('alt_allele','')}"
   if v in wanted:records.append({'source':item['name'],'query_variant':v,'rsid':wanted[v]['rsid'],'row':row})
 receipt.append({'file':item['name'],'compressed_zip_bytes':item['compressed_size'],'zip_crc_verified':True,'gzip_sha256':hashlib.sha256(data).hexdigest(),'data_rows':n,'headers':header});print(item['name'],n,flush=True)
(p/'sknsh_matches.json').write_text(json.dumps(records,indent=2)+'\n');(p/'sknsh_receipt.json').write_text(json.dumps({'url':url,'archive_size':size,'full_archive_checksum_verified':False,'files':receipt,'matches':len(records)},indent=2)+'\n');print('matches',len(records))
