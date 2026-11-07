import csv,gzip,hashlib,io,json,math,ssl,struct,urllib.request,zlib,collections
from pathlib import Path
p=Path(__file__).resolve().parent;url='https://zenodo.org/api/records/15297965/files/data_MPRA.zip/content';size=929286104;ctx=ssl.create_default_context(cafile='/etc/ssl/cert.pem')
def req(a,b):
 with urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={a}-{b}','Accept-Encoding':'identity'}),context=ctx,timeout=90) as r:
  data=r.read();assert r.status==206 and len(data)==b-a+1;return data
tail=req(size-65536,size-1);i=tail.rfind(b'PK\x05\x06');assert i>=0;e=struct.unpack_from('<4s4H2LH',tail,i);off=e[6]-(size-65536);entries=[]
for _ in range(e[4]):
 x=struct.unpack_from('<4s6H3I5H2I',tail,off);assert x[0]==b'PK\x01\x02';nl,xl,cl=x[10:13];name=tail[off+46:off+46+nl].decode();off+=46+nl+xl+cl
 if name.endswith('/DE/OL33_SKNSH_emVAR_glm_paired_20221220.out.gz'):entries.append(dict(name=name,method=x[4],crc=x[7],csize=x[8],size=x[9],offset=x[16]))
assert len(entries)==1;entry=entries[0];off=entry['offset'];h=req(off,off+29);lh=struct.unpack('<4s5H3I2H',h);assert lh[0]==b'PK\x03\x04';start=off+30+lh[9]+lh[10];raw=req(start,start+entry['csize']-1);data=zlib.decompress(raw,-15) if entry['method']==8 else raw
assert len(data)==entry['size'] and zlib.crc32(data)&0xffffffff==entry['crc'];rows=[];total=0
def num(v):
 try:
  x=float(v);return x if math.isfinite(x) else None
 except (ValueError,TypeError):return None
with gzip.open(io.BytesIO(data),'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  total+=1;chrom=r['chr'].removeprefix('chr');pos=int(r['pos']);region='chr3' if chrom=='3' and 117142005<=pos<=118142005 else ('chr5' if chrom=='5' and 103447968<=pos<=104447968 else None)
  if not region or len(r['ref_allele'])!=1 or len(r['alt_allele'])!=1 or r['ref_allele'] not in 'ACGT' or r['alt_allele'] not in 'ACGT' or r['ref_allele']==r['alt_allele']:continue
  vals=[num(r[k]) for k in ('Log2Skew','Skew_SE','Skew_logFDR','A_log2FC','B_log2FC','A_logPadj_BF','B_logPadj_BF')]
  complete=all(x is not None for x in vals)
  if complete:
   skew,se,fdr,aa,bb,ap,bp=vals;assert se>=0
   active=(aa>0 and ap>-math.log10(.05)) or (bb>0 and bp>-math.log10(.05))
   hit=active and fdr>-math.log10(.05);near=active and abs(skew)+1.96*se<.25
  else:active=hit=near=None
  rows.append(dict(region=region,**r,proxy_active=active,proxy_hit=hit,proxy_near_zero=near))
assert rows
with (p/'regional_rows.tsv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={}
for region in ('chr3','chr5'):
 subset=[r for r in rows if r['region']==region];keys=[(r['chr'],r['pos'],r['ref_allele'],r['alt_allele']) for r in subset];n=collections.Counter(keys)
 summary[region]=dict(rows=len(subset),unique_variants=len(n),duplicated_keys=sum(v>1 for v in n.values()),complete=sum(r['proxy_active'] is not None for r in subset),proxy_active=sum(r['proxy_active'] is True for r in subset),proxy_hit=sum(r['proxy_hit'] is True for r in subset),proxy_near_zero=sum(r['proxy_near_zero'] is True for r in subset))
receipt=dict(url=url,member=entry['name'],member_gzip_sha256=hashlib.sha256(data).hexdigest(),zip_crc_verified=True,whole_archive_hashed=False,total_source_rows=total,summary=summary,regional_sha256=hashlib.sha256((p/'regional_rows.tsv').read_bytes()).hexdigest());(p/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(summary,indent=2))
