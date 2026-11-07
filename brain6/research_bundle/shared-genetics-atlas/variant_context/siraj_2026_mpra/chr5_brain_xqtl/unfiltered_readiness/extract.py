import gzip,struct,zlib,hashlib,csv,json,urllib.request,io
from pathlib import Path
p=Path(__file__).resolve().parent
def spans(idx,contig,beg,end):
 data=gzip.decompress(idx);assert data[:4]==b'TBI\x01';off=4
 def read(fmt):
  nonlocal off
  v=struct.unpack_from('<'+fmt,data,off);off+=struct.calcsize('<'+fmt);return v[0] if len(v)==1 else v
 nr,fmt,cs,cb,ce,meta,skip,ln=read('8i'); names=data[off:off+ln].split(b'\0');off+=ln
 names=[x.decode() for x in names if x];assert len(names)==nr and contig in names
 bins={};linear=[]
 for j in range(nr):
  nb=read('i');d={}
  for _ in range(nb):
   binno,nc=read('Ii');d[binno]=[read('QQ') for _ in range(nc)]
  ni=read('i');lin=[read('Q') for _ in range(ni)]
  if names[j]==contig:bins=d;linear=lin
 end-=1;query={0}
 for shift,start in [(26,1),(23,9),(20,73),(17,585),(14,4681)]:query.update(range(start+(beg>>shift),start+(end>>shift)+1))
 candidates=[chunk for b in query for chunk in bins.get(b,[]) if chunk[1]>chunk[0]];assert candidates
 minoffset=linear[beg>>14] if beg>>14<len(linear) else 0
 candidates=[x for x in candidates if x[1]>minoffset];assert candidates
 a=min(max(x[0],minoffset) for x in candidates)>>16;b=max(x[1] for x in candidates)>>16
 return a,b,{'total_bins':len(bins),'chunks':len(candidates),'min_virtual_offset':minoffset}

meta=json.loads((p/'locate_receipt.json').read_text());entry=next(x for x in meta['members'] if x['member'].endswith('.bed.gz'));index=(p/'chr5_index.tbi').read_bytes();assert hashlib.md5(index).hexdigest()=='57fe8e4233f4a4ecea73bb321a20bf2e'
lo,hi=104112266,105112267
a,b,info=spans(index,'chr5',lo,hi);end=min(entry['size']-1,b+65535);start_tar=entry['data_start'];assert end-a+1<40000000
req=urllib.request.Request(meta['catalog_public_url'],headers={'Range':f'bytes={start_tar+a}-{start_tar+end}','Accept-Encoding':'identity'});r=urllib.request.urlopen(req,timeout=120);buf=r.read();assert r.status==206 and len(buf)==end-a+1
sha=hashlib.sha256(buf).hexdigest();print('compressed range',len(buf),'sha256',sha,flush=True)
pos=0;blocks=0;tail=b'';found=[];total=0;partial=True
while pos+18<=len(buf):
 assert buf[pos:pos+3]==b'\x1f\x8b\x08' and buf[pos+12:pos+14]==b'BC'
 length=struct.unpack_from('<H',buf,pos+16)[0]+1
 if pos+length>len(buf):break
 raw=zlib.decompress(buf[pos:pos+length],wbits=31);blocks+=1;pos+=length
 pieces=(tail+raw).split(b'\n');tail=pieces.pop()
 if partial:pieces=pieces[1:];partial=False
 for ln in pieces:
  total+=1
  if b'cg24905254' not in ln:continue
  x=ln.decode().split('\t');assert len(x)==20
  if x[0]=='chr5' and int(x[1])<hi and int(x[2])>lo:found.append(x)
assert found
with (p/'selected_cpg.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerows(found)
receipt={'region_grch38_bed':[lo,hi],'member':entry['member'],'member_size':entry['size'],'archive_data_start':start_tar,'compressed_range':[a,end],'range_bytes':len(buf),'range_sha256':sha,'complete_bgzf_blocks':blocks,'complete_block_bytes':pos,'rows_in_blocks':total,'selected_cpg_rows':len(found),'index':info,'index_publisher_md5_verified':True,'full_member_md5_verified':False,'archive_md5_verified':False,'selected_output_sha256':hashlib.sha256((p/'selected_cpg.tsv').read_bytes()).hexdigest(),'notes':'range reconstructed from public index; first partial line omitted; no external checksum for range; selected CpG only'}
(p/'extract_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
