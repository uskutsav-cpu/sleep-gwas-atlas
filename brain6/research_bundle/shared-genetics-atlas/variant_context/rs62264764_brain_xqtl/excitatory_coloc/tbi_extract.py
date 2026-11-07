import gzip,struct,zlib,hashlib,csv,json,subprocess,math
from pathlib import Path
import pysam
p=Path(__file__).resolve().parent; base=Path('/Volumes/Brain6Work/xqtl-chr3');lo,hi=117423157,118423158
url='https://st1.niagads.org/portal/v1/download-public/NG00184/fileset/57103';exc_tar=198894418432;exc_size=671965359
sha=lambda b:hashlib.sha256(b).hexdigest()
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
def decode(data):
 out=bytearray();i=0;blocks=0
 while i+18<=len(data):
  assert data[i:i+3]==b'\x1f\x8b\x08' and data[i+12:i+14]==b'BC'
  size=struct.unpack_from('<H',data,i+16)[0]+1
  if i+size>len(data):break
  out.extend(zlib.decompress(data[i:i+size],wbits=31));i+=size;blocks+=1
 assert blocks>0
 return out.decode(),blocks,i
def filter_lines(text):
 out=[]
 for line in text.splitlines():
  x=line.split('\t')
  if len(x)!=20 or x[0]!='chr3':continue
  try:start=int(x[1]);end=int(x[2])
  except ValueError:continue
  if start<hi and end>lo:out.append(line)
 return out
# Validate whole-region extraction using a fully authenticated local file.
a,b,meta=spans((base/'oli.bed.gz.tbi').read_bytes(),'chr3',lo,hi)
with (base/'oli.bed.gz').open('rb') as f:f.seek(a);buf=f.read(b-a+65536)
text,nblocks,used=decode(buf);test=filter_lines(text);tab=pysam.TabixFile(str(base/'oli.bed.gz'));reference=list(tab.fetch('chr3',lo,hi));tab.close();assert test==reference,(len(test),len(reference))
# Index must match publisher. Header validation ties offset to archive member name and size.
index=(base/'exc.bed.gz.tbi').read_bytes();assert hashlib.md5(index).hexdigest()=='1575ac67da022cead166bc28eaf3ca0d'
h=(Path('work/ng_chr3_exc_header.bin')).read_bytes();name=h[512:1024].split(b'\0')[0].decode();assert '.exc.' in name and '_3_snp_all.' in name and int(h[1148:1160].strip(b'\0 '),8)==exc_size
header=(Path('work/ng_chr3_exc_index_header.bin')).read_bytes();assert header[512:1024].split(b'\0')[0].decode()==name+'.tbi' and int(header[1148:1160].strip(b'\0 '),8)==len(index)
a,b,idxmeta=spans(index,'chr3',lo,hi);end=min(exc_size-1,b+65535);assert 0<=a<=end<exc_size and end-a+1<100_000_000
out=base/'exc_region.bgzf';headers=base/'exc_range_headers.txt'
subprocess.run(['curl','-fsSL','--max-time','180','-r',f'{exc_tar+a}-{exc_tar+end}','-D',str(headers),url,'-o',str(out)],check=True)
buf=out.read_bytes();assert len(buf)==end-a+1
hdr=headers.read_text();assert f'Content-Range: bytes {exc_tar+a}-{exc_tar+end}/' in hdr or f'content-range: bytes {exc_tar+a}-{exc_tar+end}/' in hdr
text,nblocks,used=decode(buf);extracted=filter_lines(text);assert extracted
# Keep complete regional rows in SSD; analysis checks exact marker and target below.
regional=base/'exc_region.tsv'
with regional.open('w') as f:
 for line in extracted:f.write(line+'\n')
ref={'input_index_md5':hashlib.md5(index).hexdigest(),'index_sha256':sha(index),'tar_member':name,'tar_member_size':exc_size,'tar_member_data_start':exc_tar,'file_byte_range':[a,end],'archive_byte_range':[exc_tar+a,exc_tar+end],'range_bytes':len(buf),'range_sha256':sha(buf),'complete_bgzf_blocks':nblocks,'complete_block_bytes':used,'region_rows':len(extracted),'region_sha256':sha(regional.read_bytes()),'local_validation':{'control':'complete checksum-verified oligodendrocyte source','region_lines_equal_tabix':True,'region_rows':len(reference)},'warning':'Excitatory full member and archive not downloaded or whole-file hashed; selected range has no publisher checksum.'}
(p/'retrieval_receipt.json').write_text(json.dumps(ref,indent=2)+'\n');print(json.dumps(ref,indent=2))
