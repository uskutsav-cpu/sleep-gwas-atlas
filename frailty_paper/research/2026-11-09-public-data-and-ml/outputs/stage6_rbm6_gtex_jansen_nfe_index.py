#!/usr/bin/env python3
"""Decode selected public Hail 0.2.22 variant-index partitions, validating every record count and index."""
from pathlib import Path
import requests,gzip,json,struct,hashlib,concurrent.futures
import lz4.block,pandas as pd
from pyliftover import LiftOver
root=Path(__file__).resolve().parents[1];out=root/'outputs';work=root/'work';base='https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/ld/gnomad.genomes.r2.1.1.nfe.common.adj.ld.variant_indices.ht/';sess=requests.Session()
def fetch(path):r=sess.get(base+path,timeout=60);r.raise_for_status();return r.content
root_raw=fetch('metadata.json.gz');rows_raw=fetch('rows/metadata.json.gz');meta=json.loads(gzip.decompress(root_raw));rowmeta=json.loads(gzip.decompress(rows_raw));counts=meta['components']['partition_counts']['counts'];assert len(counts)==len(rowmeta['partFiles'])==len(rowmeta['jRangeBounds']);assert rowmeta['codecSpec']['name']=='PackCodecSpec';
d=pd.read_csv(work/'stage6_rbm6_gtex_jansen_input.csv');lift=LiftOver(str(work/'hg38ToHg19.over.chain.gz'));bykey={}
for r in d.itertuples():
 m=lift.convert_coordinate('chr3',int(r.pos)-1)
 if len(m)!=1 or m[0][0]!='chr3' or m[0][2]!='+':continue
 key=(m[0][1]+1,tuple(sorted([r.anchor_effect_allele,r.anchor_other_allele])))
 if key in bykey:raise ValueError('duplicate allele key')
 bykey[key]=r
poslist=[k[0] for k in bykey];lo=min(poslist);hi=max(poslist);parts=[]
for i,x in enumerate(rowmeta['jRangeBounds']):
 a=x['start']['locus'];b=x['end']['locus']
 if a['contig']=='3' and b['contig']=='3' and a['position']<=hi and b['position']>=lo:parts.append(i)

def decode(part):
 raw=fetch('rows/parts/'+rowmeta['partFiles'][part]);off=0;chunks=[]
 while off<len(raw):
  if len(raw)-off<4:raise ValueError('short block header')
  n=struct.unpack_from('<I',raw,off)[0];off+=4
  if n==0:break
  if n<4 or off+n>len(raw):raise ValueError('bad block size')
  sz=struct.unpack_from('<I',raw,off)[0];chunks.append(lz4.block.decompress(raw[off+4:off+n],uncompressed_size=sz));off+=n
 if off!=len(raw):raise ValueError('unexpected trailing bytes')
 s=b''.join(chunks);i=0;rows=[]
 def byte():
  nonlocal i
  if i>=len(s):raise ValueError('unexpected row EOF')
  x=s[i];i+=1;return x
 def leb():
  val=shift=0
  while True:
   b=byte();val|=(b&127)<<shift
   if b<128:return val
   shift+=7
   if shift>63:raise ValueError('invalid varint')
 def string():
  nonlocal i
  n=leb();v=s[i:i+n].decode('utf8');i+=n;return v
 while i<len(s):
  marker=byte()
  if marker==0:break
  if marker!=1 or byte()!=0:raise ValueError('unsupported row encoding or missing fields')
  chrom=string();pos=leb();n=leb();missing=byte()
  if missing!=0 or n!=2:raise ValueError('unsupported alleles')
  a=[string() for _ in range(n)]
  if byte()!=0:raise ValueError('missing population values')
  ac=leb();af=struct.unpack_from('<d',s,i)[0];i+=8;an=leb();hom=leb();idx=leb();rows.append((chrom,pos,a[0],a[1],ac,af,an,hom,idx))
 if i!=len(s) or len(rows)!=counts[part]:raise ValueError('partition completeness mismatch')
 expected=sum(counts[:part]);
 if any(row[-1]!=expected+k for k,row in enumerate(rows)):raise ValueError('matrix indices not contiguous')
 return part,rows,hashlib.sha256(raw).hexdigest()
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:got=list(ex.map(decode,parts))
matched=[]
for part,rows,sha in got:
 for chrom,pos,ref,alt,ac,af,an,hom,idx in rows:
  r=bykey.get((pos,tuple(sorted([ref,alt]))))
  if chrom!='3' or r is None:continue
  if ref==alt or r.anchor_effect_allele not in [ref,alt]:continue
  matched.append({'rs':r.rs,'pos37':pos,'ref37':ref,'alt37':alt,'gnomad_index':idx,'gnomad_af_alt':af,'gnomad_an':an,'anchor_effect_allele':r.anchor_effect_allele,'anchor_other_allele':r.anchor_other_allele,'sign_to_anchor':1 if r.anchor_effect_allele==alt else -1})
x=pd.DataFrame(matched);assert not x.rs.duplicated().any() and not x.gnomad_index.duplicated().any();x.to_csv(out/'stage6_rbm6_gtex_jansen_nfe_variant_index.csv',index=False)
receipt={'source':base,'root_metadata_sha256':hashlib.sha256(root_raw).hexdigest(),'rows_metadata_sha256':hashlib.sha256(rows_raw).hexdigest(),'GRCh37_min':lo,'GRCh37_max':hi,'partitions':[{'index':part,'name':rowmeta['partFiles'][part],'rows':len(rows),'sha256':sha} for part,rows,sha in got],'variants_input':len(d),'variants_matched':len(x),'validation':'All partition counts matched, idx sequential from official partition-count prefix, complete decoding, exact build-converted position and allele-pair matching. Custom decoder has not been checked using Hail.'};(out/'stage6_rbm6_gtex_jansen_nfe_index_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('parts',len(parts),'input',len(d),'matched',len(x),flush=True)
