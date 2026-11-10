#!/usr/bin/env python3
"""Extract a GRCh38 locus from SingleBrain's remotely hosted BGZF file using HTTP byte ranges."""
from pathlib import Path
import csv,gzip,hashlib,json,re,struct,time,sys
import requests
root=Path(__file__).resolve().parents[1];out=root/'outputs';cell=sys.argv[1] if len(sys.argv)>1 else 'Ext'; assert cell in {'Ext','IN','OD','OPC'}; record=requests.get('https://zenodo.org/api/records/14908182',timeout=30).json(); metadata=next(f for f in record['files'] if f['key']==f'{cell}_eqtl_full_assoc.tsv.gz');url=metadata['links']['self'];size=metadata['size']
pattern=re.compile(b'\x1f\x8b\x08\x04....\x00\xff\x06\x00BC\x02\x00',re.DOTALL);sess=requests.Session(); probes=[]
def get(start,end):
 end=min(end,size-1)
 for attempt in range(3):
  try:
   r=sess.get(url,headers={'Range':f'bytes={start}-{end}'},timeout=45);r.raise_for_status()
   if r.status_code!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {start}-{end}/'):raise ValueError('Server failed range request')
   return r.content
  except Exception:
   if attempt==2:raise
   time.sleep(attempt+1)
def parse_coords(line):
 f=line.split('\t')
 if len(f)<7:return None
 try:return (int(f[2].removeprefix('chr')),int(f[3]))
 except ValueError:return None
def probe(offset):
 offset=max(0,min(offset,size-131072)); b=get(offset,offset+131071)
 for m in pattern.finditer(b):
  i=m.start();n=struct.unpack('<H',b[i+16:i+18])[0]+1
  if i+n>len(b):continue
  try: lines=gzip.decompress(b[i:i+n]).decode().split('\n')[1:]
  except Exception:continue
  for line in lines:
   c=parse_coords(line)
   if c:
    probes.append({'requested':offset,'block':offset+i,'chr':c[0],'pos':c[1]});return offset+i,c
 raise ValueError('No parseable BGZF block')
def lower(target,lo,hi):
 while hi-lo>131072:
  mid=(lo+hi)//2;block,coord=probe(mid)
  if coord<target:lo=mid
  else:hi=mid
 return lo,hi
# Anchors verified independently by probing genome order: 20% chr3 and 35% chr6.
lo=int(size*.20);hi=int(size*.35); assert probe(lo)[1] < (5,135000000);assert probe(hi)[1] > (5,137300000)
a=lower((5,135000000),lo,hi);b=lower((5,137300000),a[0],hi)
start=max(0,a[0]-131072);end=min(size-1,b[1]+131072); blob=get(start,end); rows=[];coords=[];valid=0; first=None;last=None
for m in pattern.finditer(blob):
 i=m.start()
 if i+18>len(blob):continue
 n=struct.unpack('<H',blob[i+16:i+18])[0]+1
 if i+n>len(blob):continue
 try: raw=gzip.decompress(blob[i:i+n]).decode()
 except Exception:continue
 valid+=1
 # An individual BGZF block can start/end mid-line; join after decompression below.
 if first is None:first=i
 last=i+n
# The slice contains contiguous validated blocks; gzip supports concatenated BGZF members.
if first is None or last is None:raise ValueError('No complete blocks')
text=gzip.decompress(blob[first:last]).decode();lines=text.split('\n')[1:-1];previous=None;backwards=0
near_before=None;near_after=None
for line in lines:
 f=line.split('\t');c=parse_coords(line)
 if not c:continue
 if previous and c<previous:backwards+=1
 previous=c
 if c[0]==5 and c[1]<=136202653:near_before=c[1]
 if c[0]==5 and c[1]>136202653 and near_after is None:near_after=c[1]
 if c[0]==5 and 135000000<=c[1]<=137300000 and f[0].split('.')[0]=='ENSG00000113658':rows.append(f)
assert backwards==0, f'input not sorted: {backwards}'
assert rows, 'No SMAD5 rows'
header=['feature','variant_id','chr','pos','ref','alt','Allele','beta_tissue_0','sd_tissue_0','z_tissue_0','beta_tissue_1','sd_tissue_1','z_tissue_1','beta_tissue_2','sd_tissue_2','z_tissue_2','beta_tissue_3','sd_tissue_3','z_tissue_3','fixed_beta','fixed_sd','fixed_z','Random_Z','Fixed_P','Random_P','Fixed_bonf','Random_bonf','Fixed_FDR','Random_FDR']
assert all(len(r)==len(header) for r in rows)
with gzip.open(out/f'stage5_singlebrain_{cell.lower()}_smad5_full_regional.tsv.gz','wt') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)
receipt={'url':url,'source_size_bytes':size,'source_whole_file_md5_published':metadata['checksum'],'whole_file_checksum_verified':False,'requested_byte_range':[start,end],'range_sha256':hashlib.sha256(blob).hexdigest(),'complete_bgzf_blocks':valid,'source_order_backwards':backwards,'gene':'ENSG00000113658','chr5_region':[135000000,137300000],'retained_rows':len(rows),'all_genes_last_position_lead_or_before':near_before,'all_genes_first_position_after_lead':near_after,'first_retained':rows[0][2:4],'last_retained':rows[-1][2:4],'probes':probes}
(out/f'stage5_singlebrain_{cell.lower()}_smad5_regional_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['requested_byte_range','complete_bgzf_blocks','retained_rows','first_retained','last_retained','all_genes_last_position_lead_or_before','all_genes_first_position_after_lead']},indent=2))
