import json,subprocess
from pathlib import Path
p=Path(__file__).resolve().parent
url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90565001-GCST90566000/GCST90565837/GCST90565837.tsv'
size=713768673

def read(start,n=8192):
 end=min(size-1,start+n-1)
 raw=subprocess.check_output(['curl','-fsSL','--max-time','35','-D','-','--range',f'{start}-{end}',url])
 headers,data=raw.split(b'\r\n\r\n',1)
 assert b' 206 ' in headers.split(b'\r\n',1)[0] and f'bytes {start}-{end}/{size}'.encode() in headers
 assert len(data)==end-start+1
 return data
head=read(0).split(b'\n',1)[0].decode().rstrip('\r').split('\t')
assert head==['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','p_value','n']
def row_at(offset):
 b=read(offset)
 if offset: b=b.split(b'\n',1)[1]
 line=b.split(b'\n',1)[0].decode();row=line.split('\t')
 if offset==0:row=b.split(b'\n',2)[1].decode().split('\t')
 assert len(row)==len(head)
 return int(row[0]),int(row[1])
target=(5,103947968);lo=0;hi=size-1
while hi-lo>8192:
 mid=(lo+hi)//2;key=row_at(mid)
 if key<target:lo=mid
 else:hi=mid
chunk=read(max(0,lo-8192),min(524288,size-max(0,lo-8192)))
lines=chunk.split(b'\n');lines=lines[1:-1];rows=[];previous=None
for line in lines:
 fields=line.decode().rstrip('\r').split('\t')
 if len(fields)!=len(head):raise ValueError('invalid row in range')
 key=(int(fields[0]),int(fields[1]))
 if previous is not None and key<previous:raise ValueError('source unsorted in inspected range')
 previous=key
 if key in [(5,103947968),(5,103964585)]:rows.append(dict(zip(head,fields)))
assert {int(r['base_pair_location']) for r in rows}=={103947968,103964585}
out={'url':url,'scope':'partial indexed-by-byte read; no full-source checksum rechecked','scanned_start':max(0,lo-8192),'scanned_bytes':len(chunk),'rows':rows}
(p/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
