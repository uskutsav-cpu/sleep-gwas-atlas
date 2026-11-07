#!/usr/bin/env python3
# FAILED ATTEMPT: does not support concatenated BGZF gzip members; do not run.
import csv, hashlib, json, time, subprocess, zlib
from pathlib import Path
root=Path(__file__).resolve().parent
atlas=root.parent.parent
url='https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_ADHD.gz?generation=1777989549884404'
with (atlas/'replication/six_lead_joint/six_leads.tsv').open() as f:
 leads=list(csv.DictReader(f,delimiter='\t'))
targets={(x['chr'],x['position38']):x['variant'] for x in leads}
targets[('3','117920728')]='rs62264764'
hashmd5=hashlib.md5(); hashsha=hashlib.sha256(); decomp=zlib.decompressobj(16+zlib.MAX_WBITS)
rows=[]; carry=b''; nbytes=0; nrows=0; header=None
process=subprocess.Popen(['curl','--fail','--silent','--show-error','--location','--dump-header',str(root/'response_headers.txt'),url], stdout=subprocess.PIPE)
with process.stdout as response:
 while True:
  chunk=response.read(128*1024)
  if not chunk:break
  nbytes+=len(chunk);hashmd5.update(chunk);hashsha.update(chunk)
  data=carry+decomp.decompress(chunk); lines=data.split(b'\n');carry=lines.pop()
  for line in lines:
   if header is None:
    header=line.decode().split('\t');continue
   parts=line.decode().split('\t');nrows+=1
   if (parts[0],parts[1]) in targets: rows.append(dict(zip(header,parts)))
if process.wait()!=0: raise RuntimeError('curl failed')
headers=(root/'response_headers.txt').read_text()
if not decomp.eof or decomp.unused_data:raise RuntimeError('Incomplete gzip or trailing data')
if carry:
 parts=carry.decode().split('\t');nrows+=1
 if len(parts)!=len(header):raise RuntimeError('Malformed final record')
 if (parts[0],parts[1]) in targets:rows.append(dict(zip(header,parts)))
md5=hashmd5.hexdigest();sha=hashsha.hexdigest()
if nbytes!=802195177 or md5!='e7aedc4602071840bbd5abe62d6fcdcc':raise RuntimeError('Frozen source checksum or byte count mismatch')
with (root/'source_rows.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=header,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
(root/'receipt.json').write_text(json.dumps({'url':url,'response_headers':headers,'compressed_bytes':nbytes,'md5':md5,'sha256':sha,'gzip_eof':decomp.eof,'row_count':nrows,'target_coordinates':targets and [{'chrom':k[0],'pos':k[1],'expected_rsid':v} for k,v in targets.items()],'rows_at_target_coordinates':len(rows),'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())},indent=2)+'\n')
print('complete',nbytes,'rows',nrows,'target rows',len(rows),'md5',md5,flush=True)
