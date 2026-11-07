#!/usr/bin/env python3
"""Hash exact pinned compressed object without decompressing or storing its bytes."""
import hashlib,json,subprocess,sys,time
from pathlib import Path
p=Path(__file__).resolve().parent
url='https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_ADHD.gz?generation=1777989549884404'
proc=subprocess.Popen(['curl','--fail','--silent','--show-error','--location',url],stdout=subprocess.PIPE)
md5=hashlib.md5();sha=hashlib.sha256();size=0;next_notice=128*1024*1024
with proc.stdout as stream:
 while chunk:=stream.read(256*1024):
  size+=len(chunk);md5.update(chunk);sha.update(chunk)
  if size>=next_notice:
   print(f'Hashed {size//(1024*1024)} MiB',flush=True);next_notice+=128*1024*1024
if proc.wait()!=0:raise RuntimeError('Source download failed')
result={'url':url,'compressed_bytes':size,'md5':md5.hexdigest(),'sha256':sha.hexdigest(),'whole_file_bytes_verified_against_historical_md5':size==802195177 and md5.hexdigest()=='e7aedc4602071840bbd5abe62d6fcdcc','compressed_gzip_integrity_tested':False,'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
if not result['whole_file_bytes_verified_against_historical_md5']:raise RuntimeError('Length or MD5 mismatch: '+json.dumps(result))
(p/'whole_source_hash_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print('Complete:',json.dumps(result),flush=True)
