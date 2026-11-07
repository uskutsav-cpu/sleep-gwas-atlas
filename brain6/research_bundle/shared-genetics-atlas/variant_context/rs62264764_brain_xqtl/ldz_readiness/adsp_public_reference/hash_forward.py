#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path
p=Path(__file__).resolve().parent;md=hashlib.md5();sha=hashlib.sha256();n=0;step=1024**3
while b:=sys.stdin.buffer.read(1024*1024):
 n+=len(b);md.update(b);sha.update(b);sys.stdout.buffer.write(b)
 if n>=step:print(f'streamed {n//(1024**3)} GiB',file=sys.stderr,flush=True);step+=1024**3
sys.stdout.buffer.flush()
r={'compressed_bytes':n,'md5':md.hexdigest(),'sha256':sha.hexdigest(),'expected_bytes':5043829695,'expected_md5':'481f8ba5dcb3cd30c8c31bacc675a3cf','source_file_id':55118}
if n!=r['expected_bytes'] or md.hexdigest()!=r['expected_md5']:raise RuntimeError('Source checksum or length mismatch')
(p/'stream_receipt.json').write_text(json.dumps(r,indent=2)+'\n')
