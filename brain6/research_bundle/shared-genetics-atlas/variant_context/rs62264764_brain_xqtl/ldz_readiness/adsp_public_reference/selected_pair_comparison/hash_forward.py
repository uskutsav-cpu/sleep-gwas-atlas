#!/usr/bin/env python3
import hashlib,json,sys
from pathlib import Path
p=Path(__file__).resolve().parent;h=hashlib.sha256();m=hashlib.md5();n=0;step=1024**3
while b:=sys.stdin.buffer.read(1024*1024):
 n+=len(b);h.update(b);m.update(b);sys.stdout.buffer.write(b)
 if n>=step:print(f'streamed {n//1024**3} GiB',file=sys.stderr,flush=True);step+=1024**3
sys.stdout.buffer.flush()
assert n==5043829695 and m.hexdigest()=='481f8ba5dcb3cd30c8c31bacc675a3cf' and h.hexdigest()=='e7fbc170ffdba189d844174f1658adf9e088c529a29c9e4d407c92d24c01e4e2'
(p/'stream_receipt.json').write_text(json.dumps({'compressed_bytes':n,'md5':m.hexdigest(),'sha256':h.hexdigest(),'matches_manifest_and_previous_full_stream':True},indent=2)+'\n')
