#!/usr/bin/env python3
"""Stream publisher tables, verify publisher MD5, and retain one gene's lines."""
import hashlib,json,urllib.request,zlib
from pathlib import Path
p=Path(__file__).resolve().parent
api='https://zenodo.org/api/records/19924423'
with urllib.request.urlopen(api,timeout=60) as r: meta=json.load(r)
wanted=['twas_hits.gtex-ddp.tsv.gz','twas_hits.gtex-rddp.tsv.gz','twas_hits.gtex-kdp.tsv.gz','qtls.gtex-ddp.tsv.gz']
files={f['key']:f for f in meta['files']};receipt={}
for name in wanted:
 f=files[name];md5=hashlib.md5();sha=hashlib.sha256();decomp=zlib.decompressobj(31);pending=b'';header=None;selected=[];total=0;nbytes=0
 with urllib.request.urlopen(f['links']['self'],timeout=120) as source:
  while chunk:=source.read(1024*1024):
   nbytes+=len(chunk);md5.update(chunk);sha.update(chunk);pending+=decomp.decompress(chunk)
   lines=pending.split(b'\n');pending=lines.pop()
   for line in lines:
    if header is None:header=line;continue
    total+=1
    if b'ENSG00000239268' in line or b'LINC03051' in line:selected.append(line)
  pending+=decomp.flush()
  if pending:
   total+=1
   if b'ENSG00000239268' in pending or b'LINC03051' in pending:selected.append(pending)
 assert decomp.eof and nbytes==f['size'] and f['checksum']=='md5:'+md5.hexdigest()
 target=p/(name.replace('.tsv.gz','')+'_selected.tsv')
 target.write_bytes(header+b'\n'+b''.join(line+b'\n' for line in selected))
 receipt[name]={'url':f['links']['self'],'size':nbytes,'publisher_checksum':f['checksum'],'sha256':sha.hexdigest(),
  'rows_scanned':total,'rows_selected':len(selected),'selected_sha256':hashlib.sha256(target.read_bytes()).hexdigest()}
 print(name,total,len(selected),flush=True)
(p/'receipt.json').write_text(json.dumps({'record':api,'files':receipt,'note':'Publisher hit tables are prefiltered; selection by gene identifier text.'},indent=2)+'\n')
