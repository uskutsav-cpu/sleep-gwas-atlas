#!/usr/bin/env python3
"""Combine complete, CRC-checked, ordered FFS chunks, preserving originals."""
import gzip,hashlib,json,pathlib,sys
HERE=pathlib.Path(__file__).resolve().parent
EXTERNAL=pathlib.Path('/Volumes/Extreme SSD/brain6-placo-ffs-rerun-20261027')
CHUNKS=[('FFS_1_569131.tsv.gz',569131),('FFS_569132_1138262.tsv.gz',569131)]

def combine(readers,writer,expected):
 header=None;counts=[]
 for f,n in zip(readers,expected,strict=True):
  current=f.readline()
  if not current or not current.endswith(b'\n'):raise ValueError('missing or incomplete header')
  if header is None:header=current;writer.write(header)
  elif current!=header:raise ValueError('chunk headers differ')
  count=0
  for line in f:
   if not line.endswith(b'\n'):raise ValueError('unterminated data row')
   writer.write(line);count+=1
  if count!=n:raise ValueError(f'wrong row count {count}, expected {n}')
  counts.append(count)
 return counts

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(4*1024*1024):h.update(b)
 return h.hexdigest()

def main():
 from contextlib import ExitStack
 paths=[EXTERNAL/n for n,_ in CHUNKS]; target=HERE/'FFS_screen.tsv.gz'; partial=HERE/'FFS_screen.partial.tsv.gz'
 if target.exists() or partial.exists():raise FileExistsError('preserve existing completed or partial output')
 if not all(p.is_file() for p in paths):raise FileNotFoundError('missing completed chunk; no partial chunks admitted')
 with ExitStack() as stack:
  inp=[stack.enter_context(gzip.open(p,'rb')) for p in paths]
  raw=stack.enter_context(partial.open('xb'));out=stack.enter_context(gzip.GzipFile(fileobj=raw,mode='wb',mtime=0))
  counts=combine(inp,out,[n for _,n in CHUNKS]) # reading gzip to EOF enforces CRC
 if sum(counts)!=1138262:raise ValueError('combined family support mismatch')
 # Read newly written gzip to EOF to check its own CRC before rename.
 with gzip.open(partial,'rb') as f:
  while f.read(4*1024*1024):pass
 partial.rename(target)
 result={'input_chunks':[{'path':str(x),'sha256':sha(x),'expected_rows':n} for x,(_,n) in zip(paths,CHUNKS)],'rows_per_chunk':counts,'combined_path':str(target),'combined_sha256':sha(target),'crc_checks':'passed for both source chunks and combined output','header':'single header; both input headers checked equal','note':'Full aligned row/control audit remains required.'}
 (HERE/'merge_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
