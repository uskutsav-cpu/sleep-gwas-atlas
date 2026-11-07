#!/usr/bin/env python3
"""Verify selected provenance-linked bundle entries, separately from bundle_integrity.py."""
import csv,hashlib,pathlib,sys
root=pathlib.Path(__file__).resolve().parent;manifest=root/'manifest.tsv'
with manifest.open(newline='') as f:
 reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;rows=list(reader)
assert fields==['file','bytes','sha256','source_path','source_sha256_matches','note']
problems=[]
for row in rows:
 path=root/row['file']
 if not path.is_file():problems.append((row['file'],'missing'));continue
 content=path.read_bytes();actual=hashlib.sha256(content).hexdigest()
 if str(len(content))!=row['bytes'] or actual!=row['sha256']:problems.append((row['file'],'size/hash changed'))
print(f'selected manifest: {len(rows)} entries, {len(problems)} mismatches')
for file,problem in problems:print(file,problem,sep='\t')
sys.exit(bool(problems))
