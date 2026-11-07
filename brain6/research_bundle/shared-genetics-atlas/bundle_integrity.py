#!/usr/bin/env python3
"""Inventory and verify the delivered atlas and pilot; omit live scan files."""
import argparse,csv,hashlib,pathlib,sys
atlas=pathlib.Path(__file__).resolve().parent;root=atlas.parent;manifest=atlas/'bundle_manifest.tsv'
def include(p):
 if not p.is_file() or p==manifest or '__pycache__' in p.parts or p.suffix == '.pyc':return False
 rel=p.relative_to(root).as_posix()
 if '/placo_diagnostic/' in rel and (p.name.endswith('.log') or '.partial.' in p.name):return False
 return True
def current():
 out={}
 for folder in ('shared-genetics-atlas','brain6-alphagenome-pilot'):
  for p in (root/folder).rglob('*'):
   if include(p):out[p.relative_to(root).as_posix()]=(str(p.stat().st_size),hashlib.sha256(p.read_bytes()).hexdigest())
 return out
ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();found=current()
if args.write:
 with manifest.open('w',newline='') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['path_relative_to_outputs','bytes','sha256']);w.writerows((k,*found[k]) for k in sorted(found))
 print('WROTE',len(found),'files',manifest)
else:
 with manifest.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 expected={r['path_relative_to_outputs']:(r['bytes'],r['sha256']) for r in rows};assert len(expected)==len(rows),'duplicate manifest paths'
 missing=sorted(expected.keys()-found.keys());added=sorted(found.keys()-expected.keys());changed=sorted(k for k in expected.keys()&found.keys() if expected[k]!=found[k]);print('checked',len(expected),'missing',len(missing),'unlisted',len(added),'changed',len(changed));
 for kind,values in [('MISSING',missing),('UNLISTED',added),('CHANGED',changed)]:
  for value in values:print(kind,value)
 if missing or added or changed:sys.exit(1)
