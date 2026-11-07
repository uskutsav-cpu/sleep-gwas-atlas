#!/usr/bin/env python3
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parent; base=root.parent
sources=[(base/'input_audit/download_receipt.json',Path('/Volumes/Brain6Sim/nhanes-public-2011-14-20261028')),
         (base/'item_mapping/download_receipt.json',Path('/Volumes/Brain6NHANES/source'))]
items=[]
for receipt,folder in sources:
 for r in json.loads(receipt.read_text()):
  path=folder/r['file'];h=hashlib.sha256()
  with path.open('rb') as f:
   for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
  items.append({'path':str(path),'expected':r['sha256'],'observed':h.hexdigest(),'match':h.hexdigest()==r['sha256']})
assert len(items)==32 and all(r['match'] for r in items)
(root/'input_verification.json').write_text(json.dumps({'inputs':len(items),'all_match':True,'files':items},indent=2)+'\n')
print('PASS: all 32 downloaded files match their saved SHA-256 receipts.')
