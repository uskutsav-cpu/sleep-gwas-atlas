#!/usr/bin/env python3
"""Independent deletion-mask lookup for any one-base substitutions."""
import csv,gzip,json,tarfile
from collections import defaultdict
from pathlib import Path
root=Path('/Volumes/Extreme SSD/brain6-crossdisorder-mpra'); out=Path(__file__).resolve().parent
ids={'rs185260','rs325485','rs4916723','rs62369151'}; masks=defaultdict(set); mapped=set()
with gzip.open(root/'unique_bcdat.tsv.gz','rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  if row['rsid'] not in ids:continue
  b=row['barcode'];assert len(b)==20 and set(b)<=set('ACGT');mapped.add(b)
for b in mapped:
 for i in range(20):masks[(i,b[:i]+b[i+1:])].add(b)
hits=set();n=0
with tarfile.open(root/'GSE244011_RAW.tar') as tar:
 for member in tar.getmembers():
  with gzip.open(tar.extractfile(member),'rt') as f:
   for row in csv.DictReader(f):
    b=row['barcode'];n+=1
    if len(b)!=20:continue
    for i in range(20):
     for original in masks.get((i,b[:i]+b[i+1:]),()):
      if b[i]!=original[i]:hits.add((b,original))
result={'unique_labels_all_id_rows':len(mapped),'raw_rows':n,'one_substitution_pairs':len(hits),
        'note':'All mapping rows with the four rsIDs form a superset of the coordinate/allele-restricted scan.'}
assert not hits
(out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
