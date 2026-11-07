import csv,json,gzip,collections,hashlib
from pathlib import Path
p=Path(__file__).resolve().parent;r=json.loads((p/'retrieval_receipt.json').read_text());scan=json.loads((p/'scan_receipt.json').read_text());expected={104612267:'rs2431108',104628884:'rs77960',104655775:'rs416223'};observed=collections.Counter()
for src in r['members']:
 f=Path(src['file']);h=hashlib.sha256()
 with f.open('rb') as x:
  for b in iter(lambda:x.read(1024*1024),b''):h.update(b)
 assert h.hexdigest()==src['inner_gzip_sha256']
 with gzip.open(f,'rt') as x:
  next(x)
  for line in x:
   v=line.split(None,1)[0].split(':');pos=int(v[1]);
   if pos in expected:observed[(src['cell'],expected[pos])]+=1
with (p/'candidate_rows.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==26 and collections.Counter((x['cell'],x['rsid']) for x in rows)==observed
assert min(float(x['p']) for x in rows)==0.04484256479
assert all(abs(float(x['bonferroni_observed_candidate_rows'])-min(1.,float(x['p'])*26))<1e-12 for x in rows)
(p/'verification.json').write_text(json.dumps({'rows':len(rows),'source_sha256_confirmed':True,'independent_coordinate_counts':dict((f'{k[0]}:{k[1]}',v) for k,v in observed.items()),'bonferroni_x26_arithmetic_verified':True,'minimum_source_p':min(float(x['p']) for x in rows),'note':'Data retrieval and arithmetic checks; not independent biological replication or calibration.'},indent=2)+'\n');print('verified 26 rows, source hashes, coordinate counts and correction')
