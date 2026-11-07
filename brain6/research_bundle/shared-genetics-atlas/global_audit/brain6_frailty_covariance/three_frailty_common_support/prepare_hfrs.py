import csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
old=p.parents[1]/'insomnia_frailty_same_support/receipt.json'
prior=p.parent/'paired_common_support/prepare_receipt.json'
src=Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity/HFRS.sumstats.gz')
base=Path('/Volumes/Brain6Work/brain6-frailty-paired')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as x:
  for b in iter(lambda:x.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
assert sha(src)==json.loads(old.read_text())['source_sha256'][str(src)]
receipt=json.loads(prior.read_text())
for k in ('FI','FFS','insomnia'):assert sha(base/(k+'.sumstats.gz'))==receipt['output_sha256'][str(base/(k+'.sumstats.gz'))]
key=lambda r:(r['SNP'],r['A1'],r['A2'])
out=base/'HFRS.sumstats.gz'; count=0; skipped=0
with gzip.open(src,'rt') as a,gzip.open(base/'FI.sumstats.gz','rt') as b,gzip.open(out,'wt') as o:
 raw=csv.DictReader(a,delimiter='\t');wanted=csv.DictReader(b,delimiter='\t');w=csv.writer(o,delimiter='\t',lineterminator='\n');w.writerow(['SNP','A1','A2','Z','N'])
 for f in wanted:
  for h in raw:
   if key(h)==key(f):break
   skipped+=1
  else:raise ValueError('missing or out of order: '+str(key(f)))
  assert math.isfinite(float(h['Z'])) and math.isfinite(float(h['N'])) and float(h['N'])>0
  w.writerow([*key(h),h['Z'],h['N']]);count+=1
assert count==receipt['stats']['matched']
with gzip.open(out,'rt') as a,gzip.open(base/'FI.sumstats.gz','rt') as b,gzip.open(base/'FFS.sumstats.gz','rt') as c,gzip.open(base/'insomnia.sumstats.gz','rt') as d:
 n=0
 for rows in zip(*(csv.DictReader(f,delimiter='\t') for f in (a,b,c,d)),strict=True):
  assert all(key(x)==key(rows[0]) for x in rows);n+=1
assert n==count
r={'source_sha256':sha(src),'hfrs_input_sha256':sha(out),'matching_rows':count,'source_rows_skipped_before_last_match':skipped,'ordered_four_way_keys_identical':True,'earlier_input_sha256':receipt['output_sha256'],'warning':'Exact keys and hash do not validate original phenotype, per-variant N, ancestry or overlap.'}
(p/'matching_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
