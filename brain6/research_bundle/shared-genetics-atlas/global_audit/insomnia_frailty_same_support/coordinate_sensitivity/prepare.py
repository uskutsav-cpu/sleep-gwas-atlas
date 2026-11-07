import csv,gzip,hashlib,json,pathlib
from contextlib import ExitStack
p=pathlib.Path(__file__).resolve().parent;prev=p.parent;src=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity');out=pathlib.Path('work/coordinate_ldsc_inputs');out.mkdir(parents=True,exist_ok=True)
with (prev/'hfrs_coordinate_unresolved.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
exclude={x['rsid'] for x in r};assert len(r)==len(exclude)==537
fields=['SNP','A1','A2','N','Z'];files=[src/f'{x}.sumstats.gz' for x in ('FI','FFS','HFRS','S')];dest=[out/f'{x}.sumstats.gz' for x in ('FI','FFS','HFRS','S')];seen=set();n=0
with ExitStack() as stack:
 readers=[csv.DictReader(stack.enter_context(gzip.open(f,'rt')),delimiter='\t') for f in files];writers=[]
 for f in dest:
  w=csv.DictWriter(stack.enter_context(gzip.open(f,'wt',newline='')),fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();writers.append(w)
 for rows in zip(*readers,strict=True):
  assert all((x['SNP'],x['A1'],x['A2'])==(rows[0]['SNP'],rows[0]['A1'],rows[0]['A2']) for x in rows)
  rs=rows[0]['SNP']
  if rs in exclude:assert rs not in seen;seen.add(rs);continue
  for w,row in zip(writers,rows,strict=True):w.writerow(row)
  n+=1
assert seen==exclude and n==1134416
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1<<20),b''):h.update(b)
 return h.hexdigest()
r={'input_rows':1134953,'excluded':len(seen),'remaining':n,'excluded_list_sha256':sha(prev/'hfrs_coordinate_unresolved.tsv'),'inputs':{str(x):sha(x) for x in files},'outputs':{str(x):sha(x) for x in dest},'limitations':'Post-hoc support sensitivity; exclusion does not prove markers are invalid.'};(p/'support_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print('selected support',n,'excluded',len(seen))
