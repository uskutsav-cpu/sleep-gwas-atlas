import csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent; receipt=json.loads((p/'prepare_receipt.json').read_text()); d=Path(receipt['output_directory'])
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
common=d/'walking.common.sumstats.gz'; source=d/'questionnaire.bgz'; out=d/'questionnaire.common.sumstats.gz'
assert sha(common)==receipt['output_sha256']['walking']
src=gzip.open(source,'rt'); ref=gzip.open(common,'rt'); rd=csv.DictReader(src,delimiter='\t'); targets=csv.DictReader(ref,delimiter='\t')
assert set(rd.fieldnames)==set(targets.fieldnames)=={'SNP','A1','A2','N','Z'}
target=next(targets,None); last=''; n=selected=0; complement=str.maketrans('ACGT','TGCA')
with gzip.open(out,'xt') as f:
 w=csv.DictWriter(f,fieldnames=['SNP','A1','A2','N','Z'],delimiter='\t',lineterminator='\n');w.writeheader()
 for row in rd:
  n+=1
  assert row['SNP']>last;last=row['SNP']
  assert math.isfinite(float(row['Z'])) and 0<float(row['N'])<1e9
  while target is not None and target['SNP']<row['SNP']: raise ValueError(('missing required marker',target['SNP']))
  if target is None or row['SNP']!=target['SNP']:continue
  a,b=target['A1'],target['A2'];c,e=row['A1'],row['A2']
  assert {c,e}=={a,b} or {c.translate(complement),e.translate(complement)}=={a,b}
  w.writerow({k:row[k] for k in w.fieldnames});selected+=1;target=next(targets,None)
assert target is None and selected==receipt['common_nonpalindromic_rows']
src.close();ref.close()
outdata={'url':'https://www.dropbox.com/s/zzvfqqfq9cmym9s/1200.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1','source_sha256':sha(source),'source_rows':n,'selected_rows':selected,'output_sha256':sha(out),'common_reference_sha256':sha(common),'output_path':str(out)}
(p/'questionnaire_receipt.json').write_text(json.dumps(outdata,indent=2)+'\n');print(json.dumps(outdata,indent=2))
