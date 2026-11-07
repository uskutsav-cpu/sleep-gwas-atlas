import csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent; atlas=next(q for q in p.parents if q.name=='shared-genetics-atlas'); m=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged'); f=Path('/Volumes/Extreme SSD/Codex Archive/hfrs_ukb_sensitivity')
diseases=['adhd','mdd','scz','bipolar','alz','parkinson']; targets=['FI','FFS','HFRS']; files={x:m/(x+'.sumstats.gz') for x in diseases}|{x:f/(x+'.sumstats.gz') for x in targets}
prev=next((p.parent).glob('input-hashes-*.txt')); old={Path(line.split('  ',1)[1]).name:line.split()[0] for line in prev.read_text().splitlines()};support=json.loads((f/'support.json').read_text())['output_sha256']
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as g:
  for b in iter(lambda:g.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
for k,path in files.items():
 actual=sha(path);expected=old[path.name] if k in diseases else support[k];assert actual==expected, f'{k} hash mismatch'
comp=str.maketrans('ACGT','TGCA');stats={};common=None
for name in targets+diseases:
 path=files[name]; seen=set(); dup=set();valid={}; total=bad=0
 with gzip.open(path,'rt',newline='') as h:
  for row in csv.DictReader(h,delimiter='\t'):
   total+=1;rs=row['SNP']; a=row['A1'];b=row['A2']
   if rs in seen:dup.add(rs);continue
   seen.add(rs)
   try:
    z=float(row['Z']);n=float(row['N']);ok=math.isfinite(z) and math.isfinite(n) and n>0
   except (ValueError,TypeError):ok=False
   if not ok or a not in 'ACGT' or b not in 'ACGT' or len(a)!=1 or len(b)!=1 or a==b or frozenset((a,b)) in (frozenset('AT'),frozenset('CG')):bad+=1;continue
   pair=''.join(sorted((a,b))); pair2=''.join(sorted((a.translate(comp),b.translate(comp))));key=min(pair,pair2)
   if common is None or rs in common:valid[rs]=key
 if common is None:common=valid
 else:common={rs:key for rs,key in common.items() if rs in valid and valid[rs]==key and rs not in dup}
 if dup:
  for rs in dup:common.pop(rs,None)
 stats[name]={'source_rows':total,'valid_candidates_in_current_support':len(valid),'duplicate_ids':len(dup),'invalid_or_palindromic_rows':bad,'common_remaining':len(common),'sha256':sha(path)}
 print(name,len(common),flush=True)
assert len(common)>0
(p/'receipt.json').write_text(json.dumps({'input_checks':stats,'n_nine_input_intersection':len(common),'mvp_overlap':'not assessed','ld_reference_overlap':'not assessed','alleles':'unordered, allowing strand complement; no palindromes','analysis':'feasibility only'},indent=2)+'\n')
