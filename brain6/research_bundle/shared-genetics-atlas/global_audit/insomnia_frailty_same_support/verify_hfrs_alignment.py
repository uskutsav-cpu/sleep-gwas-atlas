import csv,gzip,json,hashlib,math,pathlib,collections
p=pathlib.Path(__file__).resolve().parent;base=pathlib.Path('/Volumes/Extreme SSD/Codex Archive/hfrs_fi_ffs_common_support');parent=base/'HFRS.sumstats.gz';upstream=pathlib.Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/munged/HFRS.sumstats.gz');receipt=json.loads((base/'support.json').read_text());comp=str.maketrans('ACGT','TGCA')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as z:
  for b in iter(lambda:z.read(1<<20),b''):h.update(b)
 return h.hexdigest()
assert sha(parent)==receipt['output_sha256']['HFRS'] and sha(upstream)==receipt['input_sha256']['HFRS']
source={};duplicates=set()
with gzip.open(upstream,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  rs=r['SNP']
  try:
   z=float(r['Z']); n=float(r['N']); assert math.isfinite(z) and math.isfinite(n) and n>0
  except (ValueError, AssertionError):continue
  if rs in source:duplicates.add(rs)
  else:source[rs]=(r['A1'],r['A2'],z,n)
assert not duplicates
counts=collections.Counter();problems=[]
with gzip.open(parent,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  a,b=r['A1'],r['A2'];s=source.get(r['SNP']);counts['parent_rows']+=1
  if s is None:counts['missing']+=1;problems.append((r['SNP'],'missing'));continue
  x,y,z,n=s
  if any(len(c)!=1 or c not in 'ACGT' for c in (a,b,x,y)) or a==b or ''.join(sorted(a+b)) in ('AT','CG'):
   counts['invalid_or_palindromic']+=1;problems.append((r['SNP'],'invalid_or_palindromic'));continue
  if (x,y)==(a,b):mode='direct';expected=z
  elif (x.translate(comp),y.translate(comp))==(a,b):mode='complement';expected=z
  elif (y,x)==(a,b):mode='swapped';expected=-z
  elif (y.translate(comp),x.translate(comp))==(a,b):mode='complement_swapped';expected=-z
  else:counts['allele_mismatch']+=1;problems.append((r['SNP'],'allele_mismatch'));continue
  if float(r['Z'])!=expected or float(r['N'])!=n:
   counts['value_mismatch']+=1;problems.append((r['SNP'],'value_mismatch'));continue
  counts[mode]+=1
assert counts['parent_rows']==1137927 and not problems,(counts,problems[:10]);assert sum(counts[m] for m in ('direct','complement','swapped','complement_swapped'))==counts['parent_rows']
r={'HFRS_munged_upstream_path':str(upstream),'HFRS_upstream_sha256':sha(upstream),'aligned_parent_path':str(parent),'aligned_parent_sha256':sha(parent),'counts':dict(counts),'duplicates_in_upstream':len(duplicates),'status':'PASS','limitation':'Verifies transformation of already-munged Z to FI-oriented alleles, not correctness of original raw effect allele, GWAS, participant overlap or per-variant N.'};(p/'hfrs_alignment_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r['counts'],indent=2))
