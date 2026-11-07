import csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent; atlas=p.parents[3]; src=Path('/Volumes/Brain6Work/brain6-frailty-paired'); alt=atlas/'global_audit/input_n_audit/insomnia_n_sensitivity/matched_rawN.sumstats.gz'; old=json.loads((p.parent/'matching_receipt.json').read_text()); ref=json.loads((alt.parent/'build_receipt.json').read_text())
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as fd:
  for b in iter(lambda:fd.read(1048576),b''):h.update(b)
 return h.hexdigest()
expected={'insomnia':old['earlier_input_sha256'][str(src/'insomnia.sumstats.gz')],'FI':old['earlier_input_sha256'][str(src/'FI.sumstats.gz')],'FFS':old['earlier_input_sha256'][str(src/'FFS.sumstats.gz')],'HFRS':old['hfrs_input_sha256']}
for k,v in expected.items():assert sha(src/(k+'.sumstats.gz'))==v
assert sha(alt)==ref['derived_hashes']['alternate']
lookup={}
with gzip.open(alt,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  assert r['SNP'] not in lookup;lookup[r['SNP']]=(r['A1'],r['A2'],r['Z'],r['N'])
out=p/'insomnia.sumstats.gz';assert not out.exists();count=0
with gzip.open(src/'insomnia.sumstats.gz','rt') as f,gzip.open(out,'wt') as o:
 header=next(f);assert header.strip()=='SNP\tA1\tA2\tZ\tN';o.write(header);seen=set()
 for row in csv.DictReader(f,fieldnames=['SNP','A1','A2','Z','N'],delimiter='\t'):
  snp=row['SNP'];assert snp not in seen and snp in lookup;seen.add(snp);a,b,z,n=lookup[snp]
  assert {a,b}=={row['A1'],row['A2']} and {a,b} not in ({'A','T'},{'C','G'})
  assert abs(abs(float(z))-abs(float(row['Z'])))<1e-10 and math.isfinite(float(n)) and float(n)>0
  o.write('\t'.join([snp,row['A1'],row['A2'],row['Z'],n])+'\n');count+=1
assert count==old['matching_rows'];result={'rows':count,'source_files_sha256':expected,'raw_n_mapping_sha256':sha(alt),'alternate_insomnia_sha256':sha(out),'zero_missing_or_allele_mismatches':True,'note':'Source N semantics unverified; synthetic sensitivity only.'};(p/'build_receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
