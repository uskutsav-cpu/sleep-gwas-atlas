import csv,gzip,hashlib,io,json,math,urllib.request
from pathlib import Path
url='https://www.dropbox.com/s/foxe08tsyxi7bzq/884.ldsc.imputed_v3.both_sexes.tsv.bgz?dl=1'
with urllib.request.urlopen(url,timeout=120) as r: raw=r.read()
assert hashlib.sha256(raw).hexdigest()=='fffebe514901e7c09b242365f3a327d0a73d4df6729d64fb903f45a4eb7ad5be'
activity={}; nset=set(); invalid=duplicates=0
for d in csv.DictReader(io.StringIO(gzip.decompress(raw).decode()),delimiter='\t'):
 try:
  z=float(d['Z']); n=float(d['N']); assert math.isfinite(z) and math.isfinite(n) and n>0 and d['A1'] in 'ACGT' and d['A2'] in 'ACGT' and d['A1']!=d['A2']
 except (ValueError,AssertionError): invalid+=1; continue
 nset.add(n)
 if d['SNP'] in activity: duplicates+=1
 activity[d['SNP']]=frozenset((d['A1'],d['A2']))
complement=str.maketrans('ACGT','TGCA'); result={}
for trait in ('questionnaire','registry'):
 matched=shared=palindrome=0
 with gzip.open('/Volumes/Brain6Work/physical-function-source-sensitivity-20261028/'+trait+'.sumstats.gz','rt') as f:
  for d in csv.DictReader(f,delimiter='\t'):
   if d['SNP'] not in activity: continue
   shared+=1; pair=frozenset((d['A1'],d['A2'])); a=activity[d['SNP']]
   if a in (frozenset('AT'),frozenset('CG')): palindrome+=1
   if a==pair or a==frozenset(x.translate(complement) for x in pair): matched+=1
 result[trait]={'shared':shared,'allele_pair_or_strand_match':matched,'palindromic':palindrome}
result.update({'source_unique_ids':len(activity),'duplicates':duplicates,'invalid':invalid,'N_values':sorted(nset)})
Path(__file__).with_name('audit_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
