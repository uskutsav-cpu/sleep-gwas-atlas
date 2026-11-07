import collections,csv,gzip,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent; parent=p.parent
m=Path('/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data')
mun=m/'munged/insomnia.sumstats.gz';raw=m/'raw/insomnia.txt.gz'
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as x:
  for b in iter(lambda:x.read(1048576),b''):h.update(b)
 return h.hexdigest()
assert sha(mun)==json.loads((parent/'receipt.json').read_text())['insomnia']['sha256']
assert sha(raw)==json.loads((parent/'raw_sample.json').read_text())['insomnia']['raw_sha256']
rows={};bad_m=collections.Counter()
with gzip.open(mun,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  try:
   if not all(math.isfinite(float(row[k])) for k in ('Z','N')) or float(row['N'])<=0:raise ValueError()
  except (TypeError,ValueError):bad_m['missing_or_invalid_z_n']+=1;continue
  a,b=row['A1'],row['A2']
  if a not in 'ACGT' or b not in 'ACGT' or a==b or {a,b} in ({'A','T'},{'C','G'}):bad_m['invalid_or_ambiguous_alleles']+=1;continue
  snp=row['SNP'];assert snp not in rows;rows[snp]=row
rawhits=collections.defaultdict(list);raw_count=0
with gzip.open(raw,'rt') as f:
 for row in csv.DictReader(f,delimiter='\t'):
  raw_count+=1;snp=row['SNP']
  if snp in rows:rawhits[snp].append((row['A1'],row['A2'],row['N']))
base=p/'matched_constantN.sumstats.gz';alt=p/'matched_rawN.sumstats.gz';assert not base.exists() and not alt.exists()
stats=collections.Counter();ratios=[]
with gzip.open(base,'wt') as b,gzip.open(alt,'wt') as a:
 header='SNP\tA1\tA2\tZ\tN\n';a.write(header);b.write(header)
 for snp,row in rows.items():
  hits=rawhits[snp]
  if len(hits)!=1:stats['raw_missing' if not hits else 'raw_duplicate_id']+=1;continue
  ra,rb,n=hits[0]
  if {ra,rb}!={row['A1'],row['A2']}:stats['allele_mismatch']+=1;continue
  try:
   v=float(n)
   if not math.isfinite(v) or v<=0:raise ValueError()
  except (TypeError,ValueError):stats['raw_invalid_n']+=1;continue
  prefix='\t'.join(row[k] for k in ('SNP','A1','A2','Z'))+'\t';b.write(prefix+row['N']+'\n');a.write(prefix+n+'\n');stats['matched']+=1;ratios.append(v/float(row['N']))
with gzip.open(base,'rt') as b,gzip.open(alt,'rt') as a:
 assert next(b)==next(a);count=0
 for x,y in zip(b,a):
  X=x.rstrip().split('\t');Y=y.rstrip().split('\t');assert X[:4]==Y[:4] and float(X[4])==313750.036;count+=1
 assert count==stats['matched'] and not a.read() and not b.read()
ratios.sort();summary={'raw_rows':raw_count,'eligible_munged_rows':len(rows),'munged_exclusions':dict(bad_m),'matching':dict(stats),'ratio_rawN_to_constant':{'min':ratios[0],'median':ratios[len(ratios)//2],'max':ratios[-1]},'source_hashes':{'raw':sha(raw),'munged':sha(mun)},'derived_hashes':{'baseline':sha(base),'alternate':sha(alt)},'paired_rows_verified':count}
(p/'build_receipt.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
