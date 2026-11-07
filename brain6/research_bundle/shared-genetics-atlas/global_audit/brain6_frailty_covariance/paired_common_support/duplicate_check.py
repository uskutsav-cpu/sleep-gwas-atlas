import csv,gzip,json
from pathlib import Path
p='/Volumes/Extreme SSD/Utsav-Research-Archive/Sleep-GWAS/Sep1-local-dependencies/data/munged/insomnia.sumstats.gz'
seen=set();duplicates=set()
with gzip.open(p,'rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  snp=r['SNP']
  if snp in seen:duplicates.add(snp)
  seen.add(snp)
assert not duplicates
Path(__file__).with_name('duplicate_check.json').write_text(json.dumps({'unique_ids':len(seen),'duplicated_ids':len(duplicates),'includes_invalid_and_missing_summary_rows':True},indent=2)+'\n')
