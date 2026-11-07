#!/usr/bin/env python3
"""Compare reference-compatible locus IDs in saved FI and sleep input files."""
import csv,gzip,hashlib,json
from pathlib import Path

base=Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/lava_sensitivity_v1/inputs')
ref=Path('/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3.info')
out=Path(__file__).resolve().parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(8*1024*1024),b''):h.update(c)
 return h.hexdigest()
refs={}
with ref.open() as f:
 for row in csv.DictReader(f,delimiter='\t'):
  if 47588462 <= int(row['POS']) <= 50387742:
   assert row['SNP'] not in refs
   refs[row['SNP']]=(row['A1'],row['A2'])
sets={};sources={}
for trait in ('frailty','insomnia','sleep_apnea'):
 p=base/trait/'chr03.sumstats.tsv.gz';matched=set();bad=0;duplicates=0
 with gzip.open(p,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   snp=row['SNP']
   if snp not in refs:continue
   if set((row['A1'],row['A2']))!=set(refs[snp]):bad+=1;continue
   if snp in matched:duplicates+=1
   matched.add(snp)
 sets[trait]=matched;sources[trait]={'path':str(p),'sha256':sha(p),'matched':len(matched),'mismatch_rows':bad,'duplicate_rows':duplicates}
 assert bad==duplicates==0
a=sets['frailty'] & sets['insomnia'];b=sets['frailty'] & sets['sleep_apnea']
assert len(refs)==3065 and len(a)==2741 and len(b)==2452 and b<=a
with (out/'insomnia_only_ids.tsv').open('w') as f:
 f.write('SNP\n');f.writelines(s+'\n' for s in sorted(a-b))
result={'reference':str(ref),'reference_sha256':sha(ref),'reference_ids_in_locus':len(refs),
 'source_files':sources,'insomnia_fi_count':len(a),'apnea_fi_count':len(b),
 'shared_count':len(a&b),'insomnia_only_count':len(a-b),'apnea_only_count':len(b-a),
 'insomnia_only_ids_sha256':sha(out/'insomnia_only_ids.tsv')}
(out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
