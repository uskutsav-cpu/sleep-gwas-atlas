#!/usr/bin/env python3
"""Count source/harmonized allele fields; not a scoring or phenotype analysis."""
import csv,gzip,hashlib,json
from collections import Counter
from pathlib import Path
root=Path('/Volumes/Extreme SSD/codex-sleep-frailty-work')
out={}
for score in ('PGS005228','PGS005229'):
 path=root/f'{score}_hm38.txt.gz'; c=Counter();h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 with gzip.open(path,'rt',newline='') as f:
  reader=csv.DictReader((line for line in f if not line.startswith('#')),delimiter='\t')
  for row in reader:
   c['rows']+=1
   a=row['effect_allele'].upper(); b=row['other_allele'].upper(); v=row['hm_inferOtherAllele']
   if a==b:c['effect_equals_other']+=1
   if not b:c['other_missing']+=1
   if v:c['inferred_other_present']+=1
   if not b:
    weight=float(row['effect_weight'])
    c['missing_other_sum_abs_weights']+=abs(weight)
    c['missing_other_sum_squared_weights']+=weight*weight
    if not v:c['missing_other_no_inference']+=1
    elif '/' in v:c['missing_other_multiple_inferred']+=1
    else:
     c['missing_other_single_inferred']+=1
     if {a,v.upper()} in ({'A','T'},{'C','G'}):c['missing_other_single_inferred_palindromic']+=1
   else:
    weight=float(row['effect_weight'])
    c['reported_other_sum_abs_weights']+=abs(weight)
    c['reported_other_sum_squared_weights']+=weight*weight
   if not row['hm_chr'] or not row['hm_pos']:c['harmonized_position_missing']+=1
 out[score]={'sha256':h.hexdigest(),'bytes':path.stat().st_size,'counts':dict(c)}
Path('outputs/stage6_foote_hm38_allele_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
