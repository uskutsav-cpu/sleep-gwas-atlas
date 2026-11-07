#!/usr/bin/env python3
import csv,gzip,collections,hashlib,json,pathlib
out=pathlib.Path('outputs/shared-genetics-atlas/placo_diagnostic/completed_audit_fi_only'); fi=pathlib.Path('/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1/analysis-workspace/frailty_v1/harmonized/frailty.harmonized.tsv.gz')
with (out/'candidate_statistic_sensitivity.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
ids={x['SNP'] for x in r};found=collections.defaultdict(list)
with gzip.open(fi,'rt') as f:
 for a in csv.DictReader(f,delimiter='\t'):
  if a['SNP'] in ids:found[a['SNP']].append(a)
for x in r:
 a=found[x['SNP']];ok=[y for y in a if {y['A1'],y['A2']}=={x['A1'],x['A2']}]
 if len(a)==len(ok)==1:x.update(CHR=ok[0]['CHR'],BP_GRCH37=ok[0]['BP'],position_status='UNIQUE_RSID_ALLELE_MATCH')
 else:x.update(CHR='',BP_GRCH37='',position_status='MISSING_OR_AMBIGUOUS')
 x['max_p_source']=max(float(x['sleep_p']),float(x['frailty_p']))
with (out/'candidate_positions.tsv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(r[0]),delimiter='\t');w.writeheader();w.writerows(r)
print('positions',collections.Counter(x['position_status'] for x in r),'chromosomes',collections.Counter(x['CHR'] for x in r))
for x in sorted(r,key=lambda x:(int(x['CHR']) if x['CHR'] else 99,int(x['BP_GRCH37']) if x['BP_GRCH37'] else 0)):
 print(x['SNP'],x['CHR'],x['BP_GRCH37'],x['p_placo_plus'],x['max_p_source'])
