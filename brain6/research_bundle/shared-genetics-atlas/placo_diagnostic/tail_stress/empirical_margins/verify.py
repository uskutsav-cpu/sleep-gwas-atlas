#!/usr/bin/env python3
import collections,csv,gzip,hashlib,json,math,statistics
from pathlib import Path
p=Path(__file__).resolve().parent;base=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/common_support')
with (p/'marginal_tails.tsv').open() as f:out=list(csv.DictReader(f,delimiter='\t'))
cut=statistics.NormalDist().inv_cdf(1-.0001/2);counts=collections.Counter();ns=collections.Counter()
for pair in ('FI','FFS'):
 with gzip.open(base/'S.sumstats.gz','rt') as f,gzip.open(base/(pair+'.sumstats.gz'),'rt') as g:
  for a,b in zip(csv.DictReader(f,delimiter='\t'),csv.DictReader(g,delimiter='\t'),strict=True):
   assert (a['SNP'],a['A1'],a['A2'])==(b['SNP'],b['A1'],b['A2'])
   x=float(a['Z']);y=float(b['Z'])
   if max(x*x,y*y)>80 or (abs(x)>cut and abs(y)>cut):continue
   ns[pair]+=1
   for trait,z in (('S',x),(pair,y)):
    for t in (3,4,5,6):counts[(pair,trait,t)]+=int(abs(z)>t)
for r in out:
 key=(r['pair'],r['trait'],int(r['abs_z_threshold']));assert int(r['n'])==ns[r['pair']] and int(r['observed_count'])==counts[key]
 assert math.isclose(float(r['observed_div_expected']),int(r['observed_count'])/float(r['normal_reference_expected_count']),rel_tol=1e-12)
receipt={'threshold_rows_verified':len(out),'independent_full_rescan':True,'subset_counts':dict(ns),'cutoff':cut,'sha256':hashlib.sha256((p/'marginal_tails.tsv').read_bytes()).hexdigest(),'limits':'Confirms counts and arithmetic, not normality of null effects or calibration.'}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
