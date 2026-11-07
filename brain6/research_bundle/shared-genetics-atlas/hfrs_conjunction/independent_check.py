#!/usr/bin/env python3
"""Independent matrix-free cross-check of eligible count, minimum and BH/BY."""
import csv,json,math
import numpy as np
from pathlib import Path
base=Path('/Volumes/Extreme SSD/sleep-frailty-full-autonomous-v5-20261009/frailty_paper/full_autonomous_research_v5/genetics/new_public_validation/filtered')
out=Path('outputs/shared-genetics-atlas/hfrs_conjunction')
with (base/'R13_INSOMNIA_HM3.tsv').open() as f:
 a={}
 for r in csv.DictReader(f,delimiter='\t'):
  assert r['SNP'] not in a
  a[r['SNP']]=(r['CHR'],r['A1'],r['A2'],float(r['P']))
m=len(a); seen=set(); p=[]; strongest=(1,None)
with (base/'HFRS_HM3.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='\t'):
  name=r['SNP']
  if name not in a:continue
  assert name not in seen;seen.add(name)
  ch,x,y,q=a[name]; u,v=r['A1'],r['A2'];t=float(r['P'])
  if ch!=r['CHR'] or len({x,y})!=2 or any(len(z)!=1 or z not in 'ACGT' for z in (x,y,u,v)) or {x,y}!={u,v}:continue
  if not (math.isfinite(q) and math.isfinite(t) and 0<q<=1 and 0<t<=1):continue
  val=max(q,t);p.append(val)
  if val<strongest[0]:strongest=(val,name)
p=np.sort(np.array(p)); ranks=np.arange(1,len(p)+1); bh=float(min(1,np.min(m*p/ranks))); harmonic=sum(1/i for i in range(1,m+1)); by=min(1,bh*harmonic)
s=json.loads((out/'summary.json').read_text());assert len(p)==s['status']['eligible'] and m==s['family_size_insomnia_rows'];assert strongest[0]==s['min_conjunction_p'];assert (bh<=.05)==(s['BH_rejections_padded_family']>0);assert (by<=.05)==(s['BY_rejections_padded_family']>0)
result={'n':m,'eligible':len(p),'min_max_p':strongest[0],'min_variant':strongest[1],'minimum_BH_adjusted':bh,'minimum_BY_adjusted':by,'agreement_with_screen':True,'note':'Independent code path; assertion checks eligibility/count and minima, not upstream source validity or full allele alignment.'}
(out/'independent_check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
