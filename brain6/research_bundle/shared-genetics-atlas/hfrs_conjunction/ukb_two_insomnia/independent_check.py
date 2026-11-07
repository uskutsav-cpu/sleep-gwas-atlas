#!/usr/bin/env python3
"""Recompute unguarded upper-sensitivity p-values, no source-code imports."""
import gzip, csv, math, json
from itertools import zip_longest
from pathlib import Path
import numpy as np

base=Path('/Volumes/Brain6Work')
paths=[base/'ukb_hfrs_public_20261028.sumstats.gz',
       base/'ukb-hfrs-comparison-20261028/questionnaire.sumstats.gz',
       base/'ukb-hfrs-comparison-20261028/registry.sumstats.gz']
streams=[gzip.open(p,'rt') for p in paths]
readers=[csv.DictReader(f,delimiter='\t') for f in streams]
with Path(__file__).with_name('top100_per_source.tsv').open() as f:
    top={(x['source'],x['SNP']):x for x in csv.DictReader(f,delimiter='\t')}
assert len(top)==200
checked=set()
values=[]; minima=[(1,None),(1,None)]; n=0
try:
    for h,q,r in zip_longest(*readers):
        assert h and q and r and h['SNP']==q['SNP']==r['SNP']
        ph=math.erfc(abs(float(h['Z']))/math.sqrt(2))
        for j,x in enumerate((q,r)):
            p=max(ph,math.erfc(abs(float(x['Z']))/math.sqrt(2)))
            values.append(p)
            if p < minima[j][0]: minima[j]=(p,h['SNP'])
            key=(('questionnaire','registry')[j],h['SNP'])
            if key in top:
                t=top[key]
                guarded=max(math.erfc(max(0,abs(float(h['Z']))-1e-9)/math.sqrt(2)),
                            math.erfc(max(0,abs(float(x['Z']))-.0005)/math.sqrt(2)))
                assert abs(float(t['max_p'])-guarded)<1e-15
                assert float(t['insomnia_Z_on_hfrs_A1']) == (float(x['Z']) if h['A1']==x['A1'] else -float(x['Z']))
                checked.add(key)
        n+=1
finally:
    for f in streams: f.close()
assert n==1020621 and checked==set(top)
p=np.sort(np.array(values)); m=len(p); thresholds=.05*np.arange(1,m+1)/m
bh=int(np.max(np.flatnonzero(p<=thresholds))+1) if np.any(p<=thresholds) else 0
bon=int(np.count_nonzero(p<=.05/m))
summary={'rows':n,'slots':m,'ranked_rows_rechecked':len(checked),'bonferroni_unguarded':bon,'bh_unguarded':bh,
         'minima_unguarded':minima}
assert bon==0 and bh==0 and minima[0][1]=='rs1476535' and minima[1][1]=='rs7714971'
Path(__file__).with_name('independent_check.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
