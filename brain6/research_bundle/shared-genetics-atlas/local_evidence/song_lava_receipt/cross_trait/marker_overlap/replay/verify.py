#!/usr/bin/env python3
import hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent
a=json.loads((p/'insomnia.json').read_text());b=json.loads((p/'sleep_apnea.json').read_text())
assert a['n_snps']==2741 and a['K']==216 and a['failed']==[True,False]
assert b['n_snps']==2452 and b['K']==185 and b['failed']==[False,False]
sa,sb=set(a['snps']),set(b['snps']);assert len(sa)==2741 and len(sb)==2452 and sb<=sa
expected=set((p.parent/'insomnia_only_ids.tsv').read_text().splitlines()[1:]);assert sa-sb==expected
for d in (a,b):
 for i in range(2):
  calculated=(d['delta_squared_sums'][i]/d['K']-d['sigma'][i][i])*d['nref_scale']
  assert math.isclose(calculated,d['omega'][i][i],rel_tol=1e-10,abs_tol=1e-17)
assert a['omega'][0][0]<0 and b['omega'][0][0]>0
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
files=[p/'replay.R',p/'insomnia.json',p/'sleep_apnea.json',p/'insomnia_quoted_input.tsv',p/'sleep_apnea_quoted_input.tsv',
 Path('/Volumes/Extreme SSD/brain6-work/lava-ukb-v1.1/lava-ukb-v1.1_chr3.bcor'),
 Path('/Volumes/Extreme SSD/brain6-work/r-env-4.3.3-short/lib/R/library/LAVA/R/LAVA.rdb')]
receipt={'files_sha256':{str(x):sha(x) for x in files},'exact_nested_internal_snps':True,
 'algebra_checked':True,'insomnia_fi_omega':a['omega'][0][0],'apnea_fi_omega':b['omega'][0][0],
 'insomnia_uncapped_fi_h2':a['omega'][0][0]*a['K'],'apnea_uncapped_fi_h2':b['omega'][0][0]*b['K']}
(p/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
