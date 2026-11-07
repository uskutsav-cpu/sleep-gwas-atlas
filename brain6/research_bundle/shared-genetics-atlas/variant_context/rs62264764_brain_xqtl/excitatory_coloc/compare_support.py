import csv,hashlib,json,math
from pathlib import Path
p=Path(__file__).resolve().parent;root=Path('/Volumes/Brain6Work/xqtl-chr3');result={}
for t in ('insomnia','adhd'):
 def load(path):
  with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
 a=load(root/'exc_matched'/('exc_'+t+'.tsv'));b=load(root/'matched'/('oli_'+t+'.tsv'))
 key=lambda x:(x['variant'],x['gwas_rsid']);aa={key(x):x for x in a};bb={key(x):x for x in b};assert len(aa)==len(a) and len(bb)==len(b) and aa.keys()==bb.keys()
 for k in aa:
  for c in ('gwas_beta_alt','gwas_se'):assert math.isclose(float(aa[k][c]),float(bb[k][c]),rel_tol=1e-12,abs_tol=1e-12)
 result[t]={'same_variant_keys_and_gwas_effects':True,'n':len(aa),'note':'Different molecular datasets share donors; matching support does not prove shared causal variants.'}
(p/'support_comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
