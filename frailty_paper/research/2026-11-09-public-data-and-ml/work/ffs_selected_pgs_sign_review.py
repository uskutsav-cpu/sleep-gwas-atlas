#!/usr/bin/env python3
"""Compare marginal and published PGS signs only after exact matching."""
import csv,hashlib,json,math
from pathlib import Path
src=Path('outputs/stage6_ffs_selected_pgs_compare.csv')
rows=list(csv.DictReader(src.open(newline='')))
by_variant={}
for r in rows:
 if r['source_exact_pair_rows']!='1' or r['pgs_coord_rows']!='1' or r['pgs_effect_equals_source_A1']!='True' or r['pgs_other_valid_pair']!='True': continue
 a,b=float(r['source_beta']),float(r['pgs_weight'])
 assert math.isfinite(a) and math.isfinite(b)
 val={'variant':r['variant'],'source_beta_A1':a,'published_pgs_weight_A1':b,'sign_relation':'same' if a*b>0 else 'opposite' if a*b<0 else 'zero','score_lists':[]}
 if r['variant'] not in by_variant: by_variant[r['variant']]=val
 else:
  prev=by_variant[r['variant']]
  assert (prev['source_beta_A1'],prev['published_pgs_weight_A1'])==(a,b)
 by_variant[r['variant']]['score_lists'].append(r['list'])
out=Path('outputs/ffs_selected_pgs_sign_review.csv')
with out.open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['variant','source_beta_A1','published_pgs_weight_A1','sign_relation','score_lists']);
 for v in sorted(by_variant.values(),key=lambda v:v['variant']):w.writerow([v['variant'],v['source_beta_A1'],v['published_pgs_weight_A1'],v['sign_relation'],';'.join(v['score_lists'])])
counts={k:sum(v['sign_relation']==k for v in by_variant.values()) for k in ('same','opposite','zero')}
result={'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'distinct_variants':len(by_variant),'list_rows':sum(len(x['score_lists']) for x in by_variant.values()),'sign_counts':counts,'note':'Same GWAS; conditional published score coefficients are not independent marginal estimates; sign agreement does not validate the project score.'}
Path('outputs/ffs_selected_pgs_sign_review.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
