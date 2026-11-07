#!/usr/bin/env python3
"""Post-hoc descriptive association of audit N and residual flags; no inference."""
import collections,csv,json,statistics
from pathlib import Path
p=Path(__file__).resolve().parent;old=p.parents[2]/'rs62264764_dense_qtl/multisignal_pilot'
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
n={x['variant']:float(x['adhd_N_AUDIT']) for x in read(old/'sample_sizes.tsv')}
rows=[x for x in read(p/'repair_ldz_details.tsv') if x['trait']=='adhd']
out={}
for ref in ('eigenclip','ridge_0p001','1000g_eur'):
 r=[x for x in rows if x['reference']==ref];assert len(r)==1547
 c=collections.Counter(('below_80000' if n[x['variant']]<80000 else 'at_least_80000', 'residual_gt5' if abs(float(x['z_std_diff']))>5 else 'residual_le5') for x in r)
 out[ref]={'markers':len(r),'median_n_audit':statistics.median(n[x['variant']] for x in r),'strata':{a:{b:c[(a,b)] for b in ('residual_gt5','residual_le5')} for a in ('below_80000','at_least_80000')}}
result={'threshold_selection':'80,000 selected post hoc after inspecting top outlier effective Ns','residual_cutoff':'absolute standardized residual >5, descriptive and post hoc','strata':out,'limitations':'Counts are correlated variants, not independent observations. N_AUDIT is computed from reported case/control counts and does not capture pairwise sample overlap. The grouping cannot establish cause or a calibrated association.'}
(p/'adhd_coverage_strata.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
