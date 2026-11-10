from pathlib import Path
import json,pandas as pd
O=Path('outputs');s=pd.read_csv(O/'stage6_grip_sparse_score_summary.csv');raw=pd.read_csv(O/'stage6_grip_gwas_raw_candidates.csv');ids=pd.read_csv(O/'unvalidated_15_variant_score_prototype.csv').variant
v=raw[raw.variant_id.isin(ids)].copy();assert len(v)==15 and v.variant_id.is_unique
v['question_marks']=v.Direction.str.count(r'\?');assert v.Direction.str.len().eq(22).all()
a={}
for w,g in s.groupby('weight'):
 full=g.loc[g.omitted_variant.eq('none'),'z_approx'].item();d=g.loc[g.omitted_variant.ne('none')].copy();d['abs_change']=(d.z_approx-full).abs();q=d.loc[d.abs_change.idxmax()];a[w]={'full_z':full,'largest_change_variant':q.omitted_variant,'z_omitting_variant':float(q.z_approx)}
out={'direction_question_mark_counts':{str(k):int(n) for k,n in v.question_marks.value_counts().sort_index().items()},'leave_one_out':a,'interpretation':'Question marks suggest missing contributing study results under METAL conventions; source-specific per-SNP sample sizes and null covariance cannot be recovered from these strings.'}
(O/'stage6_grip_influence_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
