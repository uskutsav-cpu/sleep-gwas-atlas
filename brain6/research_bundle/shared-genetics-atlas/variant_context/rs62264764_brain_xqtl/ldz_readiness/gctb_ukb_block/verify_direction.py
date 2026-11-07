#!/usr/bin/env python3
import collections,csv,hashlib,json
from pathlib import Path
p=Path(__file__).resolve().parent
source=json.loads((p/'adhd_direction_audit.json').read_text())
with (p/'repair_ldz_details.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert source['matched_variants']==1547 and source['direction_length_counts']=={'13':1547}
assert sum(source['direction_character_counts'].values())==1547*13
out={}
for ref,groups in source['reference_groups'].items():
 expected=collections.Counter('gt5' if abs(float(x['z_std_diff']))>5 else 'le5' for x in rows if x['trait']=='adhd' and x['reference']==ref)
 actual=collections.Counter()
 for g in groups:
  assert len(g['question_mark_mask'])==13 and set(g['question_mark_mask'])<=set('.?')
  actual[g['residual_group']]+=g['count']
 assert expected==actual and sum(actual.values())==1547
 out[ref]={'groups_reconciled':dict(actual),'third_position_question_mark':{flag:sum(g['count'] for g in groups if g['residual_group']==flag and g['question_mark_mask'][2]=='?') for flag in ('gt5','le5')}}
receipt={'comparisons':out,'direction_audit_sha256':hashlib.sha256((p/'adhd_direction_audit.json').read_bytes()).hexdigest(),'limitations':'Reconciles aggregate group totals only; does not independently establish source encoding, subject overlap, or statistical interpretation.'}
(p/'direction_verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
