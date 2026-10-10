"""Only targeted actual-label and8+27 plan controls; no worker/body operation."""
import copy
import json
from pathlib import Path
import sys

P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'scripts'))
from core_large35_checkpoint_adoption_v1 import qc,semantic_snoring,FIRST8
from extension_replay_common_v3 import sha,write_new

plan_path=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/core_pipeline/large35_checkpoint_continuation_v1/core_large35_checkpoint_continuation_plan_v1.json')
PLAN_SHA='716b9df426dc5291296a4dbbedd9f057f615168c8c61945da4ccb7dd074f7bb0'
if sha(plan_path)!=PLAN_SHA:raise RuntimeError('FROZEN_PREPARED_PLAN_DIFFERS')
plan=json.loads(plan_path.read_text());old=json.loads(Path(plan['predecessor_plan']).read_text())
snoring=plan['adopted_members'][7];member=snoring['member']
metadata,steps=qc(member['harmonization_qc']);columns=json.loads(Path(snoring['columns_path']).read_text())
comparison=json.loads(Path(member['comparison_receipt']).read_text())
history=member['original_design']['historical_QC'];cases=[]


def test(name,trait='snoring',changes=None,reject=False):
    values=dict(history=copy.deepcopy(history),metadata=copy.deepcopy(metadata),steps=copy.deepcopy(steps),columns=copy.deepcopy(columns),comparison=copy.deepcopy(comparison))
    if changes:changes(values)
    try:semantic_snoring(trait,values['history'],values['metadata'],values['steps'],values['columns'],values['comparison'])
    except RuntimeError as e:
        if not reject:raise
        cases.append(dict(case=name,status='EXPECTED_REJECTION',reason=str(e)))
    else:
        if reject:raise RuntimeError('NEGATIVE_CONTROL_ACCEPTED: '+name)
        cases.append(dict(case=name,status='PASS'))


test('ACTUAL_SNORING_LABEL_ONLY_EQUIVALENCE')
test('OTHER_TRAIT_ALIAS_REJECTED',trait='sleep_apnea',reject=True)
test('ALTERED_DROP_COUNT_REJECTED',changes=lambda v:v['steps'][13].update(dropped=1),reject=True)
test('SOURCE_N_COLUMN_REJECTED',changes=lambda v:v['columns']['columns'].update(N='N'),reject=True)
test('ADDITIONAL_QC_REASON_CHANGE_REJECTED',changes=lambda v:v['steps'][12].update(reason='changed duplicate rule'),reject=True)
assert [x['trait_id'] for x in plan['adopted_members']]==FIRST8
assert [x['trait_id'] for x in plan['members']]==[x['trait_id'] for x in old['members'][8:]]
assert len(plan['adopted_members'])==8 and len(plan['members'])==27 and plan['expected_new_worker_count']==135
assert len(set(plan['original35_order']))==35 and FIRST8+[x['trait_id'] for x in plan['members']]==plan['original35_order']
assert plan['guard']==old['guard'] and plan['deadline_anchor_utc']==old['prepared_utc']
assert all(x['original_design']==y['original_design'] for x,y in zip(plan['members'],old['members'][8:]))
def relocated(value):
    if isinstance(value,dict):return {k:relocated(v) for k,v in value.items()}
    if isinstance(value,list):return [relocated(v) for v in value]
    if isinstance(value,str) and value.startswith(old['private_namespace']+'/'):
        return plan['private_namespace']+value[len(old['private_namespace']):]
    return value
assert plan['members']==[relocated(m) for m in old['members'][8:]]
assert plan['predecessor_failed_family_reclassified'] is False and snoring['failed_spool_preserved'] is True
cases.append(dict(case='EXACT8_ADOPTED27_NEW_ORIGINAL35_MAPPING_AND_NAMESPACE_ONLY_ARGV',status='PASS',mapping_assertions=8))
receipt=dict(status='FOCUSED_AUTHOR_DELTA_CONTROLS_PASS_NO_EXECUTION_ADMISSION',cases=cases,case_count=len(cases),
             plan_path=str(plan_path),plan_sha256=PLAN_SHA,adjudication_sha256=sha(plan['snoring_adjudication']),
             checker_sha256=sha(__file__),workers_launched=0,body_reads=0,production_mutex_operations=0,
             scope='One actual label equivalence, four focused negative controls, exact8+27 mapping; no old suite repeated.')
path=P/'reviews/core_large35_checkpoint_delta_controls_receipt_v1.json';write_new(path,receipt)
print(json.dumps(dict(receipt=str(path),sha256=sha(path),case_count=len(cases),status=receipt['status']),indent=2))
