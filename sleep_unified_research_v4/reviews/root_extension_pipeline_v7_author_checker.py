#!/usr/bin/env python3
"""Small author controls; no source, worker, runtime or estimator operations."""
import ast,copy,hashlib,json
from pathlib import Path
from types import SimpleNamespace

P=Path(__file__).resolve().parents[1]
ROOT=P.parents[1]
OUT=ROOT/'extension_v7_author_controls_v1'
assert not OUT.exists()
OUT.mkdir()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
code=P/'scripts/47_run_extension_pipeline_replay_v7.py'
tree=ast.parse(code.read_text())
fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='read_fixed_worker')
namespace={'Path':Path,'json':json,'sha':sha}
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(code),'exec'),namespace)
read=namespace['read_fixed_worker']
result=OUT/'metadata_result.json';result.write_text('{"status":"INVENTED_METADATA_ONLY"}\n')
receipt=OUT/'worker.json'
healthy={'status':'WORKER_COMPLETE_VERIFIED','plan_sha256':'invented_plan','command':['invented_no_launch'],'returncode':0,'stop_reason':None,'metadata_errors':[],'plan_unchanged':True,'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[]},'output_sha256':{str(result):sha(result)}}
controls=[]
def case(name,change,should_pass):
 r=copy.deepcopy(healthy);change(r);receipt.write_text(json.dumps(r)+'\n')
 digest=sha(receipt)
 try:read(receipt,digest,healthy['command'],healthy['plan_sha256']);accepted=True;reason=None
 except RuntimeError as e:accepted=False;reason=str(e)
 assert accepted is should_pass,(name,accepted)
 controls.append({'name':name,'accepted':accepted,'expected_accepted':should_pass,'reason':reason})
case('healthy',lambda r:None,True)
case('metadata_errors',lambda r:r.update(metadata_errors=['invented_error']),False)
case('plan_unchanged_false',lambda r:r.update(plan_unchanged=False),False)
case('cleanup_error',lambda r:r['process_group_teardown'].update(cleanup_error='invented_error'),False)
case('wrong_command',lambda r:r.update(command=['different_no_launch']),False)
case('returncode17',lambda r:r.update(returncode=17),False)
case('stop_reason_nonnull',lambda r:r.update(stop_reason='invented_stop'),False)
case('missing_stop_reason',lambda r:r.pop('stop_reason'),False)
case('unreaped',lambda r:r['process_group_teardown'].update(remaining_group_members=[123456789]),False)
case('owned_false',lambda r:r.update(owned_cleanup_verified=False),False)
receipt.write_text(json.dumps(healthy)+'\n');digest=sha(receipt)
def mutate_after_parse(text):
 r=json.loads(text);receipt.write_text(json.dumps({**healthy,'extra':'postparse_mutation'})+'\n');return r
namespace['json']=SimpleNamespace(loads=mutate_after_parse)
try:read(receipt,digest,healthy['command'],healthy['plan_sha256']);raise AssertionError('postparse mutation accepted')
except RuntimeError as e:
 controls.append({'name':'postparse_receipt_mutation','accepted':False,'expected_accepted':False,'reason':str(e)})
namespace['json']=json
receipt.write_text(json.dumps(healthy)+'\n');digest=sha(receipt)
failure=Path(str(receipt)+'.failure.json');failure.write_text('{}\n')
try:read(receipt,digest,healthy['command'],healthy['plan_sha256']);raise AssertionError('failure addendum accepted')
except RuntimeError as e:
 controls.append({'name':'failure_addendum','accepted':False,'expected_accepted':False,'reason':str(e)})
oldplan=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/extension_pipeline_replay_v6/extension_pipeline_replay_plan_v6.json')
newplan=oldplan.parent.parent/'extension_pipeline_replay_v7/extension_pipeline_replay_plan_v7.json'
assert sha(oldplan)=='d0fd223d468f26583bfd2bd477a51548646b57e91275d7733eb1ae814b94b798'
assert sha(newplan)=='ad0c08bc0b0a62c8a1bc070221689f2733b0eefe8f25de8cc1bd8689a81b1ba5'
a=json.loads(oldplan.read_text());b=json.loads(newplan.read_text())
assert len(a['members'])==len(b['members'])==100
for old,new in zip(a['members'],b['members']):
 normalized={}
 for key,value in new.items():
  if isinstance(value,str):value=value.replace(str(newplan.parent)+'/',str(oldplan.parent)+'/')
  elif isinstance(value,list):value=[v.replace(str(newplan.parent)+'/',str(oldplan.parent)+'/') if isinstance(v,str) else v for v in value]
  normalized[key]=value
 assert old==normalized,old['extension_trait_id']
assert a['guard']==b['guard'] and b['total_owned_commands']==400 and b['estimator_calls']==0
oldtree=ast.parse((P/'scripts/47_run_extension_pipeline_replay_v6.py').read_text())
unchanged=[]
for oldfn in oldtree.body:
 if isinstance(oldfn,ast.FunctionDef) and oldfn.name not in ['run','seal_master','main']:
  newfn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name==oldfn.name)
  assert ast.dump(oldfn,include_attributes=False)==ast.dump(newfn,include_attributes=False),oldfn.name
  unchanged.append(oldfn.name)
oldrun=next(x for x in oldtree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
newrun=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='run')
oldscience=next(x for x in ast.walk(oldrun) if isinstance(x,ast.FunctionDef) and x.name=='validation')
newscience=next(x for x in ast.walk(newrun) if isinstance(x,ast.FunctionDef) and x.name=='validation')
assert ast.dump(oldscience,include_attributes=False)==ast.dump(newscience,include_attributes=False)
report={'status':'AUTHOR_CONTROLS_PASS_SEPARATE_INDEPENDENT_REVIEW_REQUIRED','code_sha256':sha(code),'plan_sha256':sha(newplan),'helper_cases':controls,'helper_cases_count':len(controls),'all100_member_values_identical_after_private_namespace_normalization':True,'guard_identical':True,'unchanged_top_level_function_AST':unchanged,'original_validation_command_AST_identical':True,'scientific_validator_sha256':sha(P/'scripts/extension_replay_validate_v6.py'),'workers_launched':0,'source_reference_runtime_payloads_read':0,'estimator_calls':0,'execution_admission_granted':False}
target=P/'reviews/root_extension_pipeline_v7_author_controls_v1.json'
assert not target.exists();target.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'receipt':str(target),'sha256':sha(target),'helper_cases':len(controls),'science_members':100,'status':report['status']}))
