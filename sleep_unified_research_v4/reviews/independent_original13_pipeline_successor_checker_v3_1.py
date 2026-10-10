#!/usr/bin/env python3
"""Narrow immutable v2→v3 correction controls; no production execution."""
import ast
import copy
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
B=load('_original13_inherited_review_v2_2',R/'independent_original13_pipeline_checker_v2_2.py')
B.ROOT=P.parents[1]/'independent_original13_pipeline_successor_fixture_v3_1'
OUT=R/'independent_original13_pipeline_successor_controls_v3_1.json'
if B.ROOT.exists() or OUT.exists():raise RuntimeError('OWN_NEW_REVIEW_NAMESPACE_REQUIRED')
B.ROOT.mkdir()
EXPECTED={'101_run_original_validation_pipeline_v3.py':'97fff8092a55da411b2dab3d311d0d6d821122d62bbc39eb1a9e77d36f55c757',
 '102_prepare_original_validation_pipeline_v3.py':'28ef4f568735d6538f669d875418a3df9a030093f6dbfb0fab66ce461301efac'}
for n,h in EXPECTED.items():B.check('exact_successor_'+n,B.sha(S/n)==h)
B.C=load('_original13_actual_successor_controller',S/'101_run_original_validation_pipeline_v3.py')
renames={
 'validation_pipeline_replay_v3':'validation_pipeline_replay_v2',
 'original13_validation_pipeline_execution_v3':'original13_validation_pipeline_execution_v2',
 'frozen_original13_validation_pipeline_plan_v3':'frozen_original13_validation_pipeline_plan_v2',
 'validation_pipeline_pending_v3.json':'validation_pipeline_pending_v2.json',
 'validation_pipeline_terminal_v3.json':'validation_pipeline_terminal_v2.json',
 'validation_pipeline_execution_receipt_v3.json':'validation_pipeline_execution_receipt_v2.json',
 'validation_pipeline_plan_v3.json':'validation_pipeline_plan_v2.json',
 '101_run_original_validation_pipeline_v3.py':'99_run_original_validation_pipeline_v2.py',
 'independent_original13_pipeline_prelaunch_v3.md':'independent_original13_pipeline_prelaunch_v2.md',
 'independent_original13_pipeline_prelaunch_v3.json':'independent_original13_pipeline_prelaunch_v2.json',
 'independent_original13_pipeline_prelaunch_seal_v3.json':'independent_original13_pipeline_prelaunch_seal_v2.json'}
class OnlyDeclaredChanges(ast.NodeTransformer):
 def visit_Constant(self,node):
  if isinstance(node.value,str) and node.value in renames:node.value=renames[node.value]
  return node
 def visit_For(self,node):
  if isinstance(node.target,ast.Tuple) and [n.id for n in node.target.elts]==['oldname','oldsha']:return None
  if isinstance(node.target,ast.Tuple) and [n.id for n in node.target.elts]==['sealname','sealdigest']:return None
  return self.generic_visit(node)
 def visit_If(self,node):
  if any(isinstance(n,ast.Constant) and n.value=='CURRENT_WORKER_FAILURE_FIELDS_VETO_PIPELINE' for n in ast.walk(node)):return None
  return self.generic_visit(node)
 def visit_BoolOp(self,node):
  node=self.generic_visit(node)
  if any(isinstance(n,ast.Constant) and n.value=='metadata_errors' for n in ast.walk(node)):
   kept=[]
   for value in node.values:
    if any(isinstance(n,ast.Constant) and n.value in ['metadata_errors','plan_unchanged','cleanup_error'] for n in ast.walk(value)):continue
    if isinstance(value,ast.Compare) and isinstance(value.ops[0],ast.NotEq) and isinstance(value.comparators[0],ast.List) and not value.comparators[0].elts and any(isinstance(n,ast.Constant) and n.value=='remaining_group_members' for n in ast.walk(value)):
     value=value.left
    kept.append(value)
   node.values=kept
  return node
def ast_equal(old,new):
 left=ast.parse((S/old).read_text());right=OnlyDeclaredChanges().visit(ast.parse((S/new).read_text()))
 return ast.dump(left,include_attributes=False)==ast.dump(right,include_attributes=False)
B.check('101_entire_AST_equal_after_only_declared_version_and_failure_gates',ast_equal('99_run_original_validation_pipeline_v2.py','101_run_original_validation_pipeline_v3.py'))
B.check('102_entire_AST_equal_after_only_declared_version_and_seal_binding_loops',ast_equal('98_prepare_original_validation_pipeline_v2.py','102_prepare_original_validation_pipeline_v3.py'))
# Reuse the sealed mock harness with a transparent route adaptation. This does
# not alter either candidate or execute the production preparer/worker.
oldfixture=B.fixture
def fixture(name):
 d,ssd,oldout,plan,a=oldfixture(name);newout=ssd/'validation_pipeline_replay_v3';oldout.rename(newout)
 def route(value):
  if isinstance(value,dict):return {k:route(v) for k,v in value.items()}
  if isinstance(value,list):return [route(v) for v in value]
  if isinstance(value,str):return value.replace('validation_pipeline_replay_v2','validation_pipeline_replay_v3').replace('validation_pipeline_pending_v2.json','validation_pipeline_pending_v3.json').replace('validation_pipeline_terminal_v2.json','validation_pipeline_terminal_v3.json').replace('validation_pipeline_execution_receipt_v2.json','validation_pipeline_execution_receipt_v3.json')
  return value
 plan=route(plan)
 plan['dependencies_sha256'].update({str(S/n):h for n,h in EXPECTED.items()})
 B.rewrite(a.plan,plan);a.plan_sha=B.sha(a.plan)
 adm=json.loads(a.admission.read_text());adm.update(plan_sha256=a.plan_sha,executor_sha256=EXPECTED['101_run_original_validation_pipeline_v3.py'])
 B.rewrite(a.admission,adm);a.admission_sha=B.sha(a.admission)
 return d,ssd,newout,plan,a
B.fixture=fixture
source=inspect.getsource(B.controller_case).replace("sha(S/'99_run_original_validation_pipeline_v2.py')","sha(S/'101_run_original_validation_pipeline_v3.py')")
exec(compile(source,'<sealed review harness successor SHA route only>','exec'),B.__dict__)
committed,error=B.controller_case('healthy13_39')
B.check('healthy13_39_candidate_v3_commits',committed and error is None and B.MAP['healthy13_39']['mock_worker_count']==39)
for name in ['metadata_error','plan_unchanged_false','teardown_cleanup_error']:
 committed,error=B.controller_case(name,name)
 B.check('three_measured_v2_contradictions_now_rejected_'+name,not committed and error is not None and B.MAP[name]['mock_worker_count']==1)
# Exercise the actual final oracle separately; immediate rejections should not
# obscure that terminal consumption also checks the same evidence.
tree=ast.parse((S/'101_run_original_validation_pipeline_v3.py').read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run')
identity=next(n for n in run.body if isinstance(n,ast.FunctionDef) and n.name=='identity')
oracle=next(n for n in ast.walk(identity) if isinstance(n,ast.If) and any(isinstance(x,ast.Constant) and x.value=='EXACT_COMPLETE_REAPED39_WORKERS_REQUIRED' for x in ast.walk(n)))
compiled=compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(oracle)],type_ignores=[])),'<actual v3 final worker oracle>','exec')
command=['OWN_INVENTED'];base={'status':'WORKER_COMPLETE_VERIFIED','command':command,'plan_sha256':'a'*64,'returncode':0,'stop_reason':None,'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[]},'metadata_errors':[],'plan_unchanged':True}
def final(r):exec(compiled,{'r':r,'expected_workers':{'owned':command},'path':'owned','a':type('Args',(),{'plan_sha':'a'*64})()})
final(base);B.check('actual_final_oracle_healthy',True)
for label,key,value in [('metadata_error','metadata_errors',['x']),('null_metadata','metadata_errors',None),('plan_false','plan_unchanged',False),('plan_numeric_one','plan_unchanged',1)]:
 r=copy.deepcopy(base);r[key]=value;B.expect_fail('actual_final_oracle_'+label,lambda:final(r))
for label,value in [('null_group',None),('empty_string_group','')]:
 r=copy.deepcopy(base);r['process_group_teardown']['remaining_group_members']=value;B.expect_fail('actual_final_oracle_'+label,lambda:final(r))
r=copy.deepcopy(base);r['process_group_teardown']['cleanup_error']='x';B.expect_fail('actual_final_oracle_cleanup_error',lambda:final(r))
# Execute ONLY the two additive dependency loops, not prepare(). The callbacks
# below create metadata faults by returning altered digests or parsed values;
# no sealed evidence is mutated.
prep=ast.parse((S/'102_prepare_original_validation_pipeline_v3.py').read_text());fn=next(n for n in prep.body if isinstance(n,ast.FunctionDef) and n.name=='prepare')
loops=[n for n in fn.body if isinstance(n,ast.For) and isinstance(n.target,ast.Tuple) and [x.id for x in n.target.elts] in [['oldname','oldsha'],['sealname','sealdigest']]]
B.check('exact_two_additive_preparer_binding_loops',len(loops)==2)
loopcode=compile(ast.fix_missing_locations(ast.Module(body=copy.deepcopy(loops),type_ignores=[])),'<actual v3 additive binding loops>','exec')
def bind(fault=None):
 dependencies={};state={'post_parse':False};target=R/'independent_original13_pipeline_prelaunch_seal_v2.json'
 artifact=R/'independent_original13_pipeline_checker_v2_2.py'
 if fault=='dependency_conflict':dependencies[str(artifact)]='0'*64
 def digest(path):
  if Path(path)==target and (fault=='preparse_seal' or fault=='postparse_seal' and state['post_parse']):return '0'*64
  if Path(path)==artifact and fault=='artifact_changed':return '0'*64
  return B.sha(path)
 def parsed(text):
  value=json.loads(text)
  if value.get('schema')=='independent_original13_pipeline_candidate_review_seal_v2':state['post_parse']=True
  return value
 exec(loopcode,{'P':P,'dependencies':dependencies,'sha':digest,'regular':B.C.regular,'json':type('JSON',(),{'loads':staticmethod(parsed)})})
 return dependencies
bindings=bind();B.check('all_additive_old_rejection_and_adapter_science_bindings_accepted',True)
for fault in ['preparse_seal','postparse_seal','artifact_changed','dependency_conflict']:B.expect_fail('preparer_'+fault,lambda:bind(fault))
# Verify full old and scientific seals, all direct maps, and the unchanged
# helper/candidate hashes needed to inherit the668 previously sealed assertions.
preserved={}
for sealname,h in [('independent_original13_pipeline_prelaunch_seal_v2.json','dc3c560ada5086c778d51af4063091447846ad7269963ac6efbe2996a5cc4730'),('validation_collector_adapter_review_seal_v2.json','8938220336f9cb35a3e3396bff459c05093cf72add56dc3c8a8d757eb1a1753f')]:
 path=R/sealname;B.check('exact_preserved_'+sealname,B.sha(path)==h);preserved[str(path)]=h
 sealed=json.loads(path.read_text())
 for mapping in sealed.values():
  if isinstance(mapping,dict) and mapping and all(isinstance(v,str) and len(v)==64 for v in mapping.values()):
   for path,h in mapping.items():
    B.check('preserved_direct_map_'+Path(path).name,B.sha(path)==h);preserved[path]=h
for n,h in EXPECTED.items():B.check('successor_unchanged_after_controls_'+n,B.sha(S/n)==h)
rss,size=B.bounded()
receipt={'schema':'independent_original13_pipeline_successor_controls_v3_1','status':'PASS_NARROW_ADDITIVE_CORRECTION','checks':B.CHECKS,'check_count':len(B.CHECKS),'mock_cases':B.MAP,
 'candidate_sha256':{str(S/n):h for n,h in EXPECTED.items()},'actual_additive_preparer_dependency_map':bindings,'verified_preserved_direct_map_sha256':preserved,
 'peak_RSS_bytes':rss,'synthetic_fixture_bytes':size,'inherited_old_assertions':668,'inherited_checker_sha256':B.sha(R/'independent_original13_pipeline_checker_v2_2.py'),
 'scope':'One39-command healthy mock and three one-worker contradictions, actual final oracle AST and isolated additive dependency loops. Other668 controls inherited only for AST-unchanged functions/helper identities, with original qualifications retained.',
 'production_body_reads':0,'network_calls':0,'native_workers':0,'real_heavy_mutex_operations':0,'execution_admitted':False}
B.save(OUT,receipt)
print(json.dumps({'receipt':str(OUT),'sha256':B.sha(OUT),'checks':len(B.CHECKS),'peak_RSS_bytes':rss}))
