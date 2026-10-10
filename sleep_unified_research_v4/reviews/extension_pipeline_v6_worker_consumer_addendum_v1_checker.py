#!/usr/bin/env python3
"""Targeted new consumer concern; exact AST, no processing replay or worker."""
import ast
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import sys
import time
from types import SimpleNamespace
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';sys.path.insert(0,str(S))
BASE='extension_pipeline_v6_worker_consumer_addendum_v1'
F=P.parents[1]/(BASE+'_fixtures');RECEIPT=R/(BASE+'_controls.json')
SOURCE=S/'47_run_extension_pipeline_replay_v6.py';PIN='35ae2a55ff3523442c28d097df018771cc505e52ca29c818fe6a67b37a53a701'
OLDSEAL=R/'genomicsem_extension_pipeline_v6_preflight_seal.json';OLDPIN='ce945e2289da6ce73e75899e2c54ac5977e0ea6def639a5e48c31e20b5358831'
CHECKS=[];CASES=[];HASHES={};START=time.monotonic()
def check(x,label):
 CHECKS.append({'label':label,'passed':bool(x)})
 if not x:raise AssertionError(label)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def hashed(p,expected):
 check(Path(p).is_file() and not Path(p).is_symlink(),'regular:'+str(p));h=sha(p);check(h==expected,'SHA:'+str(p));HASHES[str(p)]=h;return h
def save(p,v):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 return sha(p)
def module():
 sp=importlib.util.spec_from_file_location('_independent_extension6_targeted',SOURCE);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
tree=ast.parse(SOURCE.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run')
regular=next(n for n in run.body if isinstance(n,ast.FunctionDef) and n.name=='regular_output_hashes')
final=next(n for n in run.body if isinstance(n,ast.FunctionDef) and n.name=='final_identity_gate')
workerloop=next(n for n in final.body if isinstance(n,ast.For) and isinstance(n.iter,ast.Call) and isinstance(n.iter.func,ast.Attribute) and isinstance(n.iter.func.value,ast.Subscript) and isinstance(n.iter.func.value.slice,ast.Constant) and n.iter.func.value.slice.value=='worker_receipt_sha256')
LOOP_SHA=hashlib.sha256(ast.dump(workerloop,include_attributes=False).encode()).hexdigest()
FUNCTION_SHA=hashlib.sha256(ast.dump(final,include_attributes=False).encode()).hexdigest()
def make_record(plan_sha,output):
 return {'status':'WORKER_COMPLETE_VERIFIED','command':['OWN_DECLARED_EXPECTED_COMMAND'],'plan_sha256':plan_sha,'returncode':0,'stop_reason':None,'plan_unchanged':True,'owned_cleanup_verified':True,'process_group_teardown':{'initial_group_members':[],'remaining_group_members':[]},'metadata_errors':[],'output_sha256':{str(output):sha(output)},'fixture_only_no_worker':True}
def contradiction(r,field):
 if field=='metadata_errors':r[field]=['EXPLICIT_FINAL_METADATA_ERROR']
 elif field=='plan_unchanged':r[field]=False
 elif field=='cleanup_error':r['process_group_teardown'][field]='EXPLICIT_TEARDOWN_ERROR'
 elif field=='command':r[field]=['WRONG_COMMAND']
 elif field=='returncode':r[field]=17
 elif field=='stop_reason':r[field]='EXPLICIT_WORKER_STOP'
 else:raise AssertionError(field)
def unit_case(name,field=None,position=0,ordinary=None):
 d=F/name;d.mkdir();output=d/'owned_stdout.fixture';output.write_bytes(b'OWN TINY OUTPUT')
 pp=d/'plan.json';save(pp,{'fixture':name});a=SimpleNamespace(plan_sha=sha(pp));mapping={};events=[]
 for index in range(2):
  path=d/('worker_'+str(index)+'.json');r=make_record(a.plan_sha,output)
  if field and index==position:contradiction(r,field)
  if ordinary and index==position:
   if ordinary=='status':r['status']='WORKER_FAILED_PRESERVED'
   if ordinary=='plan':r['plan_sha256']='f'*64
   if ordinary=='owned':r['owned_cleanup_verified']=False
   if ordinary=='remaining':r['process_group_teardown']['remaining_group_members']=[{'fixture_only':True}]
  mapping[str(path)]=save(path,r)
 master={'worker_receipt_sha256':mapping}
 def trace_sha(path):events.append(('hash',str(path)));return sha(path)
 def trace_parse(text):
  events.append(('parse',None));return json.loads(text)
 env={'Path':Path,'sha':trace_sha,'json':SimpleNamespace(loads=trace_parse),'master':master,'a':a}
 isolated=ast.FunctionDef(name='exact_worker_consumer_loop',args=ast.arguments(posonlyargs=[],args=[],kwonlyargs=[],kw_defaults=[],defaults=[]),body=[copy.deepcopy(workerloop)],decorator_list=[])
 exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(regular),isolated],type_ignores=[])),str(SOURCE),'exec'),env)
 if ordinary=='preparse_sha_drift':
  with Path(next(iter(mapping))).open('a') as f:f.write(' ')
 try:env['exact_worker_consumer_loop']();accepted=True;error=None
 except BaseException as e:accepted=False;error=type(e).__name__+': '+str(e)
 expected=ordinary is None
 check(accepted==expected,name+':measured_consumer_outcome')
 for i,event in enumerate(events):
  if event[0]=='parse':check(i>0 and events[i-1][0]=='hash',name+':fixed_receipt_hash_before_first_parse')
 CASES.append({'case':name,'field':field,'bad_record_position':position,'isolated_exact_AST_loop_accepts':accepted,'error':error,'events':events,'actual_worker':False})

def terminal_cases(runner):
 # One set of 100 invented members and 400 tiny receipts, not 400 workers.
 # Only three actual final-gate/Terminal2 paths consume this declared fixture.
 home=F/'terminal_consumed';home.mkdir();pp=home/'plan.json';ap=home/'admission.json';save(pp,{'fixture_only':True});save(ap,{'fixture_only_not_admission':True})
 a=SimpleNamespace(plan=pp,admission=ap,plan_sha=sha(pp),admission_sha=sha(ap))
 plan={'members':[]};master={'worker_receipt_sha256':{},'completed_members':[],'checkpoint_binding_sha256':{},'discovered_acquisition_receipt_sha256':{},'acquisition_family_terminal_evidence':{'explicit_mock':True},'historical190_gate_receipts':{'explicit_mock':True}}
 for index in range(100):
  trait='owned_%03d'%index;member={'extension_trait_id':trait}
  for key in ['harmonized','harmonization_qc','harmonization_receipt','munged','source_gate_receipt','comparison_receipt']:
   path=home/'outputs'/(trait+'__'+key+'.fixture');path.parent.mkdir(exist_ok=True);path.write_bytes(b'PRIVATE INVENTED METADATA\n');member[key]=str(path)
  plan['members'].append(member);bp=home/'receipts'/(trait+'.checkpoint_binding.json');master['checkpoint_binding_sha256'][str(bp)]=save(bp,{'fixture_only':True,'trait':trait})
  output_map={member[k]:sha(member[k]) for k in member if k!='extension_trait_id'}
  master['completed_members'].append({'trait':trait,'output_sha256':output_map,'comparison_sha256':sha(member['comparison_receipt']),'new_munged_sha256':sha(member['munged'])})
  for label in ['source','harmonize','munge','compare']:
   path=home/'workers'/(trait+'__'+label+'.worker.json');record=make_record(a.plan_sha,Path(member['harmonized']));master['worker_receipt_sha256'][str(path)]=save(path,record)
 first=Path(next(iter(master['worker_receipt_sha256'])));last=Path(next(reversed(master['worker_receipt_sha256'])));ownership=[True,[]]
 protected=SimpleNamespace(stage_resource_gate=lambda *a,**k:{'resource_and_science_mock_only':True},check_termination=lambda *a:None,TERMINATION_REQUEST=[],safe_print=lambda *a,**k:None,baseline_gate=lambda p:{'explicit_mock':True})
 for name,target in [('healthy_terminal',None),('contradictory_first_terminal',first),('contradictory_final_terminal',last)]:
  d=F/name;d.mkdir();pending=d/'pending.json';seal=d/'seal.json';receipt=d/'primary.json';work=copy.deepcopy(master)
  changed=None
  if target:
   changed=target.read_bytes();record=json.loads(changed)
   for field in ['metadata_errors','plan_unchanged','cleanup_error','command','returncode','stop_reason']:contradiction(record,field)
   # Same ordering as the actual immediate helper return: write the receipt,
   # then hash stamp it into the immutable master before any consumer parse.
   target.write_text(json.dumps(record,indent=2)+'\n');work['worker_receipt_sha256'][str(target)]=sha(target)
  work['status']='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS';work['termination_requests']=[]
  env={'Path':Path,'sha':sha,'json':json,'ownership':ownership,'physical_mount':lambda:None,'check_bindings':lambda *a,**k:None,'plan':plan,'a':a,'admission':{},'admission_gate':lambda *a:None,'acquisition_family_gate':lambda *a,**k:{'explicit_mock':True},'master':work,'protected':protected,'OUT':home,'checkpoint_binding_gate':lambda *a:None,'__file__':str(SOURCE)}
  exec(compile(ast.fix_missing_locations(ast.Module(body=[copy.deepcopy(regular),copy.deepcopy(final)],type_ignores=[])),str(SOURCE),'exec'),env)
  terminal=runner.TerminalCommit(pending,seal,{'fixture_only_terminal':True})
  runner.seal_master(protected,work,receipt,plan,0,ownership,terminal,env['final_identity_gate'])
  import terminal_commit_common_v2 as tc
  try:tc.require_committed(pending,seal,{'fixture_only_terminal':True},{str(receipt):sha(receipt)});consumer=True
  except BaseException:consumer=False
  committed=work['status']=='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS' and not pending.exists() and seal.is_file()
  check(committed and consumer,name+':full400_metadata_terminal_accepted')
  check(len(work['worker_receipt_sha256'])==400 and len(work['completed_members'])==100,name+':literal400and100_fixture_cardinality')
  CASES.append({'case':name,'status':work['status'],'pending_absent':not pending.exists(),'terminal_consumer_accepts':consumer,'worker_receipt_count':400,'invented_completed_member_count':100,'contradictory_path':str(target) if target else None,'all_six_fields_present_and_ignored':bool(target),'actual_pipeline_run_called':False,'actual_worker_calls':0,'mocked_dependencies':['producer','historical190','mount','source','resource','admission']})
  if changed is not None:target.write_bytes(changed)

def main():
 check(not F.exists() and not RECEIPT.exists(),'fresh_private_fixture_namespace');F.mkdir()
 hashed(SOURCE,PIN);hashed(OLDSEAL,OLDPIN);old=json.loads(OLDSEAL.read_text())
 for p,r in old['artifacts'].items():hashed(p,r['sha256']);check(Path(p).stat().st_size==r['bytes'],'old_artifact_size:'+p)
 review=json.loads((R/'genomicsem_extension_pipeline_v6_preflight.json').read_text())
 for p,h in review['source_sha256'].items():hashed(p,h)
 hashed(review['plan_path'],review['plan_sha256'])
 runner=module()
 immediate=next(n for n in ast.walk(run) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute) and n.value.func.attr=='worker')
 check(immediate.value.func.value.id=='protected','immediate_call_exact_protected_worker')
 check('metadata_errors' not in ast.unparse(final) and 'plan_unchanged' not in ast.unparse(final) and 'cleanup_error' not in ast.unparse(final),'failure_fields_absent_from_final_consumer_AST')
 for name,field,ordinary in [('healthy_unit',None,None),('wrong_status',None,'status'),('wrong_plan',None,'plan'),('cleanup_unverified',None,'owned'),('remaining_group',None,'remaining'),('receipt_preparse_drift',None,'preparse_sha_drift')]:unit_case(name,field,ordinary=ordinary)
 for field in ['metadata_errors','plan_unchanged','cleanup_error','command','returncode','stop_reason']:
  for position in [0,1]:unit_case(field+'_position'+str(position),field,position)
 terminal_cases(runner)
 for p,h in HASHES.items():check(sha(p)==h,'preserved_final_SHA:'+p)
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;elapsed=time.monotonic()-START
 files={str(p):sha(p) for p in F.rglob('*') if p.is_file() and not p.is_symlink()};size=sum(Path(p).stat().st_size for p in files)
 check(rss<128*2**20 and size<16*2**20 and elapsed<180,'bounded_review_limits')
 record={'schema':BASE,'status':'MEASURED_SIX_WORKER_CONSUMER_OMISSIONS_REJECT_UNLAUNCHED_EXTENSION6','check_count':len(CHECKS),'checks':CHECKS,'case_count':len(CASES),'cases':CASES,'exact_worker_loop_AST_sha256':LOOP_SHA,'exact_final_identity_gate_AST_sha256':FUNCTION_SHA,'reviewed_metadata_and_old_proof_sha256':HASHES,'old_review_seal_sha256':OLDPIN,'runner_sha256':PIN,'checker_sha256':sha(__file__),'fixture_regular_file_sha256':files,'private_fixture_bytes':size,'peak_rss_bytes':rss,'elapsed_seconds':elapsed,'normal_producer_contract_qualification':'Fresh actual worker sets WORKER_COMPLETE_VERIFIED only for returncode0, stop_reasonNone and unchanged plan; metadata errors call record_error and set stop, exact input command is serialized. The actual30 teardown helper returns no cleanup_error field and raises on failure. Deliberately contradictory persisted records are not shown to arise normally; completed native campaign validity is unchanged. Prior-record helper checkpoint itself checks status/plan/command but not all omitted flags.','receipt_hash_before_first_parse':True,'new_broad_replay_or_old402_controls_rerun':False,'old1729_fixture_bodies_rehashed':False,'full_pipeline_run_called':False,'actual_scientific_workers':0,'actual_mutex':0,'GWAS_reference_reads':0,'network':0,'execution_admission':False,'producer_or_source_or_candidate_mutations':0}
 save(RECEIPT,record);print(json.dumps({'receipt_sha256':sha(RECEIPT),'cases':len(CASES),'assertions':len(CHECKS),'elapsed_seconds':elapsed,'peak_rss_bytes':rss,'fixture_bytes':size}))
if __name__=='__main__':main()
