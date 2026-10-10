#!/usr/bin/env python3
"""Narrow tiny-only original47 adapter96 correction review; no production IO."""
import ast
import contextlib
import importlib.util
import json
from pathlib import Path
import sys
import types

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('inherited_fixture_helpers',HERE/'validation_collector_adapter_checker_v1_1.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
c.FIX=HERE/'validation_collector_adapter_fixtures_v2'
c.ADAPTER=c.S/'96_replay_original_validation_collector_v2.py'
c.EXPECTED={**{p:h for p,h in c.EXPECTED.items() if '92_replay' not in p},str(c.ADAPTER):'3cf4d17f86ffe1013c03045353ab2e14e1452d107bc036f0b911a1c4d49b057d'}
original_fixture=c.fixture

def bind_contract(f,family_edits=None):
 """Seal only this test's intended metadata with real Terminal2; never acquisition."""
 contract=f['plan']['acquisition_terminal'];pairs=contract['source_receipt_pairs'];root=f['root']
 # Rebind semantic negative fixtures so they reach semantic gates, not old hashes.
 for pair in pairs:
  pair['sha256']=c.sha(pair['primary_path'])
 sources=[dict(path=x['primary_path'],sha256=x['sha256']) for x in pairs]
 n=f.get('terminal_number',0);f['terminal_number']=n+1
 folder=root/('terminal_'+str(n));folder.mkdir()
 family=dict(status='ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED',completed_source_count=13,
  plan_sha256=f['plan']['acquisition_plan_sha256'],source_receipts=sources)
 if family_edits:family_edits(family)
 primary=c.save(folder/'primary.json',family);mirror=c.save(folder/'mirror.json',family)
 expected={str(primary):c.sha(primary),str(mirror):c.sha(mirror)}
 binding=dict(plan_sha256=f['plan']['acquisition_plan_sha256'],synthetic_only=True,fixture=root.name)
 terminal=c.load('synthetic_terminal_'+root.name+'_'+str(n),c.TERMINAL)
 pending=folder/'PENDING';seal=folder/'seal.json';tc=terminal.TerminalCommit(pending,seal,binding)
 c.assert_(tc.commit(expected,identity_gate=lambda:None,resource_gate=lambda:None,termination_gate=lambda:None),'real_synthetic_terminal2_'+root.name+'_'+str(n))
 paths=[x[k] for x in pairs for k in ['primary_path','mirror_path']]
 contract.update(individual_receipt_paths=paths,pending_path=str(pending),seal_path=str(seal),binding=binding,primary_receipt_sha256=expected,
  metadata_sha256={p:c.sha(p) for p in [*paths,str(primary),str(mirror),str(seal)]})

def fixture(name,kind='finn',transform=None):
 f=original_fixture(name,kind,transform);plan=f['plan'];root=f['root'];member=f['member']
 plan['output_namespace']=str(root);(root/'receipts_v4').mkdir();member['adapter_receipt']=str(root/'receipts_v4'/(member['source_id']+'.adapter.json'))
 binary=c.write(root/'runtime/python3.9','SYNTHETIC_READ_ONLY_BINARY_IDENTITY_NOT_AN_EXECUTABLE\n');binary.chmod(0o444)
 logical=binary.with_name('python');logical.symlink_to('python3.9')
 plan['python']=str(logical);plan['collector_runtime_profile']=dict(logical_executable=str(logical),resolved_binary=str(binary),resolved_binary_sha256=c.sha(binary),
  declared_symlinks=[dict(path=str(logical),target='python3.9')])
 plan['dependencies_sha256'][str(binary)]=c.sha(binary)
 pairs=[]
 for i,m in enumerate(plan['members']):
  rawmember=dict(source_id=m['source_id'],body_path=str(f['body']) if i==0 else str(root/('unused_'+str(i)+'.gz')),
   expected_bytes=member['original_source_identity']['bytes'],expected_md5=member['original_source_identity']['md5'],expected_sha256=member['original_source_identity']['sha256'])
  receipt=dict(member=rawmember,status='EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED',plan_sha256=plan['acquisition_plan_sha256'],returncode=0,stop_reason=None,
   teardown=dict(teardown_verified=True,remaining_group_members=[],cleanup_error=None))
  primary=c.save(root/'producer_receipts'/('source'+str(i)+'.json'),receipt);mirror=c.save(root/'producer_mirrors'/('source'+str(i)+'.json'),receipt)
  pairs.append(dict(member=rawmember,primary_path=str(primary),mirror_path=str(mirror),sha256=c.sha(primary)))
 plan['acquisition_terminal']['source_receipt_pairs']=pairs;bind_contract(f);c.update(f);return f
c.fixture=fixture

def reject(name,mutate=None,hook=None,expected=None,sidecar_may_exist=False,outputs_may_exist=False):
 f=fixture(name)
 if mutate:mutate(f);c.update(f)
 err,_=c.adapter_run(f,hook);c.assert_(err is not None,name+'_reject')
 if expected:c.assert_(expected in err,name+'_exact_error')
 if not sidecar_may_exist:c.assert_(not Path(f['member']['adapter_receipt']).exists(),name+'_no_sidecar')
 if not outputs_may_exist:c.assert_(not Path(f['member']['new_munged']).exists(),name+'_no_promoted_output')
 c.cases.append(dict(name=name,status='PASS_REJECTION',error=err,provisional_sidecar_present=Path(f['member']['adapter_receipt']).exists()))
 return f

def change_receipt(f,edit,index=0):
 pair=f['plan']['acquisition_terminal']['source_receipt_pairs'][index]
 receipt=json.loads(Path(pair['primary_path']).read_text());edit(receipt)
 c.save(pair['primary_path'],receipt);c.save(pair['mirror_path'],receipt);bind_contract(f)

def runtime_edit(f,edit):edit(f['plan']['collector_runtime_profile'])

def main():
 c.assert_(not c.FIX.exists(),'fresh_v2_fixture_namespace');c.FIX.mkdir();sys.path.insert(0,str(c.S))
 for p,h in c.EXPECTED.items():c.assert_(c.sha(p)==h,'pinned_before_'+p)
 prior=c.S/'92_replay_original_validation_collector_v1.py'
 c.assert_(c.sha(prior)=='15253ace8479334198f2bff06f33a9a5afd690f0c47d3f0ff31a4ad9f372e4fc','preserved_92_source')
 trees={k:ast.parse(p.read_text()) for k,p in [('v1',prior),('v2',c.ADAPTER)]}
 for name in ['hashes','regular','fixed_hashes','load']:
  values=[next(n for n in t.body if isinstance(n,ast.FunctionDef) and n.name==name) for t in trees.values()]
  c.assert_(ast.dump(values[0],include_attributes=False)==ast.dump(values[1],include_attributes=False),'unchanged_helper_AST_'+name)
 for name in ['current_body','scoped_Path','frozen_head','local_stream']:
  values=[next(n for n in ast.walk(t) if isinstance(n,ast.FunctionDef) and n.name==name) for t in trees.values()]
  c.assert_(ast.dump(values[0],include_attributes=False)==ast.dump(values[1],include_attributes=False),'unchanged_scientific_interface_AST_'+name)
 oldurlopen=c.urllib.request.urlopen
 def no_network(*args,**kwargs):raise AssertionError('NETWORK_FORBIDDEN')
 c.urllib.request.urlopen=no_network
 try:
  c.healthy('healthy_finn_v2','finn');c.healthy('healthy_mvp_v2','mvp')
  # G1: closed route before any body hashing/collector import.
  touched=[]
  def observe_body(a,f):
   original=a.hashes
   def hashes(path):
    if str(path)==str(f['body']):touched.append('body_hash')
    return original(path)
   a.hashes=hashes
  f=reject('G1_outside_wrapper_route',lambda f:f['member'].__setitem__('adapter_receipt',str(f['root']/'outside_workspace_adapter.json')),hook=observe_body,expected='EXACT_PRIVATE_ADAPTER_RECEIPT_ROUTE_REQUIRED')
  c.assert_(not touched,'G1_rejects_before_body_hash_or_parse');c.assert_(not (f['root']/'outside_workspace_adapter.json').exists(),'G1_no_outside_write')
  def prior_wrapper(f):c.write(f['member']['adapter_receipt'],'preserved_prior')
  f=reject('G1_prior_wrapper',prior_wrapper,expected='PRIOR_PRIVATE_VALIDATION_OUTPUT_PRESERVED_NO_RETRY',sidecar_may_exist=True)
  c.assert_(Path(f['member']['adapter_receipt']).read_text()=='preserved_prior','G1_prior_wrapper_unchanged')
  def sidecar_parent_symlink(f):
   p=Path(f['member']['adapter_receipt']).parent;p.rmdir();target=f['root']/'outside_receipts';target.mkdir();p.symlink_to(target)
  reject('G1_wrapper_parent_symlink',sidecar_parent_symlink,expected='PRIOR_PRIVATE_VALIDATION_OUTPUT_PRESERVED_NO_RETRY')
  # G2: measured original receipt parse swap now fails before sidecar write.
  def parse_swap(a,f):
   loads=json.loads
   def swaps(text,*args,**kw):
    obj=loads(text,*args,**kw)
    if isinstance(obj,dict) and obj.get('pipeline_status')=='STREAM_HARMONIZE_DIRECT_MUNGE_PASS':
     Path(f['member']['new_receipt']).write_text(json.dumps({**obj,'pipeline_status':'CHANGED_AFTER_PARSE'})+'\n')
    return obj
   a.json=types.SimpleNamespace(loads=swaps,dumps=json.dumps)
  reject('G2_original_receipt_parse_swap',hook=parse_swap,expected='FROZEN_VALIDATION_INPUT_CHANGED',outputs_may_exist=True)
  # G3: old post-sidecar writes now fail despite preserved provisional sidecar.
  def post_write_swap(target):
   def hook(a,f):
    native=type(Path())
    class Writer:
     def __init__(self,h):self.h=h
     def __enter__(self):return self
     def write(self,text):
      value=self.h.write(text)
      if target=='QC':p=Path(f['member']['new_qc']);p.write_text(p.read_text().replace('input_rows\t34','input_rows\t999'))
      elif target=='PLAN':f['pp'].write_text(f['pp'].read_text()+' ')
      elif target=='BODY':f['body'].write_bytes(f['body'].read_bytes()+b'extra')
      elif target=='RUNTIME':
       p=Path(f['plan']['collector_runtime_profile']['resolved_binary']);p.chmod(0o644);p.write_text('CHANGED_BINARY\n')
      elif target=='ACQUISITION':c.save(str(f['plan']['acquisition_terminal']['seal_path'])+'.failure.json',dict(failed=True))
      return value
     def __exit__(self,*args):return self.h.__exit__(*args)
    class RoutedPath(native):
     def open(self,*args,**kw):
      h=super().open(*args,**kw)
      if str(self)==f['member']['adapter_receipt'] and args and args[0]=='x':return Writer(h)
      return h
    a.Path=RoutedPath
   return hook
  for target,expected in [('QC','FROZEN_VALIDATION_INPUT_CHANGED'),('PLAN','FIXED_PLAN_REQUIRED_AFTER_SIDECAR'),('BODY','FULL_ORIGINAL_VALIDATION_BODY_CHANGED'),('RUNTIME','FROZEN_VALIDATION_INPUT_CHANGED'),('ACQUISITION','STAGE_COMMIT_CHANGED_DURING_CONSUMPTION')]:
   reject('G3_post_write_'+target,hook=post_write_swap(target),expected=expected,sidecar_may_exist=True,outputs_may_exist=True)
  def late_sidecar(a,f):
   original=a.completed_acquisition;calls=[]
   def gate(plan):
    original(plan);calls.append(1)
    if len(calls)==3:Path(f['member']['adapter_receipt']).write_text('CHANGED_FINAL_WRAPPER\n')
   a.completed_acquisition=gate
  reject('final_sidecar_drift',hook=late_sidecar,expected='FIXED_INTENDED_SIDECAR_REQUIRED_AT_RETURN',sidecar_may_exist=True,outputs_may_exist=True)
  # Receipt-count/schema/semantic and family terminal controls are metadata only.
  reject('duplicate_26_paths',lambda f:f['plan']['acquisition_terminal'].__setitem__('individual_receipt_paths',[f['plan']['acquisition_terminal']['individual_receipt_paths'][0]]*26),expected='ALL_THIRTEEN_DISTINCT_DOUBLE_RECEIPT_PATHS_REQUIRED')
  reject('duplicate_pair_members',lambda f:f['plan']['acquisition_terminal']['source_receipt_pairs'][1].__setitem__('member',f['plan']['acquisition_terminal']['source_receipt_pairs'][0]['member']),expected='EXACT_THIRTEEN_DISTINCT_ACQUISITION_MEMBERS_REQUIRED')
  reject('missing_pair',lambda f:f['plan']['acquisition_terminal']['source_receipt_pairs'].pop(),expected='EXACT_THIRTEEN_DISTINCT_ACQUISITION_MEMBERS_REQUIRED')
  reject('same_primary_mirror',lambda f:f['plan']['acquisition_terminal']['source_receipt_pairs'][0].__setitem__('mirror_path',f['plan']['acquisition_terminal']['source_receipt_pairs'][0]['primary_path']),expected='ALL_THIRTEEN_DISTINCT_DOUBLE_RECEIPT_PATHS_REQUIRED')
  reject('metadata_missing_individual',lambda f:f['plan']['acquisition_terminal']['metadata_sha256'].pop(f['plan']['acquisition_terminal']['source_receipt_pairs'][0]['primary_path']),expected='EXACT_FROZEN_INDIVIDUAL_RECEIPT_MAP_REQUIRED')
  for name,edit in [('failed_source',lambda r:r.__setitem__('status','FAILED')),('unreaped_source',lambda r:r['teardown'].__setitem__('teardown_verified',False)),('residual_group',lambda r:r['teardown'].__setitem__('remaining_group_members',[999])),('cleanup_error',lambda r:r['teardown'].__setitem__('cleanup_error','injected')),('nonzero_source',lambda r:r.__setitem__('returncode',1)),('source_stop',lambda r:r.__setitem__('stop_reason','injected')),('wrong_source_plan',lambda r:r.__setitem__('plan_sha256','0'*64)),('wrong_frozen_source_member',lambda r:r['member'].__setitem__('expected_bytes',999))]:
   reject(name,lambda f,edit=edit:change_receipt(f,edit),expected='EXACT_COMPLETE_REAPED_INDIVIDUAL_SOURCE_REQUIRED')
  reject('source_failure_addendum',lambda f:c.save(str(f['plan']['acquisition_terminal']['source_receipt_pairs'][0]['mirror_path'])+'.failure.json',dict(failed=True)),expected='SOURCE_FAILURE_VETOES_REPLAY')
  reject('producer_PENDING',lambda f:c.write(f['plan']['acquisition_terminal']['pending_path'],'PENDING'),expected='STAGE_TERMINAL_COMMIT_NOT_COMPLETE')
  reject('missing_primary_copy',lambda f:f['plan']['acquisition_terminal']['primary_receipt_sha256'].pop(next(iter(f['plan']['acquisition_terminal']['primary_receipt_sha256']))),expected='TWO_BYTE_IDENTICAL_FULL_FAMILY_PRIMARY_COPIES_REQUIRED')
  for name,edit in [('wrong_family_status',lambda x:x.__setitem__('status','FAILED')),('wrong_family_count',lambda x:x.__setitem__('completed_source_count',12)),('wrong_family_plan',lambda x:x.__setitem__('plan_sha256','0'*64)),('wrong_family_source_map',lambda x:x.__setitem__('source_receipts',x['source_receipts'][:-1]))]:
   reject(name,lambda f,edit=edit:bind_contract(f,edit),expected='EXACT_CURRENT_FULL13_ACQUISITION_FAMILY_REQUIRED')
  reject('producer_terminal_failure',lambda f:c.save(str(f['plan']['acquisition_terminal']['seal_path'])+'.failure.json',dict(failed=True)),expected='STAGE_COMMIT_CHANGED_DURING_CONSUMPTION')
  for name,edit,expected in [('wrong_runtime_logical',lambda p:p.__setitem__('logical_executable','/synthetic/wrong'),'EXACT_DECLARED_COLLECTOR_RUNTIME_ROUTE_REQUIRED'),('wrong_runtime_symlink_literal',lambda p:p.__setitem__('declared_symlinks',[]),'DECLARED_READ_ONLY_RUNTIME_IDENTITY_CHANGED'),('wrong_runtime_binary_SHA',lambda p:p.__setitem__('resolved_binary_sha256','0'*64),'DECLARED_READ_ONLY_RUNTIME_IDENTITY_CHANGED')]:
   reject(name,lambda f,edit=edit:runtime_edit(f,edit),expected=expected)
  def omit_binary_dependency(f):f['plan']['dependencies_sha256'].pop(f['plan']['collector_runtime_profile']['resolved_binary'])
  reject('unbound_runtime_binary',omit_binary_dependency,expected='DECLARED_READ_ONLY_RUNTIME_IDENTITY_CHANGED')
  def changed_runtime_route(f):
   p=Path(f['plan']['python']);p.unlink();p.symlink_to(f['body'])
  reject('changed_runtime_symlink_target',changed_runtime_route,expected='EXACT_DECLARED_COLLECTOR_RUNTIME_ROUTE_REQUIRED')
  # Reuse original EOF fixture logic without rerunning old broad component studies.
  def crc(f):
   data=bytearray(f['body'].read_bytes());data[-8]^=1;f['body'].write_bytes(data);ident=f['identity']();f['member']['original_source_identity']=ident
   f['row'].update(replication_content_length_bytes=str(ident['bytes']),replication_checksum='md5:'+ident['md5'],replication_etag=ident['md5'])
   with Path(f['plan']['queue']).open(newline='') as h:rows=list(c.csv.DictReader(h,delimiter='\t'))
   rows[0]=f['row'];c.write_tsv(Path(f['plan']['queue']),rows);f['plan']['dependencies_sha256'][f['plan']['queue']]=c.sha(f['plan']['queue']);f['member']['collector_argv'][-1]=format(ident['bytes']/(1<<30),'.17g')
  reject('actual_CRC_failure',crc,expected='CRC check failed')
  def EOF_mismatch(a,f):
   origload=a.load
   def load(name,path):
    module=origload(name,path)
    if name=='_original_validation_streaming':
     orig=module.open_verified_gzip_text
     @contextlib.contextmanager
     def reader(**kw):
      with orig(**kw) as (text,receipt):yield text,receipt
      receipt['observed_sha256']='0'*64
     module.open_verified_gzip_text=reader
    return module
   a.load=load
  reject('outer_EOF_full_SHA',hook=EOF_mismatch,expected='FULL_LOCAL_GZIP_EOF_IDENTITY_CHANGED')
 finally:c.urllib.request.urlopen=oldurlopen
 for p,h in c.EXPECTED.items():c.assert_(c.sha(p)==h,'pinned_after_'+p)
 result=dict(scope='NARROW_SYNTHETIC_ADAPTER96_BOUNDARY',status='PASS_COMPONENT_CONTROLS_NO_OPERATIONAL_ADMISSION',source_sha256=c.EXPECTED,
  case_count=len(c.cases),check_count=len(c.checks),cases=c.cases,checks=c.checks,all_expectations_pass=all(x['passed'] for x in c.checks),
  runtime=dict(python=sys.version,executable=sys.executable),synthetic_runtime_profile='Own tiny read-only nonexecutable identity file and declared symlink; not certification of actualPython39 execution.',
  no_network=True,no_production_bodies=True,no_production_mutex=True,no_production_worker=True,no_new_estimators=True,
  threshold_amendment='Only900000/700000 fixture constants become1; production original untouched.',inherited_v1_checks=202)
 c.save(HERE/'validation_collector_adapter_controls_v2.json',result)
 print(json.dumps({k:result[k] for k in ['status','case_count','check_count']}))
if __name__=='__main__':main()
