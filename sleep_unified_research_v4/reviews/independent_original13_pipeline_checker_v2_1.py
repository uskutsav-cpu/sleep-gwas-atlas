#!/usr/bin/env python3
"""Independent, synthetic-only candidate consumption/commit controls. No subprocess."""
import ast
import contextlib
import copy
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace
from unittest.mock import patch

P=Path(__file__).resolve().parents[1]
S=P/'scripts'
ROOT=P.parents[1]/'independent_original13_pipeline_fixture_v2_1'
RECEIPT=P/'reviews/independent_original13_pipeline_controls_v2_1.json'
sys.path.insert(0,str(S))
EXPECTED={'97_validate_original_validation_replay_v2.py':'d80ec493c4221a7309be3467463c3d47dbcd66ae0159c0af78ea48e2cefd33a0',
 '98_prepare_original_validation_pipeline_v2.py':'eb80c1e2185490c921bb1e87d0e5ca45fe4cc788664e9042e90caec845b87a7c',
 '99_run_original_validation_pipeline_v2.py':'c47dccbaf3e353c28a5cec4b51e9a191b4d828c93d29c3c4837cb0f327eee416'}
START=time.monotonic();CHECKS=[];MAP={}
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def save(p,v):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x') as f:f.write(json.dumps(v,indent=2,allow_nan=False)+'\n')
def rewrite(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def load(n,p):
 spec=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
V=load('_independent_original13_validator',S/'97_validate_original_validation_replay_v2.py')
A=load('_independent_original13_adapter',S/'96_replay_original_validation_collector_v2.py')
C=load('_independent_original13_controller',S/'99_run_original_validation_pipeline_v2.py')
T=load('_independent_original13_terminal',S/'terminal_commit_common_v2.py')
def check(name,value,detail=None):
 CHECKS.append({'name':name,'pass':bool(value),'detail':detail})
 if not value:raise AssertionError(name+': '+str(detail))
def expect_fail(name,fn):
 try:fn()
 except BaseException as e:check(name,True,type(e).__name__+': '+str(e));return
 check(name,False,'accepted')
def bounded():
 if time.monotonic()-START>180:raise RuntimeError('OWN_REVIEW_180_SECONDS_LIMIT')
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 if sys.platform!='darwin':rss*=1024
 if rss>128<<20:raise RuntimeError('OWN_REVIEW_128MIB_RSS_LIMIT')
 size=sum(p.stat().st_size for p in ROOT.rglob('*') if p.is_file() and not p.is_symlink())
 if size>32<<20:raise RuntimeError('OWN_REVIEW_32MIB_OUTPUT_LIMIT')
 return rss,size
def fixture(name):
 d=ROOT/name;d.mkdir(parents=True);ssd=d/'ssd';out=ssd/'validation_pipeline_replay_v2'
 for x in ['workspace','receipts_v4','logs_v4']:(out/x).mkdir(parents=True)
 acquisition=d/'acquisition';acquisition.mkdir();pairs=[];metadata={};sources=[]
 for i in range(1,14):
  member={'source_id':'fixture'+str(i),'index':i,'expected_sha256':str(i)*64}
  r={'member':member,'status':'EXACT_ORIGINAL_SOURCE_BODY_ACQUIRED','plan_sha256':'a'*64,'returncode':0,'stop_reason':None,'teardown':{'teardown_verified':True,'remaining_group_members':[]}}
  p=acquisition/(str(i)+'.primary.json');q=acquisition/(str(i)+'.mirror.json');save(p,r);save(q,r)
  digest=sha(p);metadata.update({str(p):digest,str(q):digest});sources.append({'path':str(p),'sha256':digest})
  pairs.append({'member':member,'primary_path':str(p),'mirror_path':str(q),'sha256':digest})
 family={'status':'ALL13_EXACT_ORIGINAL_SOURCE_BODIES_ACQUIRED','completed_source_count':13,'plan_sha256':'a'*64,'source_receipts':sources}
 primary={}
 for name in ['family.primary.json','family.mirror.json']:
  p=acquisition/name;save(p,family);primary[str(p)]=sha(p);metadata[str(p)]=sha(p)
 binding={'plan_sha256':'a'*64,'admission_sha256':'b'*64,'executor_sha256':'c'*64}
 pending=acquisition/'PENDING.json';seal=acquisition/'terminal.json'
 terminal=T.TerminalCommit(pending,seal,binding)
 assert terminal.commit(primary,identity_gate=lambda:None,resource_gate=lambda:None,termination_gate=lambda:None)
 metadata[str(seal)]=sha(seal)
 contract={'metadata_sha256':metadata,'source_receipt_pairs':pairs,'individual_receipt_paths':[p[k] for p in pairs for k in ['primary_path','mirror_path']],
  'primary_receipt_sha256':primary,'binding':binding,'pending_path':str(pending),'seal_path':str(seal)}
 dependencies={str(S/k):v for k,v in EXPECTED.items()}
 for n in ['queue','reference','alleles','collector','streaming']:
  p=d/(n+'.txt');p.write_text('OWN SYNTHETIC '+n+'\n');dependencies[str(p)]=sha(p)
 python=Path('/usr/bin/python3');binary=python.resolve(strict=True);dependencies[str(binary)]=sha(binary)
 profile={'logical_executable':str(python),'resolved_binary':str(binary),'resolved_binary_sha256':sha(binary),
  'declared_symlinks':[{'path':str(v),'target':str(v.readlink())} for v in [python,*python.parents] if v.is_symlink()]}
 members=[]
 for i in range(1,14):
  sid='fixture'+str(i);raw=d/(sid+'.raw.gz');raw.write_bytes(gzip.compress(b'OWN INVENTED SOURCE\n',mtime=0))
  original=d/(sid+'.original.gz');original.write_bytes(gzip.compress(b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t0.125\t123.123456789\nrs2\tG\tT\t-2.375\t123.123456789\n',mtime=0))
  oq=d/(sid+'.original.qc');oq.write_text('metric\tvalue\tnotes\nmunged_output_sha256\t'+sha(original)+'\tfrozen\neffective_sample_size\t123.123456789\tfrozen\noutput_rows\t2\tfrozen\n')
  orc=d/(sid+'.original.json');save(orc,{'OWN_SYNTHETIC':True})
  row={'replication_munged_path':'discovery_extension/results/replication/munged/'+sid+'.sumstats.gz','replication_receipt_path':'discovery_extension/results/replication/receipts/'+sid+'.json','replication_source_url':'https://invalid.example/'+sid}
  identity={'bytes':raw.stat().st_size,'md5':hashlib.md5(raw.read_bytes()).hexdigest(),'sha256':sha(raw)}
  m={'source_id':sid,'index':i,'body_path':str(raw),'original_source_identity':identity,'frozen_source_metadata':row,
   'new_munged':str(out/'workspace'/row['replication_munged_path']),'new_receipt':str(out/'workspace'/row['replication_receipt_path']),
   'new_qc':str(out/'workspace'/('discovery_extension/results/replication/qc/'+sid+'.tsv')),
   'adapter_receipt':str(out/'receipts_v4'/(sid+'.adapter.json')),'source_gate_receipt':str(out/'receipts_v4'/(sid+'.source.json')),
   'comparison_receipt':str(out/'receipts_v4'/(sid+'.comparison.json')),'original_munged':str(original),'original_munged_sha256':sha(original),
   'original_qc':str(oq),'original_qc_sha256':sha(oq),'original_receipt':str(orc),'original_receipt_sha256':sha(orc),
   'expected_output_rows':2,'effective_N_serialized_12g':'123.123456789'}
  members.append(m)
 names=['independent_spool_disposal_review_v3.md','independent_spool_disposal_receipt_v3_1.json','independent_spool_disposal_seal_v3.json']
 plan={'member_count':13,'members':members,'command_count':39,'new_estimator_fits':0,'new_network_transfers':0,'jobs':[],
  'adapter_path':str(S/'96_replay_original_validation_collector_v2.py'),'output_namespace':str(out),'workspace':str(out/'workspace'),
  'pending_path':str(out/'validation_pipeline_pending_v2.json'),'terminal_seal_path':str(out/'validation_pipeline_terminal_v2.json'),
  'primary_receipt_path':str(out/'validation_pipeline_execution_receipt_v2.json'),'python':str(python),'collector_runtime_profile':profile,
  'dependencies_sha256':dependencies,'acquisition_terminal':contract,'required_review_filenames':names,
  'queue':str(d/'queue.txt'),'reference':str(d/'reference.txt'),'alleles':str(d/'alleles.txt'),
  'original_collector':str(d/'collector.txt'),'original_streaming_io':str(d/'streaming.txt'),
  'guard':{'worker_count':1,'BLAS_threads':1,'internal_floor_bytes':3<<30,'SSD_floor_bytes':5<<30,'observed_aggregate_worker_RSS_limit_bytes':2<<30,
   'new_output_limit_bytes':2<<30,'global_reservation_bytes':300<<30,'deadline_seconds':96*3600,'per_worker_deadline_seconds':7200,
   'shared_heavy_worker_lock':str(ssd/'native_heavy_worker.lock'),'checkpoint_poll_seconds':2}}
 pp=d/'plan.json';save(pp,plan)
 admission=d/'admission.json';save(admission,{'execution_admitted':True,'plan_sha256':sha(pp),'executor_sha256':sha(S/'99_run_original_validation_pipeline_v2.py'),
  'independent_review_sha256':{str(P/'reviews'/n):sha(P/'reviews'/n) for n in names}})
 args=SimpleNamespace(plan=pp,plan_sha=sha(pp),admission=admission,admission_sha=sha(admission))
 return d,ssd,out,plan,args
def collector(plan,a,m):
 for k in ['new_munged','new_qc','new_receipt']:Path(m[k]).parent.mkdir(parents=True,exist_ok=True)
 # Same ordered content; different gzip metadata remains explicitly qualified.
 Path(m['new_munged']).write_bytes(gzip.compress(gzip.decompress(Path(m['original_munged']).read_bytes()),mtime=1))
 Path(m['new_qc']).write_text(Path(m['original_qc']).read_text().replace(m['original_munged_sha256'],sha(m['new_munged'])))
 ident=m['original_source_identity'];r={'pipeline_status':'STREAM_HARMONIZE_DIRECT_MUNGE_PASS','replication_source_id':m['source_id'],'source_url':m['frozen_source_metadata']['replication_source_url'],
  'source_verification':{'verification_status':'PASS','observed_sha256':ident['sha256'],'observed_md5':ident['md5'],'observed_size_bytes':ident['bytes'],'source':m['body_path']},
  'munged_output_sha256':sha(m['new_munged']),'harmonization_qc_sha256':sha(m['new_qc']),
  'queue_sha256':plan['dependencies_sha256'][plan['queue']],'reference_sha256':plan['dependencies_sha256'][plan['reference']],
  'hm3_alleles_sha256':plan['dependencies_sha256'][plan['alleles']],
  'pipeline_code_sha256':{p:plan['dependencies_sha256'][p] for p in [plan['original_collector'],plan['original_streaming_io']]}}
 save(m['new_receipt'],r)
 save(m['adapter_receipt'],{'plan_sha256':a.plan_sha,'source_id':m['source_id'],'original_source_identity':ident,
  'output_sha256':{m[k]:sha(m[k]) for k in ['new_munged','new_receipt','new_qc']},'full_local_source_retained':True,
  'full_source_gzip_CRC_EOF_MD5_size_SHA256_verified':True,'source_or_independent_replication_admitted':False})
def validator_fixture(name):
 d,ssd,out,plan,a=fixture(name);m=plan['members'][0]
 V.run(a.plan,a.plan_sha,m['source_id'],'source',Path(m['source_gate_receipt']));collector(plan,a,m)
 return d,ssd,out,plan,a,m
def validator_controls():
 d,ssd,out,plan,a,m=validator_fixture('validator_valid')
 V.run(a.plan,a.plan_sha,m['source_id'],'compare',Path(m['comparison_receipt']))
 r=json.loads(Path(m['comparison_receipt']).read_text());check('validator_literal_and_CRC_positive',r['stream_comparison']['data_rows']==2 and r['stream_comparison']['both_full_gzip_CRC_and_EOF_verified'])
 check('gzip_byte_difference_qualified',not r['compressed_byte_identity'] and r['compressed_metadata_difference_cause_not_inferred'])
 for label,key,hook in [('qc_changed_after_parse','new_qc','qc'),('munged_changed_after_stream','new_munged','stream'),('collector_changed_after_parse','new_receipt','qc'),('sidecar_changed_after_capture','adapter_receipt','qc')]:
  d,ssd,out,plan,a,m=validator_fixture(label);fired=[False]
  original=V.qc_rows if hook=='qc' else V.literal_streams
  def altered(*args):
   result=original(*args)
   if not fired[0] and (hook=='stream' or str(args[0])==m['new_qc']):
    fired[0]=True;Path(m[key]).write_bytes(Path(m[key]).read_bytes()+b'\nMUTATION')
   return result
  with patch.object(V,'qc_rows' if hook=='qc' else 'literal_streams',altered):
   expect_fail(label,lambda:V.run(a.plan,a.plan_sha,m['source_id'],'compare',Path(m['comparison_receipt'])))
 for label,data in [('wrong_Z',b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t0.126\t123.123456789\n'),('duplicate_SNP',b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1\t1\nrs1\tG\tT\t1\t1\n'),('ambiguous_alleles',b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tT\t1\t1\n'),('nonfinite_N',b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1\tnan\n')]:
  d=ROOT/label;d.mkdir();p=d/'a.gz';q=d/'b.gz';p.write_bytes(gzip.compress(data));q.write_bytes(p.read_bytes() if label!='wrong_Z' else gzip.compress(data.replace(b'0.126',b'0.125')))
  expect_fail(label,lambda:V.literal_streams(p,q,2 if label=='duplicate_SNP' else 1))
 d=ROOT/'CRC_fault';d.mkdir();p=d/'bad.gz';q=d/'old.gz';q.write_bytes(gzip.compress(b'SNP\tA1\tA2\tZ\tN\nrs1\tA\tC\t1\t1\n'));b=bytearray(q.read_bytes());b[-8]^=1;p.write_bytes(b)
 expect_fail('full_gzip_CRC_fault',lambda:V.literal_streams(p,q,1))
def acquisition_controls():
 d,ssd,out,plan,a=fixture('acquisition_valid');A.completed_acquisition(plan);A.runtime_gate(plan);check('two_family_and26_source_copies_positive',True)
 variants=['pending','source_failure','seal_failure','source_copy_changed','family_copy_changed','duplicate_path','cleanup_error','missing_source_copy','changed_runtime_link']
 for name in variants:
  d,ssd,out,plan,a=fixture('acquisition_'+name);c=plan['acquisition_terminal'];pair=c['source_receipt_pairs'][0]
  if name=='pending':save(c['pending_path'],{'status':'PENDING'})
  elif name=='source_failure':save(pair['primary_path']+'.failure.json',{'failed':True})
  elif name=='seal_failure':save(c['seal_path']+'.failure.json',{'failed':True})
  elif name=='source_copy_changed':Path(pair['mirror_path']).write_text('{}')
  elif name=='family_copy_changed':Path(list(c['primary_receipt_sha256'])[1]).write_text('{}')
  elif name=='duplicate_path':pair['mirror_path']=pair['primary_path'];c['individual_receipt_paths']=[p[k] for p in c['source_receipt_pairs'] for k in ['primary_path','mirror_path']]
  elif name=='missing_source_copy':Path(pair['mirror_path']).unlink()
  elif name=='changed_runtime_link':plan['collector_runtime_profile']['declared_symlinks'].append({'path':plan['python'],'target':'untrue'})
  elif name=='cleanup_error':
   r=json.loads(Path(pair['primary_path']).read_text());r['teardown']['cleanup_error']='still failed'
   for p in [pair['primary_path'],pair['mirror_path']]:rewrite(p,r);c['metadata_sha256'][p]=sha(p)
   pair['sha256']=sha(pair['primary_path'])
  expect_fail('acquisition_'+name,lambda:A.runtime_gate(plan) if name=='changed_runtime_link' else A.completed_acquisition(plan))
def controller_case(name,fault=None):
 d,ssd,out,plan,a=fixture('controller_'+name);calls=[];events=[];terminations=[];held=set();fired=[False]
 original_open=os.open;original_close=os.close;original_write=C.write_new
 def own_open(path,*args,**kw):
  fd=original_open(path,*args,**kw)
  if str(path)==plan['guard']['shared_heavy_worker_lock']:held.add(fd);events.append('mock_lock_fd_open')
  return fd
 def own_close(fd):
  if fd in held:held.remove(fd);events.append('mock_lock_fd_close')
  return original_close(fd)
 def stage(_p,start,where):
  state,reason=protected.limits(_p,start,0)
  if reason:raise RuntimeError(reason)
  if fault=='post_persistence_resource' and where=='POST_PERSISTENCE_VALIDATION_TERMINAL':raise RuntimeError('SYNTHETIC_TERMINAL_RESOURCE_FAILURE')
 def termination(where):
  if terminations:raise RuntimeError('SYNTHETIC_DEFERRED_TERMINATION')
 def cleanup(mon,ownership,plan_sha):
  events.append('cleanup_with_fd' if held else 'cleanup_no_fd');ownership[:]=[True,[]]
 def fakeworker(command,prefix,worker,_p,_path,_sha,start,monitor,ownership):
  calls.append(command);check('mock_science_runs_only_under_inherited_fd_'+name+'_'+str(len(calls)),bool(held))
  sid=worker.name.split('__')[0];m=next(x for x in plan['members'] if x['source_id']==sid);label=worker.name.split('__')[1].split('.')[0]
  if label=='collector':collector(plan,a,m)
  else:V.run(a.plan,a.plan_sha,sid,'source' if label=='source' else 'compare',Path(m['source_gate_receipt' if label=='source' else 'comparison_receipt']))
  if fault=='comparison_bad_consumed' and label=='compare' and not fired[0]:
   r=json.loads(Path(m['comparison_receipt']).read_text());r['consumed_output_sha256'][m['new_qc']]='0'*64;rewrite(m['comparison_receipt'],r);fired[0]=True
  stdout=out/'logs_v4'/(worker.stem+'.stdout.log');stdout.write_text('OWN MOCK: no native worker executed\n')
  journal=Path(str(worker)+'.failure_journal.jsonl');journal.write_text('{"event":"OWN_MOCK_PRELAUNCH"}\n{"event":"OWN_MOCK_REAPED"}\n')
  paths=[stdout,journal]+[p for p in prefix.parent.glob(prefix.name+'*') if p.is_file()]
  r={'status':'WORKER_COMPLETE_VERIFIED','command':command,'plan_sha256':a.plan_sha,'returncode':0,'stop_reason':None,'owned_cleanup_verified':True,
   'process_group_teardown':{'remaining_group_members':[]},'output_sha256':{str(p):sha(p) for p in paths},'failure_journal':str(journal),
   'metadata_errors':[],'plan_unchanged':True}
  if len(calls)==1:
   if fault=='metadata_error':r['metadata_errors']=['FAILED_METADATA_DURABILITY']
   elif fault=='plan_unchanged_false':r['plan_unchanged']=False
   elif fault=='teardown_cleanup_error':r['process_group_teardown']['cleanup_error']='exception'
   elif fault=='wrong_command':r['command']=command+['UNREGISTERED']
   elif fault=='missing_stdout':r['output_sha256'].pop(str(stdout))
   elif fault=='missing_journal':r['output_sha256'].pop(str(journal))
   elif fault=='failure_addendum':save(str(worker)+'.failure.json',{'failed':True})
  save(worker,r)
 protected=SimpleNamespace(TERMINATION_REQUEST=terminations,catchable_termination=lambda signum,frame:terminations.append(signum),
  stage_resource_gate=stage,check_termination=termination,await_owned_cleanup=cleanup,baseline_gate=lambda _p:{'OWN_MOCK_190_METADATA':'x'},safe_print=lambda *a,**k:None,worker=fakeworker)
 def dispatch(modname,path):
  if 'owned_worker' in modname:return protected
  if 'native_monitor' in modname:return SimpleNamespace()
  if 'local_adapter' in modname:return A
  if 'full_compare' in modname:return V
  raise RuntimeError('UNREGISTERED_MOCK_MODULE')
 def persist(path,value):
  original_write(path,value)
  if str(path)==plan['primary_receipt_path'] and fault in ['master_mutation','plan_mutation','output_mutation','signal_after_master']:
   if fault=='master_mutation':Path(path).write_text('{}\n')
   elif fault=='plan_mutation':a.plan.write_text('{}\n')
   elif fault=='output_mutation':Path(plan['members'][0]['new_qc']).write_text('MUTATED\n')
   elif fault=='signal_after_master':terminations.append(15)
 err=None
 with patch.object(C,'SSD',ssd),patch.object(C,'OUT',out),patch.object(C,'module',dispatch),patch.object(C,'physical_mount',lambda:None),patch.object(C,'write_new',persist),patch.object(C.shutil,'disk_usage',lambda p:SimpleNamespace(free=1<<45)),patch.object(C.os,'open',own_open),patch.object(C.os,'close',own_close),patch.object(C.fcntl,'flock',lambda fd,options:events.append('mock_flock_no_real_mutex')):
  try:C.run(a)
  except BaseException as e:err=type(e).__name__+': '+str(e)
 committed=not Path(plan['pending_path']).exists() and Path(plan['terminal_seal_path']).exists()
 if committed:
  T.require_committed(plan['pending_path'],plan['terminal_seal_path'],{'plan_sha256':a.plan_sha,'admission_sha256':a.admission_sha,'executor_sha256':sha(S/'99_run_original_validation_pipeline_v2.py')},{plan['primary_receipt_path']:sha(plan['primary_receipt_path'])})
 check('all_mock_lock_fds_closed_'+name,not held)
 check('cleanup_precedes_every_mock_fd_release_'+name,all(events[i-1]=='cleanup_with_fd' for i,v in enumerate(events) if v=='mock_lock_fd_close'))
 MAP[name]={'committed':committed,'error':err,'mock_worker_count':len(calls),'directory':str(d),'master_sha256':sha(plan['primary_receipt_path']) if Path(plan['primary_receipt_path']).exists() else None,'events':events}
 bounded();return committed,err
def controller_controls():
 committed,err=controller_case('valid');check('whole13_39_mock_controller_commits',committed and err is None and MAP['valid']['mock_worker_count']==39)
 for n in ['wrong_command','missing_stdout','missing_journal','failure_addendum','comparison_bad_consumed','master_mutation','plan_mutation','output_mutation','post_persistence_resource','signal_after_master']:
  committed,err=controller_case(n,n);check('controller_rejects_'+n,not committed and err is not None)
 # Deliberately contradictory worker records: scope is receipt-oracle robustness,
 # not an assertion that the separately bound protected producer emits them.
 for n in ['metadata_error','plan_unchanged_false','teardown_cleanup_error']:
  committed,err=controller_case(n,n);check('measured_contradictory_record_acceptance_'+n,committed and err is None,MAP[n])
def guard_controls():
 tree=ast.parse((S/'99_run_original_validation_pipeline_v2.py').read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run')
 limits=next(n for n in run.body if isinstance(n,ast.FunctionDef) and n.name=='limits')
 d,ssd,out,plan,a=fixture('limits');older=ssd/'validation_pipeline_replay_v1';older.mkdir();other=ssd/'other';other.mkdir()
 env={'plan':plan,'active':[None],'shutil':SimpleNamespace(disk_usage=lambda p:SimpleNamespace(free=1<<45)),
  'SSD':ssd,'namespace_bytes':lambda p:0,'time':SimpleNamespace(monotonic=lambda:10000.)}
 exec(compile(ast.Module(body=[limits],type_ignores=[]),'<exact candidate limits AST>','exec'),env)
 state,reason=env['limits'](plan,10000.,0);check('guard_baseline_positive',reason is None)
 guards=[('rss',None,0,(2<<30)+1,'2GIB_OWNED_VALIDATION_RSS_GUARD'),('all_epochs',lambda p:(1<<30)+1 if p.name.startswith('validation_pipeline_replay_v') else 0,0,0,'2GIB_ALL_VALIDATION_PIPELINE_GUARD'),('global',lambda p:(300<<30)+1 if p==ssd else 0,0,0,'300GIB_GLOBAL_NAMESPACE_GUARD'),('stage',None,-400000,0,'96H_FAMILY_DEADLINE'),('worker',None,0,0,'2H_WORKER_DEADLINE')]
 for name,sizes,start,rss,expected in guards:
  env['namespace_bytes']=sizes or (lambda p:0);env['active'][0]=0 if name=='worker' else None
  state,reason=env['limits'](plan,start if name=='stage' else 10000.,rss);check('guard_'+name,reason==expected,state)
 for name,path,floor in [('internal','/System/Volumes/Data',3<<30),('SSD',ssd,5<<30)]:
  env['active'][0]=None;env['namespace_bytes']=lambda p:0
  env['shutil'].disk_usage=lambda p:SimpleNamespace(free=floor-1 if p==path else 1<<45)
  state,reason=env['limits'](plan,10000.,0);check('guard_'+name+'_floor',reason==('INTERNAL_SPACE_GUARD' if name=='internal' else 'SSD_SPACE_GUARD'))
 seen=[]
 with patch.object(subprocess,'Popen',lambda *a,**kw:seen.append(kw) or 'OWN_SPY_NO_PROCESS'):
  proxy=C.InheritedMutexSubprocess(71);check('descriptor_proxy_no_process',proxy.Popen(['OWN_SPY'])=='OWN_SPY_NO_PROCESS' and seen[-1]['pass_fds']==(71,))
  expect_fail('descriptor_override_rejected',lambda:proxy.Popen(['OWN_SPY'],pass_fds=(72,)))
def main():
 if ROOT.exists() or RECEIPT.exists():raise RuntimeError('OWN_IMMUTABLE_REVIEW_TARGET_ALREADY_EXISTS')
 ROOT.mkdir()
 before={str(S/k):sha(S/k) for k in EXPECTED}
 check('exact_three_reviewed_candidate_identities',before=={str(S/k):v for k,v in EXPECTED.items()})
 validator_controls();acquisition_controls();controller_controls();guard_controls()
 rss,size=bounded();check('candidate_files_unchanged',before=={str(S/k):sha(S/k) for k in EXPECTED})
 receipt={'schema':'independent_original13_pipeline_synthetic_controls_v2_1','status':'CONTROLS_COMPLETED_WITH_THREE_CONTRADICTORY_RECORD_ACCEPTANCES',
  'candidate_sha256':before,'checks':CHECKS,'check_count':len(CHECKS),'mock_controller_cases':MAP,
  'elapsed_seconds':time.monotonic()-START,'peak_RSS_bytes':rss,'synthetic_fixture_bytes':size,
  'real_production_body_reads':0,'network_calls':0,'native_worker_launches':0,'real_heavy_mutex_operations':0,
  'qualification':'Actual adapter scientific execution, acquisition13 completion, actual preparation dependencies and plan/admission remain separate gates.'}
 save(RECEIPT,receipt)
 print(json.dumps({'receipt':str(RECEIPT),'sha256':sha(RECEIPT),'checks':len(CHECKS),'elapsed_seconds':receipt['elapsed_seconds'],'peak_RSS_bytes':rss}))
if __name__=='__main__':main()
