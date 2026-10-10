"""Essential independent v2 delta only; two bounded actual postlife frames."""
import ast
import copy
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import pickle
import re
import resource
import sys
import time
import pandas as pd
import numpy as np
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';started=time.monotonic();bindings={};checks=0
PLAN=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/core_pipeline/large35_checkpoint_continuation_v2/core_large35_checkpoint_continuation_plan_v2.json')
FIXED='8a839122ef7a9a0faefb6c3971f04abd8fa755ed630a6db1efe8f725ca54c6e3'
OUT=R/'independent_core_large35_checkpoint_continuation_prelaunch_v2.json'
def check(value):
 global checks
 assert value
 checks+=1

def read(q,expected=None):
 q=Path(q);check(q.is_file() and not q.is_symlink() and q.stat().st_size<8<<20)
 b=q.read_bytes();h=hashlib.sha256(b).hexdigest()
 if expected:check(h==expected)
 bindings[str(q)]=h
 return b

plan=json.loads(read(PLAN,FIXED));old=json.loads(read(plan['preserved_failed_checkpoint_plan'],plan['preserved_failed_checkpoint_plan_sha256']))
a=json.loads(read(R/'core_large35_checkpoint_continuation_author_receipt_v2.json','df6cf0eda538d720ceda63813edbbbbba190ef89b9bf86ede5f73ed8d3109eac'))
for q,item in a['artifacts'].items():read(q,item['sha256'])
read(R/'core_large35_checkpoint_continuation_author_v2.md','2cdaa018f059db6e12918f14c03f59bfb25003af38667e301ae97d22351ed20c')
read(R/'independent_core_large35_checkpoint_continuation_prelaunch_seal_v1.json','7247c05516afda8e50f1bbf81fc47a188811b7a6be5ada172b8e86f87931dfb3')
for name in ['core_bounded_harmonizer_v3.py','core_column_reader_v4.py','86_cleanup_core_bounded_spool_v3.py','85_validate_core_large35_bounded_replay_v3.py','core_large35_checkpoint_continuation_v1.py','108_acquire_validation_raw_sources_v4.py']:
 read(S/name,plan['dependencies_sha256'][str(S/name)])
check(pd.__version__=='2.2.3' and np.__version__=='1.26.4')

def relocated(v):
 if isinstance(v,dict):return {k:relocated(x) for k,x in v.items()}
 if isinstance(v,list):return [relocated(x) for x in v]
 if isinstance(v,str) and v.startswith(old['private_namespace']+'/'):return plan['private_namespace']+v[len(old['private_namespace']):]
 return v
expected=[relocated(x) for x in old['members']]
for m in expected:
 m['harmonize_command'][2]=str(S/'core_bounded_harmonizer_v4.py')
 if m['trait_id']=='sleep_apnea':
  cp=plan['apnea_column_checkpoint'];i=m['harmonize_command'].index('--')
  m['harmonize_command'][i:i]=['--prepared-columns',cp['manifest_path'],'--expected-columns-sha256',cp['manifest_sha256']]
  m['prepared_column_checkpoint']=cp
check(expected==plan['members'])
check(plan['adopted_members']==old['adopted_members'] and plan['adopted_worker_receipt_sha256']==old['adopted_worker_receipt_sha256'])
check(plan['guard']==old['guard'] and plan['deadline_anchor_utc']==old['deadline_anchor_utc'])
check(plan['member_count']==27 and plan['expected_new_worker_count']==135 and plan['adopted_result_count']==8 and plan['estimator_calls']==0)
check(all(plan['dependencies_sha256'].get(k)==v for k,v in old['dependencies_sha256'].items()))
deadline=(datetime.datetime.fromisoformat(plan['deadline_anchor_utc'])+datetime.timedelta(seconds=plan['guard']['deadline_seconds'])).isoformat()
check(deadline=='2026-10-14T00:11:56.230105+00:00')
check(plan['guard']['new_output_limit_bytes']==16<<30 and plan['guard']['global_reservation_bytes']==300<<30 and plan['guard']['observed_aggregate_worker_RSS_limit_bytes']==2<<30)
cp=plan['apnea_column_checkpoint'];manifest=json.loads(read(cp['manifest_path'],cp['manifest_sha256']))
check(cp['manifest_sha256']=='7ba0b854a4aada2ebc5f138310f180cb485814203133980596c479eb6bc18ff4')
check(set(manifest['typed_pieces'])=={'SNP','CHR','BP','A1','A2','FRQ','BETA','SE','P'})
raw=plan['members'][0]['original_design']['raw']
check(manifest['source']==raw['resolved_path'])
check(all(manifest[k]==raw['sealed_verified_sha256'] for k in ['source_sha256_before','source_sha256_after','expected_source_sha256']))
check(manifest['expected_rows']==20170208 and manifest['physical_shape_gate']['rows']==20170208 and manifest['physical_shape_gate']['raw_gzip_EOF_and_CRC_reached'] is True)
check(not manifest['coordinate_label_parsed'] and not manifest['original_BUILD_values'])
check(not any(manifest['preinference_token_pieces'].values()) and not any(manifest['coordinate_token_pieces'].values()))
paths=set()
for name,items in manifest['typed_pieces'].items():
 check(len(items)==404)
 for i,x in enumerate(items):
  path=Path(x['path'])
  assert path.parent==Path(cp['manifest_path']).parent and path not in paths
  assert x['start']==i*50000 and x['rows']==min(50000,20170208-i*50000)
  assert x['dtype']==manifest['global_dtype'][name] and re.fullmatch('[0-9a-f]{64}',x['sha256'])
  paths.add(path)
check(len(paths)==9*404)
source=json.loads(read(cp['source_gate_receipt'],cp['source_gate_receipt_sha256']))
failed=json.loads(read(cp['failed_harmonizer_worker'],cp['failed_harmonizer_worker_sha256']))
check(source['status']=='EXACT_CORE_RAW_SOURCE_GATE_PASS' and source['raw_sha256']==raw['sealed_verified_sha256'] and source['plan_sha256']==cp['producer_plan_sha256'])
check(failed['status']=='WORKER_FAILED_PRESERVED' and failed['returncode']==1 and failed['owned_cleanup_verified'] is True and not failed['process_group_teardown']['remaining_group_members'])
check(plan['source_checkpoint_reused_as_new_worker'] is False and plan['new_source_worker_count']==27)

# Independent exact AST normalization: only the new conditional zero-count and
# checkpoint selection alter the harmonizer; all science stages are inherited.
trees=[ast.parse((S/('core_bounded_harmonizer_v'+v+'.py')).read_text()) for v in ['3','4']]
fn=[{n.name:n for n in tree.body if isinstance(n,ast.FunctionDef)} for tree in trees]
dump=lambda n:ast.dump(n,include_attributes=False)
for name in set(fn[0])-{'original_stages','harmonize','main'}:check(dump(fn[0][name])==dump(fn[1][name]))
f=copy.deepcopy(fn[1]['original_stages']);f.body=[n for n in f.body if not (isinstance(n,ast.Assign) and any(ast.unparse(t) in ['qc_drop','removed','removed.value'] for t in n.targets)) and not (isinstance(n,ast.If) and 'ORIGINAL_DROP_REDUCTION_AST_CHANGED' in ast.unparse(n))]
for n in f.body:
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='stage_qc' for t in n.targets):
  check(ast.unparse(n.value.args[-1])=='qc_drop');n.value.args[-1]=ast.Name(id='drop_node',ctx=ast.Load())
check(dump(f)==dump(fn[0]['original_stages']))
f=copy.deepcopy(fn[1]['harmonize']);f.args.args=f.args.args[:-2];f.args.defaults=[]
old_prepared=next(n for n in fn[0]['harmonize'].body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='prepared' for t in n.targets))
new_body=[]
for n in f.body:
 if isinstance(n,ast.If) and 'BOTH_COLUMN_CHECKPOINT_CONTROLS_REQUIRED' in ast.unparse(n):continue
 if isinstance(n,ast.If) and ast.unparse(n.test)=='prepared_columns is not None':new_body.append(copy.deepcopy(old_prepared))
 else:new_body.append(n)
f.body=new_body;check(dump(f)==dump(fn[0]['harmonize']))

# Actual two small completed postlife frames, not raw/typed-column census.
sys.path.insert(0,str(S))
mods=[]
for suffix in ['3','4']:
 spec=importlib.util.spec_from_file_location('_independent_harmonizer'+suffix,S/('core_bounded_harmonizer_v'+suffix+'.py'));module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);mods.append(module)
original=P.parent/'scripts/01_harmonize.py';read(original,mods[0].ORIGINAL_SHA)
module,body=mods[0].load_original(original);stages=[m.original_stages(module,body) for m in mods]
for name in ['life','required','mapping','sample','serialize']:
 check(stages[0][name].__code__.co_code==stages[1][name].__code__.co_code and stages[0][name].__code__.co_consts==stages[1][name].__code__.co_consts)
evidence=json.loads((R/'core_apnea_empty_frame_diagnosis_controls_receipt_v1.json').read_text())['actual_small_checkpoint_evidence'];controls=[]
for path,item in evidence.items():
 if not isinstance(item,dict):continue
 frame=pickle.loads(read(path,item['sha256']));check(len(frame)==item['rows'])
 if len(frame)==0:
  try:stages[0]['qc'](frame.copy(),lambda x:x.duplicated(keep='first'))
  except TypeError as e:check(str(e)=="Cannot perform reduction 'sum' with string dtype")
  else:raise AssertionError('EXPECTED_ACTUAL_EMPTY_FAILURE_MISSING')
  new,steps=stages[1]['qc'](frame.copy(),lambda x:x.duplicated(keep='first'))
  check(len(new)==0 and len(steps)==13 and all(x[1:]==(0,0) for x in steps))
  controls.append(dict(case='actual389_zero_row_full_ordinary_QC_accounting',status='PASS',rows=0,steps=13))
 else:
  before,sa=stages[0]['qc'](frame.copy(),lambda x:x.duplicated(keep='first'))
  after,sb=stages[1]['qc'](frame.copy(),lambda x:x.duplicated(keep='first'))
  pd.testing.assert_frame_equal(before,after,check_exact=True);check(sa==sb and len(frame)==42291)
  controls.append(dict(case='actual388_nonempty_exact_values_dtypes_order_QC',status='PASS',rows=42291,retained_rows=len(after)))
 del frame
check(len(controls)==2)

# Cleanup delta: external exact-checkpoint files are validated then continue,
# never registered. Inventories/unlinks remain restricted to the new spool.
oldclean=ast.parse((S/'86_cleanup_core_bounded_spool_v3.py').read_text());newclean=ast.parse((S/'cleanup_core_bounded_spool_checkpoint_v1.py').read_text())
fa={n.name:n for n in oldclean.body if isinstance(n,ast.FunctionDef)};fb={n.name:n for n in newclean.body if isinstance(n,ast.FunctionDef)}
for name in set(fa)-{'run'}:check(dump(fa[name])==dump(fb[name]))
f=copy.deepcopy(fb['run']);clean=[]
for n in f.body:
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='external' for t in n.targets):continue
 if isinstance(n,ast.If) and ast.unparse(n.test)=="'prepared_column_checkpoint' in member":continue
 if isinstance(n,ast.For) and ast.unparse(n.target)=='item' and ast.unparse(n.iter)=='pieces':
  n.body=[x for x in n.body if not (isinstance(x,ast.If) and ast.unparse(x.test)=='path in external')]
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='record' for t in n.targets):
  pos=next(i for i,k in enumerate(n.value.keys) if isinstance(k,ast.Constant) and k.value=='external_checkpoint_files_preserved');n.value.keys.pop(pos);n.value.values.pop(pos)
 clean.append(n)
f.body=clean;check(dump(f)==dump(fa['run']))
check(all(hashlib.sha256(Path(q).read_bytes()).hexdigest()==h for q,h in bindings.items()))
check(time.monotonic()-started<60 and resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<256<<20)
r=dict(schema='independent_core_checkpoint_continuation_v2_essential_prelaunch',verdict='QUALIFIED_ACTUAL_PLAN_PRELAUNCH_PASS_REQUIRES_EXACT_ROOT_ADMISSION',plan_path=str(PLAN),plan_sha256=FIXED,
 metadata_code_and_two_small_frame_sha256=bindings,assertion_count=checks,focused_actual_frame_controls=controls,runtime=dict(pandas=pd.__version__,numpy=np.__version__),
 nonempty_scientific_AST_and_nonQC_stages_unchanged=True,external_typed_pieces_not_registered_or_unlinked=True,checkpoint_layout=dict(columns=9,pieces_per_column=404,rows=20170208,typed_piece_contents_freshly_read=0),
 preserved_failed_v1_worker_and_spool=True,new_members=27,adopted=8,new_commands=135,estimator_calls=0,deadline=deadline,original_guard_unchanged=True,
 raw_reference_runtime_full_body_hashes_or_old_suites_repeated=0,production_workers_locks_or_transfers=0,new_postlife_or_QC_caches_reused=False,
 execution_admission=False,required_root_admission_review_files=['independent_core_large35_checkpoint_continuation_prelaunch_v2.md','independent_core_large35_checkpoint_continuation_prelaunch_v2.json','independent_core_large35_checkpoint_continuation_prelaunch_seal_v2.json'],
 full85_content_QC_CRC_verification_still_required=True,actual_whole_source_RSS_storage_time_success_unobserved=True,
 elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
with OUT.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
print(json.dumps(dict(verdict=r['verdict'],assertions=checks,actual_frame_controls=len(controls),receipt_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest(),elapsed=r['elapsed_seconds'],max_RSS=r['max_RSS_bytes'])))
