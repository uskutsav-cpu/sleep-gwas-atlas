"""Narrow independent metadata/design review. No body hashing or execution."""
import ast
import copy
import csv
import datetime
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import struct
import subprocess
import sys
import time
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
PLAN=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/core_pipeline/large35_checkpoint_continuation_v1/core_large35_checkpoint_continuation_plan_v1.json')
PLAN_SHA='716b9df426dc5291296a4dbbedd9f057f615168c8c61945da4ccb7dd074f7bb0'
OUT=R/'independent_core_large35_checkpoint_continuation_prelaunch_v1.json'
started=time.monotonic();bindings={};checks=0

def check(condition):
 global checks
 assert condition
 checks+=1

def read(q,expected=None):
 q=Path(q);check(q.is_file() and not q.is_symlink() and q.stat().st_size<8<<20)
 h=hashlib.sha256(q.read_bytes()).hexdigest()
 if expected:check(h==expected)
 bindings[str(q)]=h
 return q.read_text()

plan=json.loads(read(PLAN,PLAN_SHA));old=json.loads(read(plan['predecessor_plan'],plan['predecessor_plan_sha256']))
author=json.loads(read(R/'core_large35_checkpoint_continuation_author_receipt_v1.json','49cb777a9757545e61309b34280a99b105b48ef6f569327cb8df76cd68beb729'))
read(R/'core_large35_checkpoint_continuation_author_v1.md','70e09ab9814cb0544cb190b45c1456167165d5036c475459dd681a9e4bc9f5a7')
for name in ['core_large35_checkpoint_adoption_v1.py','core_large35_checkpoint_continuation_v1.py','prepare_core_large35_checkpoint_continuation_v1.py','108_acquire_validation_raw_sources_v4.py']:
 q=S/name;read(q,plan['dependencies_sha256'][str(q)])
for q in [R/'core_large35_checkpoint_delta_controls_v1.py',R/'core_large35_checkpoint_delta_controls_receipt_v1.json']:
 read(q,author['artifact_sha256'][str(q)]['sha256'])
new_prefix=plan['private_namespace'];old_prefix=old['private_namespace']
def relocate(v):
 if isinstance(v,dict):return {k:relocate(x) for k,x in v.items()}
 if isinstance(v,list):return [relocate(x) for x in v]
 if isinstance(v,str) and v.startswith(old_prefix+'/'):return new_prefix+v[len(old_prefix):]
 return v
check(plan['members']==[relocate(m) for m in old['members'][8:]])
check(len(plan['members'])==27 and len(plan['adopted_members'])==8 and plan['member_count']==27)
check([m['trait_id'] for m in plan['adopted_members']]+[m['trait_id'] for m in plan['members']]==plan['original35_order'])
check(plan['original35_order']==[m['trait_id'] for m in old['members']] and len(set(plan['original35_order']))==35)
check(plan['guard']==old['guard'])
check(plan['guard']['new_output_limit_bytes']==16<<30 and plan['guard']['global_reservation_bytes']==300<<30 and plan['guard']['observed_aggregate_worker_RSS_limit_bytes']==2<<30)
check(plan['estimator_calls']==0 and plan['expected_new_worker_count']==27*5==135)
check(all(m['prefilter_command'] is None and m['original_design']['prefilter'] is None for m in plan['members']))
check(all(plan['dependencies_sha256'].get(k)==v for k,v in old['dependencies_sha256'].items()))
check(plan['deadline_anchor_utc']==old['prepared_utc'])
deadline=datetime.datetime.fromisoformat(plan['deadline_anchor_utc'])+datetime.timedelta(seconds=plan['guard']['deadline_seconds'])
check(deadline.isoformat()=='2026-10-14T00:11:56.230105+00:00')
check(plan['execution_admitted'] is False and plan['predecessor_failed_family_reclassified'] is False)

# Confirm named unchanged inherited helpers, not their old control suites.
for q in [S/'84_run_core_large35_bounded_replay_v3.py',S/'85_validate_core_large35_bounded_replay_v3.py',S/'86_cleanup_core_bounded_spool_v3.py',S/'terminal_commit_common_v2.py']:
 read(q,old['dependencies_sha256'][str(q)])
a=ast.parse((S/'84_run_core_large35_bounded_replay_v3.py').read_text());b=ast.parse((S/'core_large35_checkpoint_continuation_v1.py').read_text())
fa={n.name:ast.dump(n,include_attributes=False) for n in a.body if isinstance(n,ast.FunctionDef)}
fb={n.name:ast.dump(n,include_attributes=False) for n in b.body if isinstance(n,ast.FunctionDef)}
inherited=['module','admission_gate','runtime_gate','acquire_mutex','main']
for name in inherited:check(fa[name]==fb[name])
controller=(S/'core_large35_checkpoint_continuation_v1.py').read_text()
check("namespace_bytes(SSD/'core_pipeline',g['new_output_limit_bytes'])" in controller)
check("namespace_bytes(SSD,g['global_reservation_bytes'])" in controller)
check('started=time.monotonic()-age' in controller and 'TerminalCommit' in controller)
check('135' in controller and 'combined35_content_certificate' in controller)

# Actual snoring metadata; streams and retained output bodies are not re-read.
snoring=plan['adopted_members'][7];member=snoring['member']
proof=json.loads(read(plan['snoring_adjudication'],plan['checkpoint_metadata_sha256'][plan['snoring_adjudication']]))
qc_text=read(member['harmonization_qc'],plan['checkpoint_metadata_sha256'][member['harmonization_qc']])
columns=json.loads(read(snoring['columns_path'],plan['checkpoint_metadata_sha256'][snoring['columns_path']]))
comparison=json.loads(read(member['comparison_receipt'],plan['checkpoint_metadata_sha256'][member['comparison_receipt']]))
failed_worker=Path(plan['predecessor_execution_receipt']).parent/'receipts_v4/snoring__compare.worker.json'
worker=json.loads(read(failed_worker,plan['checkpoint_metadata_sha256'][str(failed_worker)]))
check(worker['status']=='WORKER_FAILED_PRESERVED' and worker['returncode']==1 and worker['owned_cleanup_verified'] is True)
check(worker['process_group_teardown']['remaining_group_members']==[] and not worker['process_group_teardown'].get('cleanup_error'))
metadata={};steps=[];inside=False
for line in qc_text.splitlines():
 f=line.split('\t')
 if f==['step','dropped','remaining']:inside=True
 elif inside and len(f)==3:steps.append(dict(reason=f[0],dropped=int(f[1]),remaining=int(f[2])))
 elif len(f)==2:metadata[f[0]]=f[1]
oldqc=member['original_design']['historical_QC']
expected=copy.deepcopy(oldqc['ordered_steps']);check(expected[13]['reason']=='N_eff below 50% of total (381,973.8) [CDG3]')
expected[13]['reason']='sample size below 50% of configured effective N (381,973.8) [CDG3]'
check(steps==expected and steps[13]['dropped']==0 and steps[13]['remaining']==7168629)
check(not (set(columns['columns']) & {'N','N_EFF','N_EFF_HALF','NCASE','NCONTROL'}))
N=4/(1/152302+1/256015);threshold=N*.5
check(struct.pack('d',N)==struct.pack('d',proof['semantic_decision']['effective_N']))
check(N==381973.77557143103 and threshold==190986.88778571552 and N!=408317)
# Primary historical and current predicate source inspect, no runtime/source read.
current=read(P.parent/'scripts/01_harmonize.py',plan['dependencies_sha256'][str(P.parent/'scripts/01_harmonize.py')])
historical=subprocess.run(['git','show','546a0464^:scripts/01_harmonize.py'],cwd=P.parent,capture_output=True,check=True).stdout
check(hashlib.sha256(historical).hexdigest()==proof['predicate_evidence']['historical_code_sha256'])
check('total_effective_n = effective_n(metadata_ncase, metadata_ncontrol)' in historical.decode())
check('out["N"] = total_effective_n' in historical.decode())
check('n_reference = effective_n(metadata_ncase, metadata_ncontrol)' in current and 'out["N"] = n_reference' in current)
check('N_EFF_MIN_FRACTION = 0.5' in current)
read(P.parent/'config/analysis_panel.tsv',plan['checkpoint_metadata_sha256'][str(P.parent/'config/analysis_panel.tsv')])
with (P.parent/'config/analysis_panel.tsv').open() as f:row=next(x for x in csv.DictReader(f,delimiter='\t') if x['trait_id']=='snoring')
check(row['type']=='binary' and float(row['ncase'])==152302 and float(row['ncontrol'])==256015)
legacy=[m['trait_id'] for m in plan['members'] if any(x['reason'].startswith('N_eff below ') for x in m['original_design']['historical_QC']['ordered_steps'])]
check(legacy==[])

# Three narrow independent logical controls; no whole-controller fixture.
sys.path.insert(0,str(S));spec=importlib.util.spec_from_file_location('_independent_core_adoption',S/'core_large35_checkpoint_adoption_v1.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
control=[]
value=mod.semantic_snoring('snoring',oldqc,metadata,steps,columns,comparison)
check(value==proof['semantic_decision']);control.append(dict(case='actual_scoped_constant_N_label',status='PASS'))
for name,change in [('wrong_finite_N_Z_count',lambda x:x['munged_comparison']['counts_new'].update(finite_N_Z=1179074)),('wrong_content_digest',lambda x:x['harmonized_comparison'].update(decompressed_sha256=['0'*64]*2))]:
 v=copy.deepcopy(comparison);change(v)
 try:mod.semantic_snoring('snoring',oldqc,metadata,steps,columns,v)
 except RuntimeError as e:control.append(dict(case=name,status='EXPECTED_REJECTION',error=str(e)))
 else:raise AssertionError('SCIENTIFIC_PROOF_CONTRADICTION_ACCEPTED')
check(len(control)==3)
check(all(hashlib.sha256(Path(q).read_bytes()).hexdigest()==h for q,h in bindings.items()))
check(time.monotonic()-started<30 and resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<128<<20)
result=dict(schema='independent_core_large35_checkpoint_continuation_prelaunch_v1',verdict='QUALIFIED_ACTUAL_PLAN_PRELAUNCH_PASS_REQUIRES_EXACT_ROOT_ADMISSION',
 plan_path=str(PLAN),plan_sha256=PLAN_SHA,metadata_code_sha256=bindings,assertion_count=checks,focused_controls=control,
 exact_adopted_traits=[x['trait_id'] for x in plan['adopted_members']],exact_remaining_traits=[x['trait_id'] for x in plan['members']],new_commands=135,new_fits=0,
 deadline=deadline.isoformat(),constant_effective_N=N,minimum_N=threshold,snoring_old_worker_and_family_remain_failed=True,
 source_output_reference_runtime_body_hashes_repeated=0,production_locks_or_workers_or_fits=0,old_suites_rerun=0,
 actual_plan_has_required_review_filenames=False,inherited_admission_accepts_nonempty_review_map=True,
 required_root_admission_review_files=['independent_core_large35_checkpoint_continuation_prelaunch_v1.md','independent_core_large35_checkpoint_continuation_prelaunch_v1.json','independent_core_large35_checkpoint_continuation_prelaunch_seal_v1.json'],
 inherited_AST_helpers=inherited,execution_admission=False,actual_remaining27_RSS_capacity_time_success_unobserved=True,
 elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
with OUT.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(verdict=result['verdict'],assertions=checks,controls=3,receipt_sha256=hashlib.sha256(OUT.read_bytes()).hexdigest())))
