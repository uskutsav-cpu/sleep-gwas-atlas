"""Narrow v1-to-v2 plan/source binding review; metadata reads only."""
import ast
import hashlib
import json
from pathlib import Path

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
Q=SSD/'core_pipeline/large35_bounded_replay_v2/core_large35_bounded_replay_plan_v2.json'
PIN='3dc24cde5ec6fe4cc2a7fcaf037cf1ea0c3d8707c14c4ec9c47531931ff86626'
CHECKS=[];IDS={};UNREAD={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bind(p,wanted=None):
 p=Path(p);h=sha(p);assert wanted is None or wanted==h,(p,wanted,h);IDS[str(p)]=h;return h
def check(n,c):
 CHECKS.append({'control':n,'pass':bool(c)});assert c,n
bind(Q,PIN);new=json.loads(Q.read_text());oldq=SSD/'core_pipeline/large35_bounded_replay_v1/core_large35_bounded_replay_plan_v1.json';bind(oldq,'eb84983bef57dfb21f06a95d9f30a3871652a7b2f3cccdca30fabf9097175ed7');old=json.loads(oldq.read_text())
check('exact35_member_count_and_order',new['member_count']==old['member_count']==len(new['members'])==35 and [m['trait_id'] for m in new['members']]==[m['trait_id'] for m in old['members']])
oldroot=old['private_namespace'];newroot=new['private_namespace']
def relocate(v):
 if isinstance(v,str):return newroot+v[len(oldroot):] if v.startswith(oldroot) else v
 if isinstance(v,list):return [relocate(x) for x in v]
 if isinstance(v,dict):return {k:relocate(x) for k,x in v.items()}
 return v
for a,b in zip(old['members'],new['members']):
 check(b['trait_id']+'_entire_member_only_private_relocation',relocate(a)==b)
 check(b['trait_id']+'_entire_original_scientific_design_unchanged',a['original_design']==b['original_design'])
 check(b['trait_id']+'_full_literal_original_and_bounded_and_munge_tokens_only_relocated',all(relocate(a[k])==b[k] for k in ['literal_original_harmonize_command','harmonize_command','munge_command']))
 check(b['trait_id']+'_exact_source_option_no_extra_or_abbreviation',b['literal_original_harmonize_command'].count('--infile')==1 and not any(x.startswith('--infil') and x!='--infile' for x in b['literal_original_harmonize_command']))
stable=['baseline_support_package','baseline_execution_plan_sha256','baseline_dependency_sha256','baseline_relocated_dependency_sha256','baseline_input_sha256','original_190_jobs','core_precision_adjudication','python','ldsc_dir','reference_prefix','environment','guard','runtime_qualification','qualified_runtime_receipt','qualified_runtime_receipt_sha256','harmonization_python','component_review_seal_sha256','bounded_adapter_review_seal_sha256','scientific_adapter_qualification','global_resource_ledger_sha256','estimator_calls','scientific_membership_or_threshold_changes','automatic_retry']
check('all_guard_baseline_runtime_scientific_qualification_keys_unchanged',all(old[k]==new[k] for k in stable))
check('original190_member_bound_maps_not_new_re_adjudication',len(new['original_190_jobs'])==190)
for p,h in new['dependencies_sha256'].items():
 f=Path(p)
 if f.name in ['hm3_grch37_variant_map.tsv.gz','hg38ToHg19.over.chain.gz','w_hm3.snplist']:
  UNREAD[p]={'expected_sha256':h,'bytes':f.stat().st_size,'body_read':False};continue
 bind(p,h)
check('211_dependencies_208_metadata_3_stat_only_reference_bodies',len(new['dependencies_sha256'])==211 and len(UNREAD)==3)
for name,pin in [('core_large35_bounded_controller_review_seal_v1.json','06ba91909321a1e7d0f481bf565f5f0ad938dabd139995f803d1fe7ebdf098e1'),('independent_spool_disposal_seal_v1.json','642fa6e1fbd12d8e36bae5d7172db184ec3fed0851cfae9c31a6103a8e8e7217')]:
 p=R/name;bind(p,pin);seal=json.loads(p.read_text())
 check(name+'_is_bound_by_successor_plan',new['dependencies_sha256'].get(str(p))==pin)
 for f,value in seal.get('artifacts',seal.get('file_sha256',{})).items():
  h=value['sha256'] if isinstance(value,dict) else value;bind(f,h);check(Path(f).name+'_preserved_rejection_artifact_bound',new['dependencies_sha256'].get(f)==h)
def functions(file):return {x.name:ast.dump(x,include_attributes=False) for x in ast.parse(file.read_text()).body if isinstance(x,ast.FunctionDef)}
f74=functions(S/'74_run_core_large35_bounded_replay_v1.py');f78=functions(S/'78_run_core_large35_bounded_replay_v2.py');f75=functions(S/'75_validate_core_large35_bounded_replay_v1.py');f79=functions(S/'79_validate_core_large35_bounded_replay_v2.py')
for n in ['module','admission_gate','runtime_gate','namespace_bytes','acquire_mutex','main']:check(n+'_inherited_controller_AST_unchanged',f74[n]==f78[n])
for n in ['compare_harmonized','qc']:check(n+'_scientific_comparison_AST_unchanged',f75[n]==f79[n])
t=(S/'78_run_core_large35_bounded_replay_v2.py').read_text()
check('immediate_harmonize_capture_before_next_worker',t.index("if label=='harmonize':")>t.index('master[\'worker_receipt_sha256\'][str(record_path)]=sha(record_path)') and t.index("master['generated_metadata_sha256'][trait]={str(f):sha(f) for f in generated}")<t.index("active[0]=None;admission_gate(plan,a)"))
check('new_generated_files_enforce_leaf_and_parent_no_symlink',"if any(not f.is_file() or f.is_symlink() or any(parent.is_symlink() for parent in f.parents) for f in generated)" in t)
for flag in ['--expected-candidate-sha256','--expected-columns-sha256','--expected-comparison-worker-sha256']:check('dynamic_cleanup_argument_'+flag,t.count(flag)==1)
check('comparison_SHA_argument_uses_immediately_frozen_worker_map',"'--expected-comparison-worker-sha256',master['worker_receipt_sha256'][str(comparison_worker)]" in t)
check('retained_copy_final_binding_uses_original_generated_map',"if item['output_sha256'][member['bounded_receipt']]!=generated[" in t and "item['output_sha256'][member['column_manifest']]!=generated[" in t)
check('175_command_count_and_postpersist_terminal_full_gates_retained',"len(master['worker_receipt_sha256'])!=175" in t and "identity_gate=final_identity" in t and "'POST_PERSISTENCE_CORE_TERMINAL'" in t)
edge=R/'core_large35_validator_identity_edge_receipt_v2.json';bind(edge);e=json.loads(edge.read_text());check('actual_validator79_false_PASS_measured',any(x['is_false_PASS_witness'] for x in e['cases']) and e['source_sha256']==bind(S/'79_validate_core_large35_bounded_replay_v2.py'))
check('all_plan_source_code_and_review_hashes_unchanged',IDS=={p:sha(p) for p in IDS})
receipt={'status':'V2_REJECT_MEASURED_QC_CONSUMPTION_BINDING_ESCAPE','plan_path':str(Q),'plan_sha256':PIN,'checks':CHECKS,'check_count':len(CHECKS),'before':IDS,'after':{p:sha(p) for p in IDS},'unread_reference_bodies':UNREAD,'actual_validator_fixture_count':len(e['cases']),'inherited22_v1_controller_scenarios_not_rerun':True,'175_command_inventory':'35 unchanged source/harmonize/munge/compare/cleanup chains;3 dynamic disposal digest arguments source-bound statically. No v2 controller simulation claimed.','real_mutex_operations':0,'actual_workers':0,'production_body_reads':0,'full_runtime_tree_rehashes':0,'fits':0,'execution_admitted':False}
with (R/'core_large35_bounded_successor_binding_receipt_v2.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({'status':receipt['status'],'checks':len(CHECKS),'metadata_code_identities':len(IDS),'actual_validator_cases':len(e['cases'])}))
