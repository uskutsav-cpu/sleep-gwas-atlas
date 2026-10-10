"""Narrow frozen-v3 path, proof-contract and predecessor-seal bindings."""
import ast
import hashlib
import json
from pathlib import Path

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
Q=SSD/'core_pipeline/large35_bounded_replay_v3/core_large35_bounded_replay_plan_v3.json'
PIN='3e8c3c7bca26b00ed337c816a723d8a42475473667d4a9c969348045fcc9a2f4'
IDS={};CHECKS=[];UNREAD={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bind(p,wanted=None):
 p=Path(p);h=sha(p);assert wanted is None or h==wanted,(p,wanted,h);IDS[str(p)]=h;return h
def check(n,c):CHECKS.append({'control':n,'pass':bool(c)});assert c,n
bind(Q,PIN);new=json.loads(Q.read_text());oldq=SSD/'core_pipeline/large35_bounded_replay_v2/core_large35_bounded_replay_plan_v2.json';bind(oldq,'3dc24cde5ec6fe4cc2a7fcaf037cf1ea0c3d8707c14c4ec9c47531931ff86626');old=json.loads(oldq.read_text())
def relocate(v):
 if isinstance(v,str):return new['private_namespace']+v[len(old['private_namespace']):] if v.startswith(old['private_namespace']) else v
 if isinstance(v,list):return [relocate(x) for x in v]
 if isinstance(v,dict):return {k:relocate(x) for k,x in v.items()}
 return v
check('all35_complete_member_objects_only_private_namespace_relocation',relocate(old['members'])==new['members'] and new['member_count']==old['member_count']==35)
check('all35_full_original_design_and_science_argv_unchanged',all(a['original_design']==b['original_design'] and all(relocate(a[k])==b[k] for k in ['literal_original_harmonize_command','harmonize_command','munge_command']) for a,b in zip(old['members'],new['members'])))
check('exact_single_source_no_abbreviation_or_additions_all35',all(m['literal_original_harmonize_command'].count('--infile')==1 and not any(x.startswith('--infil') and x!='--infile' for x in m['literal_original_harmonize_command']) for m in new['members']))
stable=['baseline_support_package','baseline_execution_plan_sha256','baseline_dependency_sha256','baseline_relocated_dependency_sha256','baseline_input_sha256','original_190_jobs','core_precision_adjudication','python','ldsc_dir','reference_prefix','environment','guard','runtime_qualification','qualified_runtime_receipt','qualified_runtime_receipt_sha256','harmonization_python','component_review_seal_sha256','bounded_adapter_review_seal_sha256','scientific_adapter_qualification','global_resource_ledger_sha256','estimator_calls','scientific_membership_or_threshold_changes','automatic_retry']
check('all_baseline190_guard_runtime_and_qualification_keys_unchanged',all(old[k]==new[k] for k in stable))
for p,h in new['dependencies_sha256'].items():
 f=Path(p)
 if f.name in ['hm3_grch37_variant_map.tsv.gz','hg38ToHg19.over.chain.gz','w_hm3.snplist']:
  UNREAD[p]={'expected_sha256':h,'bytes':f.stat().st_size,'body_read':False};continue
 bind(p,h)
check('228_deps225_code_metadata3_reference_stat_only',len(new['dependencies_sha256'])==228 and len(UNREAD)==3)
for name,pin in [('core_large35_bounded_controller_review_seal_v1.json','06ba91909321a1e7d0f481bf565f5f0ad938dabd139995f803d1fe7ebdf098e1'),('core_large35_bounded_controller_review_seal_v2.json','9ef7aa0a983033cf1ae0cac6bc232f5804f6a16e64f3de3ae32cf56d0d15e12e'),('independent_spool_disposal_seal_v1.json','642fa6e1fbd12d8e36bae5d7172db184ec3fed0851cfae9c31a6103a8e8e7217'),('independent_spool_disposal_seal_v2.json','0f892f5b5ed3897c689162b5f3a7f91ace5bf609cdfe578c59ef76c55f346fa7')]:
 p=R/name;bind(p,pin);seal=json.loads(p.read_text());check(name+'_preserved_rejection_seal_bound',new['dependencies_sha256'].get(str(p))==pin)
 check(name+'_all_artifacts_unchanged_and_bound',all(new['dependencies_sha256'].get(f)==(v['sha256'] if isinstance(v,dict) else v)==bind(f,v['sha256'] if isinstance(v,dict) else v) for f,v in seal.get('artifacts',seal.get('file_sha256',{})).items()))
oldtext=(S/'78_run_core_large35_bounded_replay_v2.py').read_text();newtext=(S/'84_run_core_large35_bounded_replay_v3.py').read_text()
for a,b in [('large35_bounded_replay_v2','large35_bounded_replay_v3'),('core_large35_pending_v2','core_large35_pending_v3'),('core_large35_terminal_seal_v2','core_large35_terminal_seal_v3'),('core_large35_execution_receipt_v2','core_large35_execution_receipt_v3'),('79_validate_core_large35_bounded_replay_v2.py','85_validate_core_large35_bounded_replay_v3.py'),('80_cleanup_core_bounded_spool_v2.py','86_cleanup_core_bounded_spool_v3.py')]:oldtext=oldtext.replace(a,b)
check('entire84_controller_exact78_only_prescribed_filename_replacements',oldtext==newtext)
check('175_exact_command_and_three_dynamic_digest_arg_logic_unchanged',"len(master['worker_receipt_sha256'])!=175" in newtext and all(newtext.count(k)==1 for k in ['--expected-candidate-sha256','--expected-columns-sha256','--expected-comparison-worker-sha256']))
check('delayed_metadata_regular_capture_and_immutable_copy_final_gate_inherited',"master['generated_metadata_sha256'][trait]={str(f):sha(f) for f in generated}" in newtext and 'CORE_FROZEN_GENERATED_METADATA_COPY_CHANGED' in newtext)
def funcs(p):return {n.name:ast.dump(n,include_attributes=False) for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
a=funcs(S/'79_validate_core_large35_bounded_replay_v2.py');b=funcs(S/'85_validate_core_large35_bounded_replay_v3.py')
check('QC_and_entire_harmonized_science_functions_AST_identical',a['qc']==b['qc'] and a['compare_harmonized']==b['compare_harmonized'])
v=(S/'85_validate_core_large35_bounded_replay_v3.py').read_text()
check('four_current_identity_capture_before_source_parse_or_QC_or_streams',v.index('consumed_outputs=regular_hashes')<v.index("source=json.loads(Path(member['source_gate_receipt']).read_text())")<v.index("r['harmonized_comparison']=compare_harmonized")<v.index("observed,steps=qc("))
check('two_frozen_sealed_baseline_hashes_before_comparison',v.index('consumed_baselines=')<v.index('if regular_hashes(consumed_baselines)!=consumed_baselines:')<v.index("r['harmonized_comparison']=compare_harmonized"))
check('leaf_and_parent_regular_no_symlink_for_all_captured_files',"path.is_symlink() or not path.is_file() or any(parent.is_symlink() for parent in path.parents)" in v)
check('consumed_current_and_baseline_hashes_after_comparison_and_before_persistence',v.count('regular_hashes(consumed_outputs)!=consumed_outputs or regular_hashes(consumed_baselines)!=consumed_baselines')==2 and v.index('COMPARED_OUTPUT_SOURCE_QC_OR_BASELINE_CHANGED_AFTER_CONSUMPTION')<v.index('write_new(a.out,r)'))
check('result_names_immutable_consumed_map',"r['output_identity_sha256']=consumed_outputs" in v)
edge=R/'core_large35_validator_identity_edge_receipt_v3.json';bind(edge);e=json.loads(edge.read_text());check('16_actual85_cases_no_false_accepts',len(e['cases'])==16 and e['source_sha256']==bind(S/'85_validate_core_large35_bounded_replay_v3.py') and all(c['validator_accepts']==(c['case']=='healthy') for c in e['cases']))
prep=P/'logs/core_large35_preparation_harness_correction_v3.json';bind(prep);check('unfrozen_preparation_failure_evidence_retained',prep.is_file())
check('all_bound_plan_sources_metadata_seals_unchanged',IDS=={p:sha(p) for p in IDS})
receipt={'status':'PASS_NARROW_V3_CORRECTION_PENDING_INDEPENDENT_DISPOSAL_DEPENDENCY_NOT_EXECUTION_ADMISSION','plan_path':str(Q),'plan_sha256':PIN,'check_count':len(CHECKS),'checks':CHECKS,'before':IDS,'after':{p:sha(p) for p in IDS},'unread_reference_bodies':UNREAD,'16_actual_validator_cases':str(edge),'inherited_v2_186_checks_not_broadly_rerun':True,'inherited_v1_22_controller_scenarios_not_rerun':True,'actual_workers_or_mutex_or_fits':0,'production_body_reads':0,'full_runtime_tree_rehashes':0,'execution_admitted':False,'disposal86_independent_qualification':'Must bind separate component PASS before overall eligibility.'}
with (R/'core_large35_bounded_successor_binding_receipt_v3.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({'status':receipt['status'],'checks':len(CHECKS),'bound_code_metadata_identities':len(IDS)}))
