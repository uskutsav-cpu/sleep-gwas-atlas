"""Narrow frozen v5 command/runtime bindings, metadata and private tiny fixtures."""
import ast
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import tempfile
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
OLD=SSD/'extension_pipeline_replay_v4/extension_pipeline_replay_plan_v3.json'
PLAN=SSD/'extension_pipeline_replay_v5/extension_pipeline_replay_plan_v5.json'
PROVENANCE=P/'source_provenance/historical_munge_header_text_compatibility_v1.json'
PINNED={str(PLAN):'766245a2fc3bef811fdd1217cb1fbff47175c7d518b4337900ca9e5fb17c7866',str(OLD):'7c941a230287e04b42526b215b1f13373b850885fd9c8abacd2f0ab83731c690',str(PROVENANCE):'8e9388c92d6801c299a1d61b7d59fa91da789fd2c2422f4f4983e721dbcdb41a',str(R/'genomicsem_extension_pipeline_v4_preflight_seal.json'):'fe2d11116b038511f997a90f2c1dfa509df49651b4eba3e3682ddba95df394fc',str(S/'sensitivity_executor_v4_4.py'):'7ac802596f29adfacb017f406676e6a30ee00bf890904976afda7f88a3ffa9c5',str(S/'terminal_commit_common_v2.py'):'9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd',str(S/'extension_replay_common_v4.py'):'319fd5da9dbac1ef197aeb107c566394d8a43bff25241af0c63e6053ef19c2e7'}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
checks=[]
def check(condition,label,detail=None):
    checks.append({'control':label,'pass':bool(condition),'detail':detail})
    assert condition,label
for path,h in PINNED.items():check(sha(path)==h,'pin_'+Path(path).name)
old,plan,compat=[json.loads(x.read_text()) for x in [OLD,PLAN,PROVENANCE]]
oldseal=json.loads((R/'genomicsem_extension_pipeline_v4_preflight_seal.json').read_text())
for path,meta in oldseal['artifacts'].items():check(sha(path)==meta['sha256'],'preserved_v4_artifact_'+Path(path).name)
compat_seal=R/'independent_historical_munge_text_v1.sha256'
compat_artifacts={}
for line in compat_seal.read_text().splitlines():
    h,path=line.split('  ',1);compat_artifacts[path]=h
    check(sha(path)==h,'independent_compatibility_artifact_'+Path(path).name)
check(len(compat_artifacts)==5,'exact_five_independent_compatibility_artifacts')
compat_artifacts[str(compat_seal)]=sha(compat_seal)
review=json.loads((R/'independent_historical_munge_text_receipt_v1.json').read_text())
check(review['passed_control_count']==43 and review['production_GWAS_body_reads']==0 and review['estimator_calls']==0,'sealed43_qualification_not_actual_replay')
check(review['scientific_AST_identical_except_read_header'] is True and review['historical_per_trait_binary_attestation'] is False and review['scientific_full_replay_certified'] is False,'compatibility_scope_not_overclaimed')
check(plan['python']==old['python'] and review['consumed_metadata_and_code_before_after_sha256'][plan['python']]['sha256']==plan['dependencies_sha256'][plan['python']],'same_declared_munge_runtime')
check(compat['candidate_is_newer_LDSC_implementation'] is False and compat['candidate_is_historical_per_trait_binary_attestation'] is False and compat['scientific_replay_certified'] is False,'preexisting_compatibility_not_newer_LDSC')
check(plan['historical_munge_compatibility_provenance']==str(PROVENANCE) and plan['historical_munge_compatibility_provenance_sha256']==PINNED[str(PROVENANCE)] and plan['historical_munge_qualification']==compat['qualification'],'exact_compatibility_provenance_and_qualification')

check(len(plan['members'])==len(old['members'])==plan['member_count']==100,'exact_ordered100')
check(plan['total_owned_commands']==400 and plan['heavy_preprocessing_commands']==200 and plan['owned_stdlib_validation_commands']==200 and plan['estimator_calls']==0,'exact_original400_commands_no_estimator')
for a,b in zip(old['members'],plan['members']):
    expected=copy.deepcopy(a)
    for k,v in list(expected.items()):
        if isinstance(v,str) and v.startswith(str(OLD.parent)+'/'):expected[k]=str(PLAN.parent)+v[len(str(OLD.parent)):]
        elif k in ['harmonize_command','munge_command']:
            expected[k]=[str(PLAN.parent)+x[len(str(OLD.parent)):] if x.startswith(str(OLD.parent)+'/') else x for x in v]
    check(expected['munge_command'].count(compat['pinned_unmodified_munger'])==1,'one_original_munger_path_'+str(b['index']))
    expected['munge_command']=[compat['candidate_executable'] if x==compat['pinned_unmodified_munger'] else x for x in expected['munge_command']]
    check(expected==b,'exact_namespace_and_compat_munger_only_member_'+str(b['index']))
    check(b['munge_command'][:4]==[plan['python'],'-u','-B',compat['candidate_executable']],'exact_declared_munge_command_prefix_'+str(b['index']))
changed={'schema','prepared_utc','members','dependencies_sha256','execution_preconditions','resource_preflight'}
for k,v in old.items():
    if k not in changed:check(plan[k]==v,'unchanged_plan_field_'+k)
check(set(plan)-set(old)=={'preserved_v4_plan','preserved_v4_plan_sha256','historical_munge_compatibility_provenance','historical_munge_compatibility_provenance_sha256','historical_munge_qualification'},'only_prespecified_new_plan_fields')
check(all(plan['dependencies_sha256'].get(k)==v for k,v in old['dependencies_sha256'].items()),'all437_preserved_v4_declarations')
check(plan['preserved_v4_plan']==str(OLD) and plan['preserved_v4_plan_sha256']==PINNED[str(OLD)],'preserved_blocked_v4_plan_exact')

metadata=dict(PINNED);metadata.update(compat_artifacts)
newdeps={k:v for k,v in plan['dependencies_sha256'].items() if k not in old['dependencies_sha256']}
metadata.update(newdeps)
for path,identity in compat['source_and_copy_identity'].items():
    check(plan['dependencies_sha256'].get(path)==identity['sha256'],'declared_exact_historical_copy_'+Path(path).name)
    metadata[path]=identity['sha256']
check(len(compat['source_and_copy_identity'])==22,'exact11_source_copy_pairs')
for path,h in compat['unchanged_pinned_imported_module_sha256'].items():check(plan['baseline_dependency_sha256'].get(path)==h,'unchanged_estimator_module_binding_'+Path(path).name)
for path,h in metadata.items():
    check(Path(path).suffix not in ['.gz','.bgz','.bed','.bim','.fam'] and path not in plan['archived_input_sha256'] and path not in {m['raw'] for m in plan['members']},'non_production_body_metadata_'+Path(path).name)
    check(sha(path)==h,'fresh_code_metadata_'+Path(path).name)

class VersionPaths(ast.NodeTransformer):
    def visit_Constant(self,n):
        if isinstance(n.value,str):
            for a,b in [('extension_pipeline_replay_v5','extension_pipeline_replay_v4'),('pipeline_attempt_v5.json','pipeline_attempt_v4.json'),('pipeline_pending_v5.json','pipeline_pending_v4.json'),('pipeline_terminal_seal_v5.json','pipeline_terminal_seal_v4.json'),('pipeline_execution_receipt_v5.json','pipeline_execution_receipt_v4.json'),('extension_replay_validate_v5.py','extension_replay_validate_v4.py'),('extension_pipeline_replay_plan_v5.json','extension_pipeline_replay_plan_v4.json')]:n.value=n.value.replace(a,b)
        return n
oldast=ast.parse((S/'47_run_extension_pipeline_replay_v4.py').read_text())
newast=VersionPaths().visit(ast.parse((S/'47_run_extension_pipeline_replay_v5.py').read_text()))
check(ast.dump(oldast,include_attributes=False)==ast.dump(newast,include_attributes=False),'whole_runner_AST_only_prespecified_v5_paths_default')
check((S/'extension_replay_validate_v5.py').read_bytes()==(S/'extension_replay_validate_v4.py').read_bytes(),'validator_byte_identical_common4_import_retained')
tree=ast.parse((S/'47_run_extension_pipeline_replay_v5.py').read_text())
default=next(x for x in ast.walk(tree) if isinstance(x,ast.Call) and isinstance(x.func,ast.Attribute) and x.func.attr=='add_argument' and x.args and isinstance(x.args[0],ast.Constant) and x.args[0].value=='--plan')
value=next(x.value for x in default.keywords if x.arg=='default')
check(isinstance(value,ast.BinOp) and value.right.value==PLAN.name and 'extension_pipeline_replay_v5' in (S/'47_run_extension_pipeline_replay_v5.py').read_text(),'prepared_default_path_now_matches')

sys.path.insert(0,str(S))
sp=importlib.util.spec_from_file_location('_v5_pipeline_review',S/'47_run_extension_pipeline_replay_v5.py')
runner=importlib.util.module_from_spec(sp);sp.loader.exec_module(runner)
with tempfile.TemporaryDirectory(prefix='extension-v5-binding-fixtures-') as tmp:
    out=Path(tmp)
    tiny={'members':[{'harmonized':str(out/'h.gz'),'harmonization_qc':str(out/'qc.tsv'),'harmonization_receipt':str(out/'h.json'),'munged':str(out/'m.gz'),'munged_prefix':str(out/'m'),'source_gate_receipt':str(out/'s.json'),'comparison_receipt':str(out/'c.json')}]}
    with patch.object(runner,'OUT',out):
        runner.no_preserved_outputs_gate(tiny);check(True,'tiny_v5_private_namespace_accepts')
        wrong=copy.deepcopy(tiny);wrong['members'][0]['munged']=str(out.parent/'escape.gz')
        try:runner.no_preserved_outputs_gate(wrong)
        except RuntimeError as e:check('OUTPUT_ESCAPES' in str(e),'tiny_v5_namespace_escape_rejects')
        else:raise AssertionError('namespace escape accepted')
        (out/'m.partial').write_text('own fixture')
        try:runner.no_preserved_outputs_gate(tiny)
        except RuntimeError as e:check('PREFIX_PRESERVED' in str(e),'tiny_v5_existing_prefix_rejects')
        else:raise AssertionError('existing prefix accepted')
    admission={'execution_admitted':True,'plan_sha256':PINNED[str(PLAN)],'executor_sha256':sha(S/'47_run_extension_pipeline_replay_v5.py'),'checkpoint_receipt_binding_policy':plan['checkpoint_receipt_binding_policy'],'acquisition_execution_plan_sha256':plan['acquisition_plan_sha256'],'acquisition_executor_sha256':plan['acquisition_execution_binding']['executor_sha256'],'acquisition_root_admission_sha256':plan['acquisition_execution_binding']['sha256'][plan['acquisition_execution_binding']['root_admission']],'independent_review_sha256':compat_artifacts}
    with patch.object(runner,'acquisition_operational_gate',lambda p:None):
        runner.admission_gate(plan,admission,PINNED[str(PLAN)],sha(S/'47_run_extension_pipeline_replay_v5.py'));check(True,'tiny_admission_accepts_exact_five_compatibility_artifact_and_seal_binding')
        bad=copy.deepcopy(admission);bad['independent_review_sha256'][str(compat_seal)]='0'*64
        try:runner.admission_gate(plan,bad,PINNED[str(PLAN)],sha(S/'47_run_extension_pipeline_replay_v5.py'))
        except RuntimeError as e:check('REVIEW_CHANGED' in str(e),'tiny_admission_changed_compatibility_seal_rejects')
        else:raise AssertionError('Changed seal accepted')

for path,h in metadata.items():check(sha(path)==h,'code_metadata_unchanged_after_'+Path(path).name)
receipt={'status':'PASS_NARROW_EXTENSION_V5_BINDING_PREFLIGHT_NOT_ROOT_ADMISSION','check_count':len(checks),'checks':checks,'fresh_metadata_before_after':metadata,'fresh_metadata_count':len(metadata),'new_dependency_declarations':newdeps,'all_dependency_declarations':plan['dependencies_sha256'],'inherited_v4_controls':333,'inherited_v3_controls':304,'compatibility43_control_seal_path':str(compat_seal),'compatibility43_control_seal_sha256':sha(compat_seal),'required_root_compatibility_review_bindings':compat_artifacts,'plan_does_not_itself_bind_later43_review':not any(x in plan['dependencies_sha256'] for x in compat_artifacts),'root_must_bind_review_artifacts_in_independent_review_sha256':True,'fixture_only':True,'production_raw_body_reads':0,'reference_body_reads':0,'production_mutex_operations':0,'real_worker_or_fit_operations':0,'execution_admitted':False,'full_actual_replay_certified':False,'peak_reviewer_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(R/'genomicsem_extension_pipeline_v5_controls_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'fresh_code_metadata':len(metadata),'dependency_declarations':len(plan['dependencies_sha256']),'compatibility_seal':sha(compat_seal)}))
