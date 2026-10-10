"""Narrow v6 metadata, exact new terminal callbacks and invented tiny files."""
import ast
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
from types import SimpleNamespace
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN=SSD/'extension_pipeline_replay_v6/extension_pipeline_replay_plan_v6.json'
OLD=SSD/'extension_pipeline_replay_v5/extension_pipeline_replay_plan_v5.json'
RAW=P/'manifests/extension_raw_acquisition_plan_v4_8.json'
ROOT=R/'genomicsem_extension_pipeline_v6_fixture_controls';ROOT.mkdir(exist_ok=False)
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
checks=[]
def check(value,label,detail=None):
    checks.append({'control':label,'pass':bool(value),'detail':detail});assert value,label
def module(path,name):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
PINS={str(PLAN):'d0fd223d468f26583bfd2bd477a51548646b57e91275d7733eb1ae814b94b798',str(OLD):'766245a2fc3bef811fdd1217cb1fbff47175c7d518b4337900ca9e5fb17c7866',str(RAW):'77bf2da98b851ac1aa9d414b2b695fb04a65b0c4ca987f7d85db2097e9e7ef8b',str(S/'sensitivity_executor_v4_4.py'):'7ac802596f29adfacb017f406676e6a30ee00bf890904976afda7f88a3ffa9c5',str(S/'terminal_commit_common_v2.py'):'9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd'}
for path,h in PINS.items():check(sha(path)==h,'pin_'+Path(path).name)
old,plan,raw=[json.loads(x.read_text()) for x in [OLD,PLAN,RAW]]
for sealpath in [R/'genomicsem_extension_pipeline_v5_preflight_seal.json',R/'genomicsem_acquisition_v8_preflight_seal.json']:
    for path,meta in json.loads(sealpath.read_text())['artifacts'].items():check(sha(path)==meta['sha256'],'prior_sealed_'+Path(path).name)
compat_artifacts={}
for line in (R/'independent_historical_munge_text_v1.sha256').read_text().splitlines():
    h,path=line.split('  ',1);compat_artifacts[path]=h;check(sha(path)==h,'compatibility_sealed_'+Path(path).name)
check(plan['dependencies_sha256'][str(R/'independent_historical_munge_text_v1.sha256')]=='8a801dd0160bbd0a86c4a19e373d0da39224a32d1550c816cc8d92a9702d3171','later43_compatibility_seal_now_bound_by_plan')
check(len(old['members'])==len(plan['members'])==len(raw['members'])==100 and plan['total_owned_commands']==400 and plan['estimator_calls']==0,'same100_members400_commands_no_fit')
for a,b,m in zip(old['members'],plan['members'],raw['members']):
    expected=copy.deepcopy(a)
    for k,v in list(expected.items()):
        if isinstance(v,str) and v.startswith(str(OLD.parent)+'/'):expected[k]=str(PLAN.parent)+v[len(str(OLD.parent)):]
        elif k in ['harmonize_command','munge_command']:
            expected[k]=[str(PLAN.parent)+x[len(str(OLD.parent)):] if x.startswith(str(OLD.parent)+'/') else m['body_path'] if x==a['raw'] else x for x in v]
    expected['raw']=m['body_path'];expected['acquisition_member']=m
    fixed=raw['reused_checkpoints'].get(str(m['index']))
    expected['acquisition_receipt']=fixed['receipt_path'] if fixed else str(Path(raw['pending_path']).parent/'receipts'/(m['extension_trait_id']+'.json'))
    expected['acquisition_origin_plan_sha256']=fixed['original_acquisition_plan_sha256'] if fixed else PINS[str(RAW)]
    check(expected==b,'exact_raw_origin_and_namespace_only_member_'+str(b['index']))
    check({k:v for k,v in a['acquisition_member'].items() if k!='body_path'}=={k:v for k,v in m.items() if k!='body_path'},'original_source_tuple_'+str(b['index']))
changed={'schema','prepared_utc','members','dependencies_sha256','acquisition_plan','acquisition_plan_sha256','reservation_arithmetic','execution_preconditions','resource_preflight','acquisition_execution_binding','acquisition_family_receipt','acquisition_operational_identity_by_origin','acquisition_family_terminal_contract','global_resource_ledger_path','global_resource_ledger_sha256'}
for k,v in old.items():
    if k not in changed:check(plan[k]==v,'unchanged_plan_field_'+k)
check(all(plan['dependencies_sha256'].get(k)==h for k,h in old['dependencies_sha256'].items()),'all468_old_declarations_retained')
check(all(plan['dependencies_sha256'].get(k)==h for k,h in raw['bound_sources'].items()),'all169_v8_static_bindings_retained')
contract=plan['acquisition_family_terminal_contract']
check((contract['reused_exact_source_count'],contract['new_exact_source_count'])==(4,96) and contract['pending_path']==raw['pending_path'] and contract['terminal_seal_path']==raw['terminal_seal_path'],'exact_v8_4plus96_contract')
check(plan['acquisition_operational_identity_by_origin'][PINS[str(RAW)]]=={'source_folder':str(Path(raw['pending_path']).parent),'curl':'/usr/bin/curl','per_body_seconds_limit':7200,'resume_offset_by_index':{'5':559470675}},'exact_v8_source5_origin_context')
for origin,v in old['acquisition_operational_identity_by_origin'].items():check(plan['acquisition_operational_identity_by_origin'][origin]==v,'preserved_individual_origin_'+origin)
check(plan['global_resource_ledger_sha256']==raw['global_resource_ledger_sha256'] and plan['global_resource_ledger_path']==raw['global_resource_ledger_path'],'exact_current_ledger_binding')
metadata=dict(PINS);metadata.update(plan['acquisition_execution_binding']['sha256']);metadata.update({k:h for k,h in plan['dependencies_sha256'].items() if k not in old['dependencies_sha256']});metadata.update(compat_artifacts)
for n in ['46_prepare_extension_pipeline_replay_v6.py','47_run_extension_pipeline_replay_v6.py','extension_replay_common_v5.py','extension_replay_validate_v6.py']:metadata[str(S/n)]=plan['dependencies_sha256'][str(S/n)]
for path,h in metadata.items():check(Path(path).suffix not in ['.gz','.bgz','.bed','.bim','.fam'] and path not in plan['archived_input_sha256'] and sha(path)==h,'fresh_static_nonbody_'+Path(path).name)

funcs=lambda p:{n.name:n for n in ast.parse(p.read_text()).body if isinstance(n,ast.FunctionDef)}
common_old,common_new=funcs(S/'extension_replay_common_v4.py'),funcs(S/'extension_replay_common_v5.py')
for name in common_old:
    if name!='acquisition_family_gate':check(ast.dump(common_old[name],include_attributes=False)==ast.dump(common_new[name],include_attributes=False),'exact_common_science_or_origin_'+name)
class CommonCounts(ast.NodeTransformer):
    def visit_Constant(self,n):
        if isinstance(n.value,str):n.value=n.value.replace('v8 four-reused/96-new','v7 three-reused/97-new').replace('FROZEN4_PLUS96','FROZEN3_PLUS97')
        return n
    def visit_Tuple(self,n):
        if all(isinstance(x,ast.Constant) for x in n.elts) and [x.value for x in n.elts]==[4,96]:n.elts=[ast.Constant(3),ast.Constant(97)]
        return self.generic_visit(n)
check(ast.dump(common_old['acquisition_family_gate'],include_attributes=False)==ast.dump(CommonCounts().visit(copy.deepcopy(common_new['acquisition_family_gate'])),include_attributes=False),'family_gate_only_strict4_96_and_doc')
oldrun,newrun=funcs(S/'47_run_extension_pipeline_replay_v5.py'),funcs(S/'47_run_extension_pipeline_replay_v6.py')
for name in oldrun:
    if name not in ['run','main']:check(ast.dump(oldrun[name],include_attributes=False)==ast.dump(newrun[name],include_attributes=False),'unchanged_runner_helper_'+name)
class VersionOnly(ast.NodeTransformer):
    def visit_Constant(self,n):
        if isinstance(n.value,str):n.value=n.value.replace('extension_replay_validate_v6.py','extension_replay_validate_v5.py').replace('extension_pipeline_replay_v6','extension_pipeline_replay_v5').replace('pipeline_attempt_v6','pipeline_attempt_v5').replace('pipeline_pending_v6','pipeline_pending_v5').replace('pipeline_terminal_seal_v6','pipeline_terminal_seal_v5').replace('pipeline_execution_receipt_v6','pipeline_execution_receipt_v5')
        return n
    def visit_ImportFrom(self,n):
        if n.module=='extension_replay_common_v5':n.module='extension_replay_common_v4'
        return n
check(ast.dump(ast.parse((S/'extension_replay_validate_v5.py').read_text()),include_attributes=False)==ast.dump(VersionOnly().visit(ast.parse((S/'extension_replay_validate_v6.py').read_text())),include_attributes=False),'validator_entire_AST_only_common_import')
check(ast.dump(oldrun['main'],include_attributes=False)==ast.dump(VersionOnly().visit(copy.deepcopy(newrun['main'])),include_attributes=False),'executor_default_plan_matches_v6')
oldjobs=next(n for n in ast.walk(oldrun['run']) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='jobs' for t in n.targets))
newjobs=next(n for n in ast.walk(newrun['run']) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='jobs' for t in n.targets))
check(ast.dump(oldjobs,include_attributes=False)==ast.dump(newjobs,include_attributes=False),'four_original_job_roles_and_prefixes_unchanged')
for name in ['limits','validation']:
    a=next(n for n in ast.walk(oldrun['run']) if isinstance(n,ast.FunctionDef) and n.name==name);b=next(n for n in ast.walk(newrun['run']) if isinstance(n,ast.FunctionDef) and n.name==name)
    check(ast.dump(a,include_attributes=False)==ast.dump(VersionOnly().visit(copy.deepcopy(b)),include_attributes=False),'unchanged_inner_science_guard_'+name)

sys.path.insert(0,str(S));runner=module(S/'47_run_extension_pipeline_replay_v6.py','_pipeline6_fixture');common=module(S/'extension_replay_common_v5.py','_common5_fixture')
# Compile exact new nested callback functions without invoking run/worker/locks.
inner=[n for n in newrun['run'].body if isinstance(n,ast.FunctionDef) and n.name in ['regular_output_hashes','final_identity_gate']]
home=ROOT/'consumed';home.mkdir();(home/'receipts').mkdir();(home/'workers').mkdir();(home/'outputs').mkdir();(home/'worker_artifacts').mkdir()
def save(path,value):path.write_text(json.dumps(value,sort_keys=True)+'\n')
binding_plan=home/'plan.json';admit=home/'admission.json';save(binding_plan,{'synthetic':True});save(admit,{'synthetic':True})
a=SimpleNamespace(plan=binding_plan,plan_sha=sha(binding_plan),admission=admit,admission_sha=sha(admit))
tiny={'members':[]};master={'worker_receipt_sha256':{},'completed_members':[],'discovered_acquisition_receipt_sha256':{},'checkpoint_binding_sha256':{},'acquisition_family_terminal_evidence':{'fixture_only':True},'historical190_gate_receipts':{'fixture_only':True}}
worker_files=[];output_files=[]
for j in range(100):
    trait='fixture_%03d'%j;m={'extension_trait_id':trait};output_map={}
    for label in ['harmonized','harmonization_qc','harmonization_receipt','munged','source_gate_receipt','comparison_receipt']:
        file=home/'outputs'/(trait+'__'+label+'.txt');file.write_text('only invented '+label+'\n');m[label]=str(file);output_map[str(file)]=sha(file);output_files.append(file)
    tiny['members'].append(m)
    bp=home/'receipts'/(trait+'.checkpoint_binding.json');save(bp,{'fixture_only':True});master['checkpoint_binding_sha256'][str(bp)]=sha(bp)
    master['completed_members'].append({'trait':trait,'output_sha256':output_map,'comparison_sha256':sha(m['comparison_receipt']),'new_munged_sha256':sha(m['munged'])})
    for label in ['source','harmonize','munge','compare']:
        artifact=home/'worker_artifacts'/(trait+'__'+label+'.stdout');artifact.write_text('synthetic consumed stdout\n')
        file=home/'workers'/(trait+'__'+label+'.worker.json');save(file,{'status':'WORKER_COMPLETE_VERIFIED','plan_sha256':a.plan_sha,'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[]},'output_sha256':{str(artifact):sha(artifact)}})
        master['worker_receipt_sha256'][str(file)]=sha(file);worker_files.append(file)
env={'Path':Path,'sha':sha,'json':json,'ownership':[True,[]],'physical_mount':lambda:None,'check_bindings':lambda *x,**k:None,'plan':tiny,'a':a,'admission':{},'admission_gate':lambda *x:None,'acquisition_family_gate':lambda *x,**k:{'fixture_only':True},'master':master,'protected':SimpleNamespace(baseline_gate=lambda p:{'fixture_only':True}),'OUT':home,'checkpoint_binding_gate':lambda *x:None,'__file__':str(S/'47_run_extension_pipeline_replay_v6.py')}
exec(compile(ast.fix_missing_locations(ast.Module(body=copy.deepcopy(inner),type_ignores=[])),str(S/'47_run_extension_pipeline_replay_v6.py'),'exec'),env)
final=env['final_identity_gate'];final();check(True,'healthy_full400_worker600_consumed_output_callback')
def reject(label,mutation,restore,error=None):
    mutation()
    try:
        try:final()
        except (RuntimeError,FileNotFoundError,KeyError) as e:check(error is None or error in str(e),label,type(e).__name__+': '+str(e))
        else:raise AssertionError('Consumed callback accepted '+label)
    finally:restore()
file=worker_files[0];original=file.read_bytes()
reject('missing_old_worker_receipt_rejects',lambda:file.unlink(),lambda:file.write_bytes(original))
reject('changed_old_worker_receipt_rejects',lambda:file.write_bytes(original+b' '),lambda:file.write_bytes(original),'REGULAR_OUTPUT_CHANGED')
file=output_files[1];original=file.read_bytes()
reject('missing_old_QC_consumed_output_rejects',lambda:file.unlink(),lambda:file.write_bytes(original))
reject('changed_old_QC_consumed_output_rejects',lambda:file.write_bytes(original+b'changed'),lambda:file.write_bytes(original),'REGULAR_OUTPUT_CHANGED')
backup=home/'same_bytes_backup.txt';backup.write_bytes(original)
def make_symlink():file.unlink();file.symlink_to(backup)
def restore_file():file.unlink();file.write_bytes(original)
reject('same_bytes_output_symlink_rejects',make_symlink,restore_file,'REGULAR_OUTPUT_CHANGED')
first=worker_files[0];bytes0=first.read_bytes();backup2=home/'same_worker_backup.json';backup2.write_bytes(bytes0)
reject('same_bytes_worker_receipt_symlink_rejects',lambda:(first.unlink(),first.symlink_to(backup2)),lambda:(first.unlink(),first.write_bytes(bytes0)),'REGULAR_OUTPUT_CHANGED')
artifact=home/'worker_artifacts/fixture_000__source.stdout';art=artifact.read_bytes()
reject('missing_old_worker_consumed_stdout_rejects',lambda:artifact.unlink(),lambda:artifact.write_bytes(art))
key=str(worker_files[0]);digest=master['worker_receipt_sha256'][key]
reject('399_worker_map_rejects',lambda:master['worker_receipt_sha256'].pop(key),lambda:master['worker_receipt_sha256'].__setitem__(key,digest),'EXACT400')
extra=home/'workers/extra.worker.json';extra.write_bytes(bytes0)
reject('401_worker_map_rejects',lambda:master['worker_receipt_sha256'].__setitem__(str(extra),sha(extra)),lambda:master['worker_receipt_sha256'].pop(str(extra)),'EXACT400')
parent=home/'worker_artifacts';parent_real=home/'worker_artifacts_real'
def symlink_parent():parent.rename(parent_real);parent.symlink_to(parent_real,target_is_directory=True)
def restore_parent():parent.unlink();parent_real.rename(parent)
reject('same_bytes_parent_directory_symlink_rejects',symlink_parent,restore_parent,'REGULAR_OUTPUT_CHANGED')
for label,field,value in [('failed_worker_status','status','WORKER_FAILED_PRESERVED'),('wrong_worker_plan','plan_sha256','wrong'),('worker_cleanup_unverified','owned_cleanup_verified',False),('worker_group_member_survives','process_group_teardown',{'remaining_group_members':[123]})]:
    file=worker_files[0];oldbytes=file.read_bytes();oldhash=master['worker_receipt_sha256'][str(file)];data=json.loads(oldbytes);data[field]=value
    def change(file=file,data=data):save(file,data);master['worker_receipt_sha256'][str(file)]=sha(file)
    def restore(file=file,oldbytes=oldbytes,oldhash=oldhash):file.write_bytes(oldbytes);master['worker_receipt_sha256'][str(file)]=oldhash
    reject(label+'_rejects',change,restore,'NOT_COMPLETE_OR_REAPED')
env['ownership'][0]=False
try:final()
except RuntimeError as e:check('OWNED_WORKERS_NOT_CLEANED' in str(e),'uncleared_outer_ownership_rejects')
finally:env['ownership'][0]=True
final();check(True,'restored_full400_callback_passes')

fixture=module(P/'statistical_validation/genomicsem_extension_checkpoint_identity_controls_v3.py','_family6_fixture')
for good in [True,False]:
    f=fixture.family(ROOT/('family4_96' if good else 'family_wrong3_97'))
    f['plan']['acquisition_family_terminal_contract'].update(reused_exact_source_count=4,new_exact_source_count=96)
    f['primary'].update(reused_exact_source_count=4 if good else 3,new_exact_source_count=96 if good else 97)
    f['s']['primary_receipt_sha256']={str(x):fixture.save(x,f['primary']) for x in f['copies']};fixture.save(f['seal'],f['s'])
    if good:
        with patch.object(common,'acquisition_operational_gate',lambda p:None):evidence=common.acquisition_family_gate(f['plan'],f['map'],require_complete=True)
        check(len(evidence['source_receipt_sha256'])==100,'correct4_96_full100_family_accepts')
        f['plan']['acquisition_family_terminal_contract'].update(reused_exact_source_count=3,new_exact_source_count=97)
        try:common.acquisition_family_gate(f['plan'],f['map'],require_complete=True)
        except RuntimeError as e:check('FROZEN4_PLUS96' in str(e),'old3_97_contract_rejects')
        else:raise AssertionError('old contract accepted')
    else:
        try:common.acquisition_family_gate(f['plan'],f['map'],require_complete=True)
        except RuntimeError as e:check('PRIMARY_FULL100' in str(e),'old3_97_primary_counts_reject')
        else:raise AssertionError('old primary counts accepted')
for path,h in metadata.items():check(sha(path)==h,'metadata_unchanged_after_'+Path(path).name)
result={'status':'PASS_NARROW_V6_NEW_CONSUMED_OUTPUT_AND_WORKER_TERMINAL_GATE','checks':checks,'check_count':len(checks),'fresh_metadata_before_after':metadata,'fresh_metadata_count':len(metadata),'source_plan_sha256':PINS[str(PLAN)],'all_dependency_declarations':plan['dependencies_sha256'],'new_dependency_declarations':{k:h for k,h in plan['dependencies_sha256'].items() if k not in old['dependencies_sha256']},'inherited_v5_controls':541,'inherited_compatibility_controls':43,'real_producer_v8_terminal_not_read_or_credited':True,'terminal_callback_fixture_only':'Exact extracted regular_output_hashes and final_identity_gate; acquisition/baseline/identity/mount callbacks labelled synthetic mocks.400worker/600consumed outputs and100checkpoint bindings were own invented tiny files.','execution_admitted':False,'actual_pipeline_or_full_raw_family_certified':False,'production_raw_body_reads':0,'reference_body_reads':0,'production_mutex_operations':0,'real_worker_or_fit_operations':0,'peak_reviewer_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
with (R/'genomicsem_extension_pipeline_v6_controls_receipt.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':result['status'],'checks':len(checks),'fresh_metadata':len(metadata),'peak_RSS_bytes':result['peak_reviewer_RSS_bytes']}))
