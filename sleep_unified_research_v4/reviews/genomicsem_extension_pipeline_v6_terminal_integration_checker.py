"""Exact seal_master/TerminalCommit with new callback, only invented own files."""
import ast
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch

P=Path(__file__).resolve().parents[1];S=P/'scripts';R=P/'reviews'
ROOT=R/'genomicsem_extension_pipeline_v6_terminal_integration_fixtures';ROOT.mkdir(exist_ok=False)
HOME=R/'genomicsem_extension_pipeline_v6_fixture_controls_v1_1/consumed'
sha=lambda x:hashlib.sha256(Path(x).read_bytes()).hexdigest()
source=S/'47_run_extension_pipeline_replay_v6.py';before=sha(source)
sys.path.insert(0,str(S));sp=importlib.util.spec_from_file_location('_v6_integration',source);runner=importlib.util.module_from_spec(sp);sp.loader.exec_module(runner)
tree=ast.parse(source.read_text());run=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='run');inner=[n for n in run.body if isinstance(n,ast.FunctionDef) and n.name in ['regular_output_hashes','final_identity_gate']]
a=SimpleNamespace(plan=HOME/'plan.json',plan_sha=sha(HOME/'plan.json'),admission=HOME/'admission.json',admission_sha=sha(HOME/'admission.json'))
plan={'members':[]};master={'worker_receipt_sha256':{},'completed_members':[],'checkpoint_binding_sha256':{},'discovered_acquisition_receipt_sha256':{},'acquisition_family_terminal_evidence':{'fixture_only':True},'historical190_gate_receipts':{'fixture_only':True}}
for j in range(100):
    t='fixture_%03d'%j;m={'extension_trait_id':t}
    for key in ['harmonized','harmonization_qc','harmonization_receipt','munged','source_gate_receipt','comparison_receipt']:m[key]=str(HOME/'outputs'/(t+'__'+key+'.txt'))
    plan['members'].append(m)
    bp=HOME/'receipts'/(t+'.checkpoint_binding.json');master['checkpoint_binding_sha256'][str(bp)]=sha(bp)
    master['completed_members'].append({'trait':t,'output_sha256':{m[k]:sha(m[k]) for k in m if k!='extension_trait_id'},'comparison_sha256':sha(m['comparison_receipt']),'new_munged_sha256':sha(m['munged'])})
    for label in ['source','harmonize','munge','compare']:
        file=HOME/'workers'/(t+'__'+label+'.worker.json');master['worker_receipt_sha256'][str(file)]=sha(file)
ownership=[True,[]]
env={'Path':Path,'sha':sha,'json':json,'ownership':ownership,'physical_mount':lambda:None,'check_bindings':lambda *x,**k:None,'plan':plan,'a':a,'admission':{},'admission_gate':lambda *x:None,'acquisition_family_gate':lambda *x,**k:{'fixture_only':True},'master':master,'protected':SimpleNamespace(baseline_gate=lambda p:{'fixture_only':True}),'OUT':HOME,'checkpoint_binding_gate':lambda *x:None,'__file__':str(source)}
exec(compile(ast.fix_missing_locations(ast.Module(body=inner,type_ignores=[])),str(source),'exec'),env)
checks=[];evidence=[]
def check(value,label):checks.append({'control':label,'pass':bool(value)});assert value,label
for label in ['healthy','missing_old_worker_after_provisional','missing_old_QC_after_provisional','samebytes_QC_symlink_after_provisional']:
    folder=ROOT/label;folder.mkdir();pending=folder/'pending.json';seal=folder/'seal.json';receipt=folder/'primary.json'
    terminal=runner.TerminalCommit(pending,seal,{'fixture_only':True})
    target=copy.deepcopy(master);target['status']='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS';target['termination_requests']=[]
    protected=SimpleNamespace(stage_resource_gate=lambda *args,**kw:{'fixture_only_resource_PASS':True},check_termination=lambda *args:None,TERMINATION_REQUEST=[],safe_print=lambda *args,**kw:None)
    file=HOME/'workers/fixture_000__source.worker.json' if label=='missing_old_worker_after_provisional' else HOME/'outputs/fixture_000__harmonization_qc.txt';content=file.read_bytes();backup=folder/'same_bytes_backup.txt';backup.write_bytes(content)
    original_write=runner.write_new;mutated=[]
    def save_then_fault(path,value):
        original_write(path,value)
        if Path(path)==receipt and label!='healthy':
            file.unlink();mutated.append(True)
            if label=='samebytes_QC_symlink_after_provisional':file.symlink_to(backup)
    try:
        with patch.object(runner,'write_new',save_then_fault),patch('subprocess.Popen',side_effect=AssertionError('NO_WORKER_IN_REVIEW')):
            runner.seal_master(protected,target,receipt,plan,0,ownership,terminal,env['final_identity_gate'])
        if label=='healthy':
            check(target['status']=='ALL100_RAW_TO_MUNGED_TEMPLATE_REPLAY_PASS' and not pending.exists() and seal.is_file(),label+'_commits')
        else:
            check(mutated==[True],label+'_fault_after_primary_persistence')
            check(target['status']=='FAILED_PRESERVED_NO_AUTOMATIC_RETRY' and pending.is_file(),label+'_durable_PENDING_veto')
            check(Path(str(receipt)+'.failure.json').is_file(),label+'_failure_addendum')
        evidence.append({'case':label,'result_status':target['status'],'pending_retained':pending.exists(),'failure_addendum_exists':Path(str(receipt)+'.failure.json').is_file(),'primary_sha256':sha(receipt),'no_real_workers':True})
    finally:
        if file.is_symlink():file.unlink()
        if not file.exists():file.write_bytes(content)
env['final_identity_gate']();check(True,'all_invented_consumed_files_restored_exact_callback_pass')
check(sha(source)==before,'production_callback_source_unchanged')
result={'status':'PASS_V6_POST_PROVISIONAL_CONSUMED_FILE_FAULT_INTEGRATION','check_count':len(checks),'checks':checks,'observations':evidence,'exact_source_sha256':before,'fixture_scientific_identity_mount_resource_callbacks_mocked':True,'real_worker_or_fit_operations':0,'production_raw_body_reads':0,'reference_body_reads':0,'production_mutex_operations':0,'execution_admitted':False}
with (R/'genomicsem_extension_pipeline_v6_terminal_integration_receipt.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({'status':result['status'],'checks':len(checks)}))
