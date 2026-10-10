"""Independent metadata and controller stubs; no real mutex or worker calls."""
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

P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';sys.path.insert(0,str(S));sys.dont_write_bytecode=True
SSD=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4')
PLAN=SSD/'core_pipeline/large35_bounded_replay_v1/core_large35_bounded_replay_plan_v1.json'
PIN='eb84983bef57dfb21f06a95d9f30a3871652a7b2f3cccdca30fabf9097175ed7'
ROOT=R/'core_large35_bounded_controller_controls_v1';ROOT.mkdir(exist_ok=False)
CHECKS=[];BINDINGS={};CASES=[];START=time.monotonic()
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def bind(path,wanted=None):
    path=Path(path);value=sha(path)
    assert wanted is None or value==wanted,(str(path),value,wanted)
    BINDINGS[str(path)]=value;return value
def check(name,condition,detail=None):
    CHECKS.append({'control':name,'pass':bool(condition),'detail':detail})
    if not condition:raise AssertionError(name)
def fresh(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2);f.write('\n')
def module(path):
    spec=importlib.util.spec_from_file_location('owned_controller_fixture_'+str(len(CASES)),path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
plan=json.loads(PLAN.read_text());bind(PLAN,PIN)
designpath=R/'independent_core_pipeline_replay_manifest_v4.json';bind(designpath,'451ea6d16263a467b8ac587e2069d66dfc7f736c9e694e23801efb5ba307413b')
design=json.loads(designpath.read_text());rt=json.loads(Path(plan['qualified_runtime_receipt']).read_text());compat=json.loads(Path(plan['historical_munge_compatibility_provenance']).read_text())
selected=[r for r in design['rows'] if r['prefilter'] is None and r['trait_id']!='bmi']
check('exact_original35_members_and_order',[m['original_design'] for m in plan['members']]==selected and len(selected)==plan['member_count']==35 and len(design['rows'])==45)
for m in plan['members']:
    d=m['original_design'];old=d['new_private_output_namespace'];new=str(PLAN.parent/m['trait_id'])
    def relocate(value):
        if value=='{FROZEN_HARMONIZATION_PYTHON}':return rt['command'][0]
        return new+value[len(old):] if value.startswith(old+'/') else value
    h=[relocate(x) for x in d['command_templates']['harmonize_original']]
    mu=[compat['candidate_executable'] if x==compat['pinned_unmodified_munger'] else relocate(x) for x in d['command_templates']['munge_stock']]
    check(m['trait_id']+'_complete_literal_original_argv_only_relocated',m['literal_original_harmonize_command']==h and m['munge_command']==mu and m['prefilter_command'] is None)
    bounded=[rt['command'][0],'-B',str(S/'core_bounded_harmonizer_v3.py'),'--original',h[2],'--spool',m['ephemeral_spool'],'--expected-rows',str(d['harmonizer_input_rows']),'--expected-source-sha256',d['raw']['sealed_verified_sha256'],'--',*h[3:]]
    check(m['trait_id']+'_bounded_wrapper_exact_tokens',m['harmonize_command']==bounded)
    flags=h[3:]
    check(m['trait_id']+'_exact_single_frozen_source_no_abbreviation',flags.count('--infile')==1 and flags[flags.index('--infile')+1]==d['raw']['resolved_path'] and not any(x.startswith('--infil') and x!='--infile' for x in flags))
    check(m['trait_id']+'_global_original_QC_counts_preserved',d['historical_QC']['metadata']['rows_in']==str(d['harmonizer_input_rows']) and d['historical_QC']['metadata']['rows_out']==str(d['harmonized_rows']))
unread={}
for path,wanted in plan['dependencies_sha256'].items():
    f=Path(path)
    if f.name in ['hm3_grch37_variant_map.tsv.gz','hg38ToHg19.over.chain.gz','w_hm3.snplist']:
        unread[path]={'expected_sha256':wanted,'bytes':f.stat().st_size,'body_opened':False};continue
    assert f.stat().st_size<=12*(1<<20),('unexpected large dependency',path)
    bind(path,wanted)
check('dependency_scope_184_metadata_code_3_unread_reference_bodies',len(plan['dependencies_sha256'])==187 and len(unread)==3)
for key in ['component_review_seal','bounded_adapter_review_seal']:
    sealpath=Path(plan[key]);bind(sealpath,plan[key+'_sha256']);seal=json.loads(sealpath.read_text())
    for path,item in seal['artifacts'].items():bind(path,item['sha256'])
    if key=='component_review_seal':check('reader4_sealed_component_PASS',seal['candidate_sha256']=='1e5a7f86ab51cd9cec2d705144ca7623bd8e57e218c3889cf1a94b4db7c5b375' and seal['status']=='SEALED_COMPONENT_PREPARATION_PASS_NOT_EXECUTION_ADMISSION')
    else:check('adapter3_exact_reader4_joint_binding',seal['candidate_sha256']=='3b87e849e3001a2a2d30fc0971dd3bd97847736abb9bed0dff521d1b5e1d0210' and seal['component4_sha256']=='1e5a7f86ab51cd9cec2d705144ca7623bd8e57e218c3889cf1a94b4db7c5b375')
oldplan=json.loads((SSD/'core_pipeline/original_small_input_replay_v2/core_original_small_replay_plan_v2.json').read_text());bind(SSD/'core_pipeline/original_small_input_replay_v2/core_original_small_replay_plan_v2.json','9bc40795f552ea203a1a307cdef1e6c5c5f43442407be2e1a559a22325389463')
stable=['baseline_support_package','baseline_execution_plan_sha256','baseline_dependency_sha256','baseline_relocated_dependency_sha256','baseline_input_sha256','original_190_jobs','core_precision_adjudication','python','ldsc_dir','reference_prefix','environment','guard','runtime_qualification','qualified_runtime_receipt','qualified_runtime_receipt_sha256','harmonization_python']
check('baseline190_runtime_science_and_all_guard_values_unchanged',all(plan[k]==oldplan[k] for k in stable) and len(plan['original_190_jobs'])==190)
check('scientific_attempt_and_memory_qualifications_explicit',not plan['automatic_retry'] and plan['estimator_calls']==0 and plan['scientific_membership_or_threshold_changes'] is False and 'not previously measured' in plan['scientific_adapter_qualification'] and 'source-specific attempts' in plan['scientific_adapter_qualification'])
check('full175_preprocessing_only_commands',35*5==175 and all(m['prefilter_command'] is None for m in plan['members']))
check('validator75_byte_identical_old72v2',bind(S/'75_validate_core_large35_bounded_replay_v1.py')==bind(S/'72_validate_core_original_small_replay_v2.py'))
for name in ['module','admission_gate','runtime_gate','namespace_bytes','acquire_mutex','main']:
    trees=[ast.parse((S/f).read_text()) for f in ['71_run_core_original_small_replay_v2.py','74_run_core_large35_bounded_replay_v1.py']]
    functions=[next(x for x in t.body if isinstance(x,ast.FunctionDef) and x.name==name) for t in trees]
    check('unchanged_controller_function_AST_'+name,ast.dump(functions[0],include_attributes=False)==ast.dump(functions[1],include_attributes=False))
for sealname in ['core_original_small_replay_independent_review_seal_v1.json','terminal_commit_common_review_seal_v2.json']:
    sealpath=R/sealname;bind(sealpath);seal=json.loads(sealpath.read_text())
    for path,value in seal['file_sha256'].items():bind(path,value)
for line in (R/'independent_core_small_v2_v1.sha256').read_text().splitlines():
    wanted,path=line.split('  ',1);bind(path,wanted)
bind(R/'independent_core_small_v2_v1.sha256')
prior=json.loads((R/'core_original_small_replay_controls_receipt_v1_2.json').read_text())
check('inherited41_small_controller_validator_controls_PASSED',prior['all_control_expectations_pass'] and prior['control_count']==41 and prior['actual_workers']==prior['fits']==0)
ledger=json.loads(Path(plan['global_resource_ledger_path']).read_text())
check('ledger4_4_exact_total_16GiB_allcore_300GiB',sum(ledger['component_bytes'].values())==ledger['reserved_total_bytes'] and ledger['reserved_total_bytes']+ledger['unallocated_margin_bytes']==ledger['ceiling_bytes']==300*(1<<30) and ledger['component_bytes']['proposed_core_derivative_cap_bytes']==16*(1<<30) and ledger['core16GiB_cap_includes_final_outputs_spools_indexes_prefilters_transient_and_failed_attempts'])

OUTPUTS=['harmonized','munged','source_gate_receipt','comparison_receipt','harmonization_qc','bounded_receipt','column_manifest','spool_inventory','spool_cleanup_receipt']
def controller(case):
    root=ROOT/case;root.mkdir();m=module(S/'74_run_core_large35_bounded_replay_v1.py')
    m.SSD=root/'SSD';m.OUT=m.SSD/'core_pipeline/large35_bounded_replay_v1';m.OUT.mkdir(parents=True)
    for name in ['receipts_v4','logs_v4','tmp','cache']:(m.OUT/name).mkdir()
    dep=root/'dependency_metadata.json';fresh(dep,{'fixture_only':True})
    runtimefile=root/'runtime_identity_metadata';runtimefile.write_text('not an executable\n');runtime=root/'runtime_receipt.json';fresh(runtime,{'regular_files':{str(runtimefile):{'bytes':runtimefile.stat().st_size,'sha256':sha(runtimefile)}},'symlinks':{}})
    fixture={'guard':copy.deepcopy(plan['guard']),'qualified_runtime_receipt':str(runtime),'qualified_runtime_receipt_sha256':sha(runtime),'harmonization_python':'METADATA_STUB_NEVER_EXECUTED','dependencies_sha256':{str(dep):sha(dep)},'archived_input_sha256':{},'member_count':35,'private_namespace':str(m.OUT),'members':[]}
    fixture['guard']['shared_heavy_worker_lock']=str(m.SSD/'native_heavy_worker.lock')
    for actual in plan['members']:
        trait=actual['trait_id'];td=m.OUT/trait;td.mkdir()
        for folder in ['harmonized','munged','receipts']:(td/folder).mkdir()
        raw=td/'raw_identity_metadata.json';fresh(raw,{'trait':trait,'not_a_GWAS':True})
        member={'trait_id':trait,'prefilter_command':None,'ephemeral_spool':str(td/'ephemeral_bounded_spool'),'munged_prefix':str(td/'munged'/trait),'original_design':{'raw':{'resolved_path':str(raw),'sealed_verified_sha256':sha(raw)},'prefilter':None}}
        for key in OUTPUTS:member[key]=str(td/('harmonized' if key in ['harmonized','harmonization_qc'] else 'munged' if key=='munged' else 'receipts')/(key+'.json'))
        member.update(harmonize_command=['METADATA_ONLY','harmonize',member['harmonized']],munge_command=['METADATA_ONLY','munge',member['munged']]);fixture['members'].append(member)
    planpath=root/'fixture_plan.json';fresh(planpath,fixture);ph=sha(planpath)
    admission=root/'fixture_admission.json';fresh(admission,{'execution_admitted':True,'plan_sha256':ph,'executor_sha256':sha(S/'74_run_core_large35_bounded_replay_v1.py'),'independent_review_sha256':{str(dep):sha(dep)},'NOT_ACTUAL_ROOT_ADMISSION':True})
    args=SimpleNamespace(plan=planpath,plan_sha=ph,admission=admission,admission_sha=sha(admission))
    events=[];locks=[];held=[None];clock=[0];resource_fault=[None];normalwrite=m.write_new
    class Stub:
        TERMINATION_REQUEST=[]
        def __init__(self):self.commands=[];self.gates=[];self.cleanup_calls=0
        def catchable_termination(self,*args):self.TERMINATION_REQUEST.append('metadata_signal')
        def check_termination(self,where):
            if self.TERMINATION_REQUEST:raise RuntimeError('METADATA_DEFERRED_TERMINATION')
        def baseline_gate(self,p):return {'INVENTED_190_METADATA_NOT_NATIVE_PROOF':'stable'}
        def safe_print(self,*args,**kwargs):pass
        def stage_resource_gate(self,p,start,where):
            self.gates.append(where);self.check_termination(where);state,reason=self.limits(p,start,0)
            if reason:raise RuntimeError(reason)
            return state
        def await_owned_cleanup(self,monitor,ownership,planhash):
            check(case+'_cleanup_mutex_still_symbolically_held',held[0] is not None)
            events.append('OWNED_GROUP_CLEANUP_VERIFIED');ownership[:]=[True,[]];self.cleanup_calls+=1;return []
        def worker(self,command,prefix,record,p,pp,planhash,start,monitor,ownership):
            assert held[0] is not None;self.commands.append(command)
            member=next(x for x in p['members'] if str(prefix) in [x['source_gate_receipt'],x['harmonized'],x['munged_prefix'],x['comparison_receipt'],x['spool_cleanup_receipt']]);trait=member['trait_id']
            label='source' if str(prefix)==member['source_gate_receipt'] else 'harmonize' if str(prefix)==member['harmonized'] else 'munge' if str(prefix)==member['munged_prefix'] else 'compare' if str(prefix)==member['comparison_receipt'] else 'cleanup'
            events.append(trait+':'+label+':HELD');journal=Path(str(record)+'.failure_journal.jsonl');journal.write_text('{"METADATA_STUB":true}\n');stdout=m.OUT/'logs_v4'/(record.stem+'.stdout.log');stdout.write_text('NO WORKER\n')
            if label=='harmonize':
                spool=Path(member['ephemeral_spool']);spool.mkdir();(spool/'own_generated_metadata').write_text('OWN_SCRATCH_NOT_DATA\n')
            if case in ['harmonize_failure_owned_cleanup','compare_failure','cleanup_failure'] and trait==p['members'][0]['trait_id'] and label==case.split('_')[0]:
                if case=='harmonize_failure_owned_cleanup':ownership[:]=[False,['OWNED_METADATA_NO_PROCESS']]
                fresh(record,{'status':'WORKER_FAILED_PRESERVED','command':command,'plan_sha256':planhash,'output_sha256':{str(journal):sha(journal)}});raise RuntimeError('INVENTED_'+label+'_FAILURE')
            if case=='RSS_over_2GiB' and label=='harmonize':
                state,reason=self.limits(p,start,2*(1<<30)+1);assert reason=='AGGREGATE_CORE_RSS_GUARD';raise RuntimeError(reason)
            if case=='worker_over_2h' and label=='harmonize':clock[0]=7201;state,reason=self.limits(p,start,0);assert reason=='2H_CORE_WORKER_DEADLINE';raise RuntimeError(reason)
            outputs=[]
            if label=='source':fresh(member['source_gate_receipt'],{'status':'EXACT_CORE_RAW_SOURCE_GATE_PASS','trait':trait,'plan_sha256':planhash});outputs=[member['source_gate_receipt']]
            elif label=='harmonize':
                fresh(member['harmonized'],{'metadata_fixture':True});fresh(member['harmonization_qc'],{'metadata_fixture':True});outputs=[member['harmonized']]
            elif label=='munge':fresh(member['munged'],{'metadata_fixture':True});outputs=[member['munged']]
            elif label=='compare':fresh(member['comparison_receipt'],{'status':'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS','trait':trait,'plan_sha256':planhash});outputs=[member['comparison_receipt']]
            else:
                self.cleanup_calls+=1;assert self.commands[-2][self.commands[-2].index('--mode')+1]=='compare'
                fresh(member['bounded_receipt'],{'COPIED_SYNTHETIC_CANDIDATE':True});fresh(member['column_manifest'],{'COPIED_SYNTHETIC_COLUMNS':True});fresh(member['spool_inventory'],{'OWN_GENERATED_METADATA_INVENTORY':True})
                proof={member[k]:sha(member[k]) for k in ['comparison_receipt','source_gate_receipt','harmonized','munged','harmonization_qc']}
                (Path(member['ephemeral_spool'])/'own_generated_metadata').unlink();Path(member['ephemeral_spool']).rmdir()
                fresh(member['spool_cleanup_receipt'],{'status':'OWN_GENERATED_SCRATCH_DISPOSED_AFTER_EXACT_FULL_CONTENT_QC_PASS','plan_sha256':planhash,'trait':trait,'spool_absent':True,'inventory_sha256':sha(member['spool_inventory']),'copied_candidate_receipt_sha256':sha(member['bounded_receipt']),'copied_column_manifest_sha256':sha(member['column_manifest']),'scientific_success_proof_sha256':proof});outputs=[member['spool_cleanup_receipt']]
            fresh(record,{'status':'WORKER_COMPLETE_VERIFIED','command':command+['WRONG'] if case=='wrong_recorded_command' and len(self.commands)==1 else command,'plan_sha256':planhash,'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[]},'output_sha256':{path:sha(path) for path in [*outputs,str(journal),str(stdout)]}})
    stub=Stub();m.module=lambda name,path:stub if name=='_core_owned_worker' else SimpleNamespace();m.physical_mount=lambda:{'OWN_METADATA_MOUNT':True}
    def mock_acquire(p,protected,start):
        protected.stage_resource_gate(p,start,'MOCK_WAITING_FOR_MUTEX');assert held[0] is None;held[0]=1000+len(locks);locks.append(held[0]);events.append('SYMBOLIC_LOCK_ACQUIRE');return held[0]
    m.acquire_mutex=mock_acquire
    class OSProxy:
        def close(self,fd):assert held[0]==fd;events.append('SYMBOLIC_LOCK_RELEASE');held[0]=None
        def readlink(self,path):return os.readlink(path)
    m.os=OSProxy();m.signal=SimpleNamespace(SIGINT=2,SIGTERM=15,SIGHUP=1,signal=lambda *args:None);m.time=SimpleNamespace(monotonic=lambda:clock[0])
    m.check_bindings=lambda p:check(case+'_own_dependency_identity',sha(dep)==p['dependencies_sha256'][str(dep)])
    class Disk:
        def disk_usage(self,path):return SimpleNamespace(free=0 if resource_fault[0]=='internal' and str(path)=='/System/Volumes/Data' else 1<<50)
    m.shutil=Disk();native_meter=m.namespace_bytes
    m.namespace_bytes=lambda path:17*(1<<30) if resource_fault[0]=='core_cap' and str(path).endswith('/core_pipeline') else 301*(1<<30) if resource_fault[0]=='campaign_cap' and Path(path)==m.SSD else native_meter(path)
    def save(path,value):
        if Path(path).name!='core_large35_execution_receipt_v1.json':return normalwrite(path,value)
        if case=='primary_persistence_failure':Path(path).write_text('{"PARTIAL":');raise OSError('OWN_PRIMARY_PERSISTENCE_FAULT')
        normalwrite(path,value);events.append('PRIMARY_STAGE_RECEIPT_PERSISTED')
        first=p['members'][0] if (p:=fixture) else None
        firstworker=m.OUT/'receipts_v4'/(first['trait_id']+'__source.worker.json')
        if case=='missing_old_worker':firstworker.unlink()
        if case=='missing_old_output':Path(first['harmonization_qc']).unlink()
        if case=='same_bytes_output_symlink':
            old=Path(first['harmonized']);backup=root/'same_bytes_owned_backup';backup.write_bytes(old.read_bytes());old.unlink();old.symlink_to(backup)
        if case=='same_bytes_worker_parent_symlink':
            old=m.OUT/'receipts_v4';backup=m.OUT/'same_bytes_owned_receipt_backup';old.rename(backup);old.symlink_to(backup,target_is_directory=True)
        if case=='residual_spool':Path(first['ephemeral_spool']).mkdir()
        if case=='cleanup_proof_tampered':
            path=Path(first['spool_cleanup_receipt']);v=json.loads(path.read_text());v['inventory_sha256']='0'*64;path.write_text(json.dumps(v))
        if case=='receipt_count174':value['worker_receipt_sha256'].pop(str(firstworker))
        if case=='receipt_count176':value['worker_receipt_sha256']['OWN_EXTRA_METADATA']='0'*64
        if case=='completed_order_changed':value['completed_members'].reverse()
        if case=='postpersist_internal_floor':resource_fault[0]='internal'
        if case=='postpersist_core_cap':resource_fault[0]='core_cap'
        if case=='postpersist_campaign_cap':resource_fault[0]='campaign_cap'
        if case=='postpersist_stage_deadline':clock[0]=345601
        if case=='deferred_signal_postpersist':stub.TERMINATION_REQUEST.append('METADATA_ONLY_SIGTERM')
    m.write_new=save
    returned=True;error=None
    try:m.run(args)
    except (Exception,SystemExit) as e:returned=False;error=type(e).__name__+': '+str(e)
    import terminal_commit_common_v2 as terminal
    receipt=m.OUT/'core_large35_execution_receipt_v1.json';pending=m.OUT/'core_large35_pending_v1.json';sealed=m.OUT/'core_large35_terminal_seal_v1.json'
    consumed=False
    try:terminal.require_committed(pending,sealed,{'plan_sha256':ph,'admission_sha256':args.admission_sha,'executor_sha256':sha(S/'74_run_core_large35_bounded_replay_v1.py')},{str(receipt):sha(receipt)});consumed=True
    except BaseException:pass
    healthy=case=='success'
    check(case+'_controller_and_terminal_consumer_agree',returned==consumed==healthy,{'returned':returned,'consumer_accepts':consumed,'error':error})
    check(case+'_symbolic_mutex_never_survives_verified_cleanup',held[0] is None)
    if healthy:
        check('all175_commands_35_cleanup_35_locks',len(stub.commands)==175 and stub.cleanup_calls==len(locks)==35)
        check('postpersistence_resource_gate_after_primary_receipt','POST_PERSISTENCE_CORE_TERMINAL' in stub.gates and events.index('PRIMARY_STAGE_RECEIPT_PERSISTED')>events.index('SYMBOLIC_LOCK_RELEASE'))
        for trait in [x['trait_id'] for x in fixture['members']]:check(trait+'_five_stages_lock_through_disposal',[x for x in events if x.startswith(trait+':')]==[trait+':'+label+':HELD' for label in ['source','harmonize','munge','compare','cleanup']])
    else:
        check(case+'_durable_PENDING_veto',pending.is_file())
        if case in ['harmonize_failure_owned_cleanup','compare_failure','cleanup_failure','RSS_over_2GiB','worker_over_2h']:
            check(case+'_first_failed_scratch_preserved',Path(fixture['members'][0]['ephemeral_spool']).is_dir())
            check(case+'_no_second_trait_started',len(locks)==1)
    CASES.append({'case':case,'returned_success':returned,'terminal_consumer_accepts':consumed,'error':error,'stub_commands':len(stub.commands),'symbolic_locks':len(locks),'events':events,'gate_sequence':stub.gates,'production_mutex_calls':0,'actual_worker_calls':0})
    print(json.dumps({'case_complete':case,'checks_so_far':len(CHECKS),'metadata_stub_commands':len(stub.commands)}),flush=True)

for case in ['success','harmonize_failure_owned_cleanup','compare_failure','cleanup_failure','RSS_over_2GiB','worker_over_2h','wrong_recorded_command','missing_old_worker','missing_old_output','same_bytes_output_symlink','same_bytes_worker_parent_symlink','residual_spool','cleanup_proof_tampered','receipt_count174','receipt_count176','completed_order_changed','postpersist_internal_floor','postpersist_core_cap','postpersist_campaign_cap','postpersist_stage_deadline','deferred_signal_postpersist','primary_persistence_failure']:controller(case)
AFTER={path:sha(path) for path in BINDINGS};check('all_code_plan_metadata_and_seal_identities_unchanged',BINDINGS==AFTER)
receipt={'status':'PASS_PRELAUNCH_PLAN_AND_CONTROLLER_STUBS_NOT_EXECUTION_ADMISSION','check_count':len(CHECKS),'checks':CHECKS,'case_count':len(CASES),'cases':CASES,'plan_path':str(PLAN),'plan_sha256':PIN,'source_identities_before':BINDINGS,'source_identities_after':AFTER,'unread_reference_bodies':unread,'inherited_41_controller_validator_controls_not_rerun':True,'reader4_joint_source_binding_verified':True,'memory_eligibility':'Production whole-column2GiB eligibility NOT certified; guarded source-specific attempts only.','cleanup_helper_review':'76 cleanup internals pending separate numerical review; controller stub tests model its contract without running76.','peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'elapsed_seconds':time.monotonic()-START,'production_body_reads':0,'reference_body_reads':0,'full_runtime_tree_rehashes':0,'real_mutex_operations':0,'actual_workers_or_fits':0,'execution_admitted':False}
with (R/'core_large35_bounded_controller_controls_receipt_v1.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({key:receipt[key] for key in ['status','check_count','case_count','peak_rss_bytes','elapsed_seconds']}))
