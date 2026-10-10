"""Partial actual84 controller handoff fixtures; zero real workers/mutex/disposals."""
import copy,hashlib,importlib.util,json,os,sys
from pathlib import Path
from types import SimpleNamespace
P=Path(__file__).resolve().parents[1];R=P/'reviews';S=P/'scripts';sys.path.insert(0,str(S));sys.dont_write_bytecode=True
ROOT=R/'core_large35_controller_handoff_controls_v3';ROOT.mkdir(exist_ok=False)
Q=Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/core_pipeline/large35_bounded_replay_v3/core_large35_bounded_replay_plan_v3.json');plan=json.loads(Q.read_text())
CHECKS=[];CASES=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fresh(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2);f.write('\n')
def check(n,c,d=None):CHECKS.append({'control':n,'pass':bool(c),'detail':d});assert c,n
def module(path):
 spec=importlib.util.spec_from_file_location('own_partial_controller84',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def controller(case):
    root=ROOT/case;root.mkdir();m=module(S/'84_run_core_large35_bounded_replay_v3.py')
    m.SSD=root/'SSD';m.OUT=m.SSD/'core_pipeline/large35_bounded_replay_v3';m.OUT.mkdir(parents=True)
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
    admission=root/'fixture_admission.json';fresh(admission,{'execution_admitted':True,'plan_sha256':ph,'executor_sha256':sha(S/'84_run_core_large35_bounded_replay_v3.py'),'independent_review_sha256':{str(dep):sha(dep)},'NOT_ACTUAL_ROOT_ADMISSION':True})
    args=SimpleNamespace(plan=planpath,plan_sha=ph,admission=admission,admission_sha=sha(admission))
    events=[];locks=[];held=[None];clock=[0];resource_fault=[None];normalwrite=m.write_new;original_generated={};original_comparator={};handoff=[]
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
            if case=='delayed_comparator_worker_swap' and where=='AFTER_CORE_COMPARE':
                path=m.OUT/'receipts_v4'/(fixture['members'][0]['trait_id']+'__compare.worker.json');v=json.loads(path.read_text());v['OWN_CHANGED_AFTER_WORKER_FREEZE']=True;path.write_text(json.dumps(v))
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
                spool=Path(member['ephemeral_spool']);spool.mkdir();(spool/'own_generated_metadata').write_text('OWN_SCRATCH_NOT_DATA\n');(spool/'columns').mkdir()
                candidate=spool/'global_harmonization_candidate_receipt.json';columns=spool/'columns/column_preparation_manifest.json'
                fresh(candidate,{'OWN_METADATA_CANDIDATE':trait});fresh(columns,{'OWN_METADATA_COLUMNS':trait})
                original_generated[trait]=(sha(candidate),sha(columns))
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
            elif label=='munge':
                if case=='delayed_candidate_swap':(Path(member['ephemeral_spool'])/'global_harmonization_candidate_receipt.json').write_text('{"OWN_REPLACED_AFTER_HARMONIZE":true}\n')
                fresh(member['munged'],{'metadata_fixture':True});outputs=[member['munged']]
            elif label=='compare':
                if case=='delayed_columns_swap':(Path(member['ephemeral_spool'])/'columns/column_preparation_manifest.json').write_text('{"OWN_REPLACED_AFTER_HARMONIZE":true}\n')
                fresh(member['comparison_receipt'],{'status':'EXACT_CORE_HARMONIZED_MUNGED_CONTENT_AND_QC_REPLAY_PASS','trait':trait,'plan_sha256':planhash});outputs=[member['comparison_receipt']]
            else:
                self.cleanup_calls+=1
                comparison=m.OUT/'receipts_v4'/(trait+'__compare.worker.json')
                expected=[original_generated[trait][0],original_generated[trait][1],original_comparator[trait]]
                flags=['--expected-candidate-sha256','--expected-columns-sha256','--expected-comparison-worker-sha256']
                actual=[command[command.index(flag)+1] for flag in flags]
                check(case+'_three_immutable_digest_arguments_exact',actual==expected)
                check(case+'_exact_dynamic_cleanup_command',command==[p['harmonization_python'],'-B',str(S/'86_cleanup_core_bounded_spool_v3.py'),'--plan',str(pp),'--plan-sha',planhash,'--trait',trait,*sum(([flag,value] for flag,value in zip(flags,expected)),[])])
                now=[sha(Path(member['ephemeral_spool'])/'global_harmonization_candidate_receipt.json'),sha(Path(member['ephemeral_spool'])/'columns/column_preparation_manifest.json'),sha(comparison)]
                check(case+'_handoff_preserves_pre_swap_digest',now==expected if case=='healthy_handoff_stop_before_disposal' else now!=expected)
                handoff.append({'command':command,'expected_frozen_digests':expected,'current_digests':now,'disposal_NOT_CALLED':True})
                fresh(record,{'status':'WORKER_FAILED_PRESERVED','command':command,'plan_sha256':planhash,'fixture_deliberate_stop_before_disposal':True})
                raise RuntimeError('OWN_FIXTURE_STOP_BEFORE_DISPOSAL_AFTER_DIGEST_HANDOFF')
            fresh(record,{'status':'WORKER_COMPLETE_VERIFIED','command':command+['WRONG'] if case=='wrong_recorded_command' and len(self.commands)==1 else command,'plan_sha256':planhash,'owned_cleanup_verified':True,'process_group_teardown':{'remaining_group_members':[]},'output_sha256':{path:sha(path) for path in [*outputs,str(journal),str(stdout)]}})
            if label=='compare':original_comparator[trait]=sha(record)
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
        if Path(path).name!='core_large35_execution_receipt_v3.json':return normalwrite(path,value)
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
    receipt=m.OUT/'core_large35_execution_receipt_v3.json';pending=m.OUT/'core_large35_pending_v3.json';sealed=m.OUT/'core_large35_terminal_seal_v3.json'
    consumed=False
    try:terminal.require_committed(pending,sealed,{'plan_sha256':ph,'admission_sha256':args.admission_sha,'executor_sha256':sha(S/'84_run_core_large35_bounded_replay_v3.py')},{str(receipt):sha(receipt)});consumed=True
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
        if case in ['healthy_handoff_stop_before_disposal','delayed_candidate_swap','delayed_columns_swap','delayed_comparator_worker_swap']:
            check(case+'_first_failed_scratch_preserved',Path(fixture['members'][0]['ephemeral_spool']).is_dir())
            check(case+'_no_second_trait_started',len(locks)==1)
    CASES.append({'handoff':handoff,'case':case,'returned_success':returned,'terminal_consumer_accepts':consumed,'error':error,'stub_commands':len(stub.commands),'symbolic_locks':len(locks),'events':events,'gate_sequence':stub.gates,'production_mutex_calls':0,'actual_worker_calls':0})
    print(json.dumps({'case_complete':case,'checks_so_far':len(CHECKS),'metadata_stub_commands':len(stub.commands)}),flush=True)

for case in ['healthy_handoff_stop_before_disposal','delayed_candidate_swap','delayed_columns_swap','delayed_comparator_worker_swap']:controller(case)
assert len(CASES)==4 and all(len(c['handoff'])==1 and c['stub_commands']==5 and c['symbolic_locks']==1 for c in CASES)
receipt={'status':'PASS_NARROW_ACTUAL84_IMMUTABLE_HANDOFF_PARTIAL_STUBS_NOT_EXECUTION_ADMISSION','source_sha256':sha(S/'84_run_core_large35_bounded_replay_v3.py'),'checks':CHECKS,'check_count':len(CHECKS),'cases':CASES,'full_controller_run_claimed':False,'helper86_disposals_called':0,'production_body_reads':0,'real_workers_mutex_or_fits':0,'execution_admitted':False}
with (R/'core_large35_controller_handoff_receipt_v3.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
print(json.dumps({'status':receipt['status'],'checks':len(CHECKS),'cases':len(CASES)}))
