#!/usr/bin/env python3
"""Small private metadata-only controls; imports frozen helper, never launches workers."""
import hashlib
import importlib.util
import json
import os
import stat
import sys
from pathlib import Path

sys.dont_write_bytecode = True
P = Path('/Users/swethasunilkumar/Documents/Codex/2026-10-08/go/work/sleep-gwas-atlas/sleep_unified_research_v4')
HELPER = P/'scripts/terminal_commit_common_v2.py'
EXPECTED = '9ef14eda37d07b6a3442181820f8909bedd7764aae01fbcaafc966791dd5a7dd'
ROOT = P/'reviews/terminal_commit_controls_v2'
RECORDS = []


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def module():
    s = importlib.util.spec_from_file_location('_private_terminal_control', HELPER)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


def consume(m, pending, seal, binding, receipts):
    try:
        return {'accepted':True, 'seal_sha256':m.require_committed(pending,seal,binding,receipts)}
    except BaseException as e:
        return {'accepted':False,'error':type(e).__name__+': '+str(e)}


def emit(case, wanted, **data):
    good = all(data.get(k) == v for k,v in wanted.items())
    RECORDS.append({'case':case,'expected':wanted,'control_pass':good,**data})
    if not good:
        raise AssertionError(json.dumps(RECORDS[-1]))


def setup(case, binding=None):
    d=ROOT/case;d.mkdir()
    m=module();pending=d/'PENDING';seal=d/'seal.json';a=d/'receipt_a.json';b=d/'receipt_b.json'
    a.write_text('{"status":"PROVISIONAL_SUCCESS","copy":"a"}\n')
    b.write_text('{"status":"PROVISIONAL_SUCCESS","copy":"b"}\n')
    binding={'stage':'metadata_only_fixture','plan_sha256':'a'*64} if binding is None else binding
    return d,m,pending,seal,a,b,binding


def suppress_supplemental(m):
    real=m.save_new
    def save(path,value):
        if str(path).endswith('.failure.json'):
            raise OSError('INJECTED_SUPPLEMENTAL_WRITE_FAILURE')
        return real(path,value)
    m.save_new=save


def commit_case(case, mode=None, suppress=False):
    d,m,pending,seal,a,b,binding=setup(case)
    if mode in {'unlink_before_failure','unlink_after_effect_failure'}:
        Base=type(Path())
        class PrivatePath(Base):
            def unlink(self,*args,**kwargs):
                if str(self)==str(pending):
                    if mode=='unlink_after_effect_failure':
                        super().unlink(*args,**kwargs)
                    raise OSError('INJECTED_MARKER_UNLINK_FAILURE')
                return super().unlink(*args,**kwargs)
        m.Path=PrivatePath
    t=m.TerminalCommit(pending,seal,binding)
    receipts={str(a):digest(a),str(b):digest(b)}
    if suppress:suppress_supplemental(m)
    if mode=='empty_receipts':receipts={}
    if mode=='wrong_receipt_hash':receipts[str(a)]='b'*64
    if mode=='missing_receipt':a.unlink()
    if mode=='receipt_symlink':
        clone=d/'receipt_clone.json';clone.write_bytes(a.read_bytes());a.unlink();a.symlink_to(clone)
    if mode=='pending_changed':pending.write_text('CHANGED\n')
    if mode=='pending_symlink':
        clone=d/'pending_clone';clone.write_bytes(pending.read_bytes());pending.unlink();pending.symlink_to(clone)
    if mode=='pending_deleted_before':pending.unlink()
    if mode in {'receipt_failure_regular','seal_failure_regular','receipt_failure_dangling','seal_failure_dangling'}:
        target=Path(str(a if mode.startswith('receipt') else seal)+'.failure.json')
        if mode.endswith('dangling'):target.symlink_to(d/'absent_failure_target')
        else:target.write_text('{"status":"FAILED"}\n')
    if mode in {'seal_write_failure','seal_partial_write_failure'}:
        old=m.save_new
        def save(path,value):
            if str(path)==str(seal):
                if mode=='seal_partial_write_failure':Path(path).write_text('{"partial":')
                raise OSError('INJECTED_SEAL_WRITE_FAILURE')
            return old(path,value)
        m.save_new=save
    if mode=='seal_directory_sync_failure':
        old=m.sync_directory
        m.sync_directory=lambda path:(_ for _ in ()).throw(OSError('INJECTED_DIRECTORY_SYNC_FAILURE'))
    if mode=='seal_file_fsync_failure':
        class PrivateOS:
            def __getattr__(self,k):return getattr(os,k)
            def fsync(self,fd):
                if stat.S_ISREG(os.fstat(fd).st_mode):raise OSError('INJECTED_SEAL_FILE_FSYNC_FAILURE')
                return os.fsync(fd)
        m.os=PrivateOS()
    if mode in {'seal_hash_readback_failure','receipt_hash_readback_failure'}:
        old=m.sha;counts={}
        def hashfile(path):
            key=str(path);counts[key]=counts.get(key,0)+1
            if mode=='seal_hash_readback_failure' and key==str(seal):raise OSError('INJECTED_SEAL_HASH_READ_FAILURE')
            if mode=='receipt_hash_readback_failure' and key==str(a) and counts[key]>=2:raise OSError('INJECTED_RECEIPT_HASH_READ_FAILURE')
            return old(path)
        m.sha=hashfile
    call={'identity':0,'resource':0,'termination':0}
    def identity():
        call['identity']+=1
        n=call['identity']
        if mode=='identity_first_failure' and n==1:raise RuntimeError('INJECTED_IDENTITY_FAILURE')
        if mode=='identity_post_failure' and n==2:raise RuntimeError('INJECTED_POST_IDENTITY_FAILURE')
        if n==2:
            if mode=='receipt_post_changed':a.write_text('CHANGED_AFTER_PERSISTENCE\n')
            if mode=='seal_post_changed':seal.write_text('CHANGED_AFTER_PERSISTENCE\n')
            if mode=='pending_post_changed':pending.write_text('CHANGED_AFTER_PERSISTENCE\n')
            if mode=='pending_deleted_after_persistence':pending.unlink()
            if mode in {'receipt_post_symlink','seal_post_symlink','pending_post_symlink'}:
                target={'receipt_post_symlink':a,'seal_post_symlink':seal,'pending_post_symlink':pending}[mode]
                clone=d/'post_symlink_clone';clone.write_bytes(target.read_bytes());target.unlink();target.symlink_to(clone)
    def resource():
        call['resource']+=1
        if mode=='resource_first_failure' and call['resource']==1:raise RuntimeError('INJECTED_RESOURCE_FAILURE')
        if mode=='resource_post_failure' and call['resource']==2:raise RuntimeError('INJECTED_POST_PERSISTENCE_RESOURCE_FAILURE')
    def termination():
        call['termination']+=1
        if mode=='signal_first_failure' and call['termination']==1:raise KeyboardInterrupt('INJECTED_DEFERRED_SIGNAL')
        if mode=='signal_post_failure' and call['termination']==2:raise KeyboardInterrupt('INJECTED_POST_PERSISTENCE_DEFERRED_SIGNAL')
    committed=t.commit(receipts,identity_gate=identity,resource_gate=resource,termination_gate=termination)
    # Restore only the private module's instrumented interfaces for independent consumption.
    reader=module()
    c=consume(reader,pending,seal,binding,receipts)
    expected={'commit_success':mode is None,'consumer_accepts':mode is None}
    # These witnesses are outside the documented immutable/deferred ownership boundary,
    # but record the exact behavior rather than inventing protection the helper lacks.
    limitation=mode in {'unlink_after_effect_failure','pending_deleted_after_persistence'} and suppress
    if limitation:expected={'commit_success':False,'consumer_accepts':True}
    emit(case,expected,commit_success=committed,consumer_accepts=c['accepted'],
         consumer=c,pending_present_or_symlink=pending.exists() or pending.is_symlink(),
         seal_present_or_symlink=seal.exists() or seal.is_symlink(),callback_counts=call,
         supplemental_write_failure_injected=suppress,
         documented_ownership_boundary_violation=limitation)


def initial_case(case,mode):
    d,m,pending,seal,a,b,binding=setup(case,{} if mode=='empty_stage_binding' else None)
    if mode=='initial_pending_write_failure':
        m.save_new=lambda *args:(_ for _ in ()).throw(OSError('INJECTED_INITIAL_WRITE_FAILURE'))
    if mode=='initial_directory_sync_failure':
        m.sync_directory=lambda *args:(_ for _ in ()).throw(OSError('INJECTED_INITIAL_SYNC_FAILURE'))
    if mode=='initial_hash_read_failure':
        m.sha=lambda *args:(_ for _ in ()).throw(OSError('INJECTED_INITIAL_HASH_FAILURE'))
    if mode=='existing_pending':pending.write_text('PRIOR_ATTEMPT\n')
    if mode=='existing_seal':seal.write_text('PRIOR_ATTEMPT\n')
    if mode=='existing_pending_dangling':pending.symlink_to(d/'missing')
    if mode=='existing_seal_dangling':seal.symlink_to(d/'missing')
    try:m.TerminalCommit(pending,seal,binding);created=True;error=None
    except BaseException as e:created=False;error=type(e).__name__+': '+str(e)
    c=consume(module(),pending,seal,binding,{str(a):digest(a),str(b):digest(b)})
    emit(case,{'created':False,'consumer_accepts':False},created=created,consumer_accepts=c['accepted'],error=error,
         pending_present_or_symlink=pending.exists() or pending.is_symlink())


def consumer_case(case,mode):
    d,m,pending,seal,a,b,binding=setup(case)
    t=m.TerminalCommit(pending,seal,binding);receipts={str(a):digest(a),str(b):digest(b)}
    assert t.commit(receipts,identity_gate=lambda:None,resource_gate=lambda:None,termination_gate=lambda:None)
    if mode=='receipt_mutated':a.write_text('POSTCOMMIT_MUTATION\n')
    if mode in {'receipt_symlink','seal_symlink'}:
        target=a if mode.startswith('receipt') else seal;clone=d/'clone';clone.write_bytes(target.read_bytes());target.unlink();target.symlink_to(clone)
    if mode=='pending_dangling':pending.symlink_to(d/'missing')
    if mode=='pending_regular':pending.write_text('REAPPEARED_PENDING\n')
    if mode in {'receipt_failure_regular','receipt_failure_dangling','seal_failure_regular','seal_failure_dangling'}:
        target=Path(str(a if mode.startswith('receipt') else seal)+'.failure.json')
        if mode.endswith('dangling'):target.symlink_to(d/'missing')
        else:target.write_text('{"status":"FAILED"}\n')
    if mode in {'seal_status','seal_binding','seal_receipts','seal_pending_path'}:
        x=json.loads(seal.read_text())
        key={'seal_status':'status','seal_binding':'binding','seal_receipts':'result_receipt_sha256','seal_pending_path':'pending_path'}[mode]
        x[key]='MUTATED';seal.write_text(json.dumps(x)+'\n')
    if mode=='binding_empty':binding={}
    c=consume(module(),pending,seal,binding,receipts)
    emit(case,{'consumer_accepts':False},consumer_accepts=c['accepted'],consumer=c)


def main():
    assert digest(HELPER)==EXPECTED
    ROOT.mkdir(exist_ok=False)
    commit_case('successful_commit_and_consumer')
    for mode in ['empty_receipts','wrong_receipt_hash','missing_receipt','receipt_symlink','pending_changed','pending_symlink',
                 'pending_deleted_before','receipt_failure_regular','seal_failure_regular','receipt_failure_dangling','seal_failure_dangling',
                 'seal_write_failure','seal_partial_write_failure','seal_directory_sync_failure','seal_file_fsync_failure',
                 'seal_hash_readback_failure','receipt_hash_readback_failure','identity_first_failure','identity_post_failure',
                 'resource_first_failure','resource_post_failure','signal_first_failure','signal_post_failure','receipt_post_changed',
                 'seal_post_changed','pending_post_changed','receipt_post_symlink','seal_post_symlink','pending_post_symlink','unlink_before_failure']:
        commit_case(mode,mode,suppress=mode in {'seal_write_failure','seal_partial_write_failure','seal_directory_sync_failure',
                   'seal_file_fsync_failure','seal_hash_readback_failure','receipt_hash_readback_failure','identity_first_failure',
                   'identity_post_failure','resource_first_failure','resource_post_failure','signal_first_failure','signal_post_failure',
                   'receipt_post_changed','seal_post_changed','pending_post_changed','receipt_post_symlink','seal_post_symlink',
                   'pending_post_symlink','unlink_before_failure'})
    for mode in ['empty_stage_binding','initial_pending_write_failure','initial_directory_sync_failure','initial_hash_read_failure',
                 'existing_pending','existing_seal','existing_pending_dangling','existing_seal_dangling']:
        initial_case(mode,mode)
    for mode in ['receipt_mutated','receipt_symlink','seal_symlink','pending_dangling','pending_regular',
                 'receipt_failure_regular','receipt_failure_dangling','seal_failure_regular','seal_failure_dangling',
                 'seal_status','seal_binding','seal_receipts','seal_pending_path','binding_empty']:
        consumer_case('consumer_'+mode,mode)
    # Explicit boundary witnesses: a caller must own/protect PENDING and defer asynchronous raises.
    commit_case('boundary_unlink_after_effect_failure','unlink_after_effect_failure',suppress=True)
    commit_case('boundary_pending_deleted_after_persistence','pending_deleted_after_persistence',suppress=True)
    # A failed second-copy write occurs before commit, so initial PENDING remains the veto.
    d,m,pending,seal,a,b,binding=setup('second_copy_write_failure')
    t=m.TerminalCommit(pending,seal,binding);b.write_text('{"partial":')
    c=consume(module(),pending,seal,binding,{str(a):digest(a)})
    emit('second_copy_write_failure',{'consumer_accepts':False,'pending_present':True},consumer_accepts=c['accepted'],pending_present=pending.is_file())
    # Recheck the successful fixture's exact durable pending metadata evidence before deletion
    # through a separate initialization-only instance, with real file/directory fsync.
    d,m,pending,seal,a,b,binding=setup('initial_pending_durability')
    events=[];save=m.save_new;sync=m.sync_directory
    def wrapped_save(path,v):events.append('save_new');return save(path,v)
    def wrapped_sync(path):events.append('sync_directory');return sync(path)
    m.save_new=wrapped_save;m.sync_directory=wrapped_sync;t=m.TerminalCommit(pending,seal,binding)
    q=json.loads(pending.read_text())
    emit('initial_pending_durability',{'durable_metadata_match':True},durable_metadata_match=
         events==['save_new','sync_directory'] and t.pending_sha==digest(pending) and q['binding']==binding
         and q['status']=='PENDING_NOT_ADMISSIBLE' and q['owner_pid']==os.getpid(),events=events)
    assert digest(HELPER)==EXPECTED
    files={};links={}
    for f in sorted(ROOT.rglob('*')):
        if f.is_symlink():links[str(f)]=os.readlink(f)
        elif f.is_file():files[str(f)]=digest(f)
    receipt={'schema':'independent_terminal_commit_v2_metadata_controls','helper_sha256':EXPECTED,
             'control_script_sha256':digest(__file__),'controls':RECORDS,'control_count':len(RECORDS),
             'all_control_expectations_pass':all(x['control_pass'] for x in RECORDS),
             'qualified_boundary_false_admission_witnesses':[x['case'] for x in RECORDS if x.get('documented_ownership_boundary_violation')],
             'GWAS_body_reads':0,'reference_body_reads':0,'decompressions':0,'workers':0,'fits':0,
             'fixture_regular_file_sha256':files,'fixture_symlink_literal_targets':links,
             'execution_admission_granted':False}
    out=P/'reviews/terminal_commit_controls_receipt_v2.json'
    with out.open('x') as f:json.dump(receipt,f,indent=2,sort_keys=True);f.write('\n')
    print(json.dumps({'controls':len(RECORDS),'receipt_sha256':digest(out),
                      'boundary_witnesses':receipt['qualified_boundary_false_admission_witnesses'],
                      'helper_sha256_after':digest(HELPER)},indent=2))


if __name__=='__main__':main()
