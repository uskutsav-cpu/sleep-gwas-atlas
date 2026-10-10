"""One guarded, separately admitted source-only attempt. Preparation is elsewhere."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import time

from signed_reference_source_common_v3 import (PLAN,ADMISSION,load_frozen,save,Guard,
    identity_gate,resource_snapshot,exclusive_heavy_lock,shared_raw_family_lock,
    deferred_termination_signals,assert_no_termination,regular_sha,verify_outputs,
    freeze_registered,freeze_companions,safe_diagnostic,terminal_intended_sha)
from canonical_calibration_common_v4_8 import assert_monitor_error_fields_clear
from terminal_commit_common_v2 import TerminalCommit
from signed_reference_compact_v3 import SCHEMA, SCHEMA_SHA256, static_rows, genotype_rows


def monitor_helper(plan):
    spec=importlib.util.spec_from_file_location('source_owned_helpers',plan['owned_monitor_helper'])
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def group_exists(pid):
    try:os.killpg(pid,0);return True
    except ProcessLookupError:return False


def cleanup(plan,ownership):
    while ownership:
        entry=ownership[0];proc=entry['proc']
        try:
            helper=monitor_helper(plan)
            proof=helper.terminate_owned(proc)
            if group_exists(proc.pid):raise RuntimeError('OWNED_SOURCE_GROUP_STILL_PRESENT')
        except BaseException as error:
            entry['cleanup_errors'].append(repr(error))
            try:
                if group_exists(proc.pid):os.killpg(proc.pid,signal.SIGKILL)
                proc.wait(timeout=10)
                if group_exists(proc.pid):raise RuntimeError('OWNED_SOURCE_GROUP_AFTER_KILL')
                proof=dict(helper_error=repr(error),remaining_group_members=[],signal_zero_group_absent=True)
            except BaseException as recovery:
                entry['cleanup_errors'].append(repr(recovery))
                safe_diagnostic(lambda:'SOURCE_BOTH_LOCKS_RETAINED_UNVERIFIED_GROUP: '+repr(recovery))
                try:time.sleep(1)
                except BaseException:pass
                continue
        entry['final_cleanup']=proof;ownership.pop(0)


def runtime_profile_gate(plan,admission):
    # The existing A identity was sealed before this work; full current reads
    # occur only in the separately admitted stage. B's complete current runtime
    # audit is an explicit independent prerequisite, not fabricated at preparation.
    profiles=[(plan['runtime_A_profile_path'],plan['runtime_A_profile_sha256'],plan['runtime_A_environment_root'])]
    bp=plan['runtime_B_required_profile_path'];bh=admission.get('runtime_B_profile_sha256')
    if admission.get('independent_complete_B_runtime_profile_review_pass') is not True or bh!=plan['runtime_B_profile_sha256']:raise RuntimeError('SOURCE_DECODER_B_COMPLETE_RUNTIME_IDENTITY_PENDING')
    profiles.append((bp,bh,plan['runtime_B_environment_root']))
    for path,digest,expected_root in profiles:
        profile=load_frozen(path,digest)
        if profile['environment_root']!=expected_root:raise RuntimeError('DECLARED_READ_ONLY_RUNTIME_ROOT_DIFFERS')
        root=Path(profile['environment_root']);directory_links=profile.get('symlink_directories',{})
        expected=set(profile['regular_files'])|set(profile['symlinks'])|set(directory_links)
        observed={str(p) for p in root.rglob('*') if p.is_symlink() or p.is_file()}
        if observed!=expected:raise RuntimeError('READ_ONLY_RUNTIME_EXACT_INVENTORY_CHANGED')
        if path==bp:
            physical_dirs=sorted([str(root),*[str(p) for p in root.rglob('*') if p.is_dir() and not p.is_symlink()]])
            if profile.get('physical_directories')!=physical_dirs:raise RuntimeError('B_RUNTIME_COMPLETE_PHYSICAL_DIRECTORY_CENSUS_REQUIRED')
        for p,item in profile['regular_files'].items():
            if Path(p).stat().st_size!=item['bytes'] or regular_sha(p)!=item['sha256']:raise RuntimeError('READ_ONLY_RUNTIME_REGULAR_FILE_CHANGED')
        for p,item in profile['symlinks'].items():
            p=Path(p)
            if not p.is_symlink() or os.readlink(p)!=item['literal_link'] or str(p.resolve())!=item['resolved_target']:raise RuntimeError('DECLARED_RUNTIME_SYMLINK_CHANGED')
            if item.get('target_is_file') is False:
                # Legacy independently sealed A profile recorded contained
                # directory links without a separate directory census. Its
                # complete physical regular-file/symlink closure is retained;
                # this is not the new B directory/tree identity qualification.
                if path==bp or not p.resolve().is_dir() or not p.resolve().is_relative_to(root):raise RuntimeError('LEGACY_A_DIRECTORY_LINK_CONTRACT_CHANGED')
            elif regular_sha(p.resolve())!=item['resolved_sha256']:raise RuntimeError('DECLARED_RUNTIME_SYMLINK_TARGET_CHANGED')
        for p,item in directory_links.items():
            p=Path(p);target=p.resolve()
            if not p.is_symlink() or os.readlink(p)!=item['literal_link'] or str(target)!=item['resolved_target'] or not target.is_dir() or target.is_symlink() or not target.is_relative_to(root):raise RuntimeError('B_RUNTIME_DIRECTORY_LINK_NOT_EXACT_CONTAINED')
            def below(text):return Path(text)==target or Path(text).is_relative_to(target)
            tree=dict(regular_files={k:{'bytes':v['bytes'],'sha256':v['sha256']} for k,v in profile['regular_files'].items() if below(k)},
                physical_directories=[k for k in profile['physical_directories'] if below(k)],
                file_symlinks={k:{'literal_link':v['literal_link'],'resolved_target':v['resolved_target'],'resolved_sha256':v['resolved_sha256']} for k,v in profile['symlinks'].items() if below(k)},
                directory_symlinks={k:{'literal_link':v['literal_link'],'resolved_target':v['resolved_target']} for k,v in directory_links.items() if below(k)})
            th=hashlib.sha256(json.dumps(tree,sort_keys=True,separators=(',',':')).encode()).hexdigest()
            if item.get('tree_sha256')!=th:raise RuntimeError('B_RUNTIME_DIRECTORY_LINK_TREE_IDENTITY_FAILED')
        if regular_sha(path)!=digest:raise RuntimeError('READ_ONLY_RUNTIME_PROFILE_CHANGED_AFTER_USE')


def monitor_worker(plan,command,fds,started,ownership,gates):
    root=Path(plan['output_root']);helper=monitor_helper(plan);entry=None;proc=None;first=None;reason=None;error=None;peak=0
    monitor=dict(status='FAILED_PRESERVED',returncode=None,stop_reason=None,cleanup_errors=[],run_error=None,
      process_group_teardown=None,first_teardown_attempt=None,command=command,plan_sha256=gates['plan_sha'])
    try:
        gates['identity']();resource_snapshot(plan,started,launch=True)
        env={**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1',
          'MKL_NUM_THREADS':'1','VECLIB_MAXIMUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1',
          'TMPDIR':str(root/'tmp'),'XDG_CACHE_HOME':str(root/'tmp'),'MPLCONFIGDIR':str(root/'tmp')}
        with (root/'worker_stdout.log').open('xb') as log:
            proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True,pass_fds=fds)
            entry=dict(proc=proc,cleanup_errors=[]);ownership.append(entry)
            while True:
                rc=proc.poll();rss,pids=helper.owned_rss(proc.pid);peak=max(peak,rss)
                resource_snapshot(plan,started)
                if rss>plan['resource_contract']['aggregate_RSS_stop_bytes']:raise RuntimeError('AGGREGATE_SOURCE_OWNED2GIB_RSS')
                if rc is not None:break
                time.sleep(1)
            monitor['returncode']=rc
    except BaseException as e:error=dict(type=type(e).__name__,repr=repr(e));reason='SOURCE_MONITOR_OR_WORKER_FAILED'
    finally:
        if proc is not None:
            if entry is None:entry=dict(proc=proc,cleanup_errors=[]);ownership.append(entry)
            try:
                first=helper.terminate_owned(proc)
                if first.get('initial_group_members'):reason=reason or 'UNEXPECTED_DESCENDANTS_AT_SOURCE_TERMINATION'
                if group_exists(proc.pid):raise RuntimeError('SOURCE_GROUP_EXISTS_AFTER_INITIAL_TEARDOWN')
            except BaseException as e:
                first=dict(exception=repr(e));entry['cleanup_errors'].append(repr(e));reason=reason or 'SOURCE_TEARDOWN_EXCEPTION'
            cleanup(plan,ownership)
            monitor.update(returncode=proc.returncode,owned_process_group_id=proc.pid,
              reaped_immediate_child=True,process_group_teardown=entry.get('final_cleanup'),
              first_teardown_attempt=first,cleanup_errors=list(entry['cleanup_errors']),signal_zero_group_absent=True)
        monitor.update(peak_aggregate_RSS_bytes=peak,run_error=error,stop_reason=reason)
        if monitor['returncode']!=0:monitor['stop_reason']=monitor['stop_reason'] or 'SOURCE_WORKER_NONZERO'
        try:gates['identity']();resource_snapshot(plan,started)
        except BaseException as e:
            monitor['post_teardown_error']=repr(e);monitor['stop_reason']=monitor['stop_reason'] or 'SOURCE_POST_TEARDOWN_GATE'
        if not monitor['stop_reason']:monitor['status']='COMPLETE_REAPED_SOURCE_ONLY'
        path=root/'monitor_receipt.json'
        try:save(path,monitor)
        except BaseException as e:
            reason=reason or 'SOURCE_MONITOR_PERSISTENCE_FAILED'
            try:save(str(path)+'.failure.json',dict(status='FAILED_PRESERVED',first_stop=reason,error=repr(e)))
            except BaseException as fallback:safe_diagnostic(lambda:'SOURCE_MONITOR_AND_FAILURE_WRITE_FAILED: '+repr(fallback))
            raise
    if monitor['stop_reason']:raise RuntimeError(monitor['stop_reason'])
    assert_monitor_error_fields_clear(monitor)
    return path,regular_sha(path)


def expected_worker_outputs(root):
    root=Path(root);q=root/'qc';c=q/'controls'
    names=[root/'archive/1000G_Phase3_plinkfiles.tgz',root/'archive/1000G_Phase3_plinkfiles.tgz.curl_status',
      q/'fresh_zenodo_metadata.json',q/'fresh_zenodo_metadata.json.curl_status',q/'fresh_archive_HEAD.txt',q/'fresh_archive_HEAD.txt.curl_status',
      q/'G2_authentication.json',q/'archive_members.json',q/'private_sample_order.json',q/'G3_format.json',root/'extracted/source_index.sqlite',q/'reference_index.sqlite',
      q/'compact_schema.json',q/'master_static.bin',q/'master_SNP_utf8.bin',q/'G4_static_seal.json',q/'genotype_QC.bin',q/'G4_J_seal.json',
      c/'fixed_G5_selections.json',q/'G5_operator_controls.json']
    names.extend(root/'extracted/members'/f'1000G.EUR.QC.{chrom}.{ext}' for chrom in range(1,23) for ext in ['bed','bim','fam'])
    names.extend(c/f'chr{chrom}.{ext}' for chrom in range(1,23) for ext in ['bed','metadata.json','B.npz','B.receipt.json','B.log'])
    return {str(p) for p in names}


def consume_worker_result(plan,monitor_path,monitor_sha,plan_sha,expected_command,immutable_result_sha=None,resource_guard=lambda:None):
    resource_guard()
    monitor=load_frozen(monitor_path,monitor_sha);assert_monitor_error_fields_clear(monitor)
    if monitor.get('command')!=list(expected_command):raise RuntimeError('SOURCE_EXACT_CONTROLLER_DERIVED_COMMAND_REQUIRED')
    if monitor.get('status')!='COMPLETE_REAPED_SOURCE_ONLY' or monitor.get('returncode')!=0 or monitor.get('stop_reason') is not None or monitor.get('reaped_immediate_child') is not True or monitor.get('signal_zero_group_absent') is not True or monitor.get('plan_sha256')!=plan_sha or monitor.get('process_group_teardown',{}).get('remaining_group_members')!=[]:raise RuntimeError('SOURCE_MONITOR_NOT_ADMISSIBLE')
    root=Path(plan['output_root']);path=root/'qc/source_worker_result.json';digest=immutable_result_sha or regular_sha(path);result=load_frozen(path,digest)
    if result.get('status')!='SOURCE_G2_G5_PASS_PRIVATE_PROVISIONAL' or result.get('plan_sha256')!=plan_sha or not result.get('ordinary_output_sha256'):raise RuntimeError('SOURCE_WORKER_RESULT_BINDING_FAILED')
    for field in ['S_simulations_run','C_calibration_credit','estimator_calls']:
        if result.get(field)!=0:raise RuntimeError('SOURCE_SCOPE_CONTRADICTION')
    if result.get('biological41Cov_release_admitted') is not False:raise RuntimeError('SOURCE_SCOPE_CONTRADICTION')
    mapping=result['ordinary_output_sha256']
    if set(mapping)!=expected_worker_outputs(root):raise RuntimeError('SOURCE_EXACT_REGISTERED_G2_G5_OUTPUT_SET_REQUIRED')
    for p in mapping:
        if Path(p).name.startswith('._') or not Path(p).is_relative_to(root):raise RuntimeError('SOURCE_RESULT_OUTPUT_ROUTE_INVALID')
    verify_outputs(mapping)
    def phase(name):
        p=root/'qc'/name;return load_frozen(p,mapping[str(p)])
    g2=phase('G2_authentication.json');g3=phase('G3_format.json');g4s=phase('G4_static_seal.json');g4j=phase('G4_J_seal.json');g5=phase('G5_operator_controls.json')
    archive=root/'archive/1000G_Phase3_plinkfiles.tgz';ai=result.get('archive_identity',{})
    if g2.get('status')!='AUTHENTICATED_NEW_CANDIDATE_NOT_OLD_COMPACT_IDENTITY' or g2.get('candidate')!=plan['candidate'] or g2.get('archive_identity')!=ai or ai.get('bytes')!=plan['candidate']['bytes'] or 'md5:'+str(ai.get('md5'))!=plan['candidate']['upstream_md5'] or ai.get('sha256')!=mapping[str(archive)]:raise RuntimeError('SOURCE_UPSTREAM_ARCHIVE_AUTHENTICATION_RECEIPT_CONTRADICTION')
    if g3.get('status')!='EXACT_SOURCE_FORMAT_PASS' or g3.get('samples')!=503 or g3.get('all22_fam_byte_identical') is not True or g3.get('padding_all_rows_checked') is not True or set(g3.get('chromosome_variant_counts',{}))!={str(c) for c in range(1,23)} or g3.get('source_index_sha256')!=result.get('source_index_sha256') or result.get('source_index_sha256')!=mapping[str(root/'extracted/source_index.sqlite')]:raise RuntimeError('SOURCE503_FORMAT_RECEIPT_CONTRADICTION')
    sb=g4s.get('compact_static',{});qb=g4j.get('compact_genotype',{})
    if phase('compact_schema.json')!=SCHEMA or g4s.get('compact_schema_file_sha256')!=mapping[str(root/'qc/compact_schema.json')]:raise RuntimeError('SOURCE_EXACT_COMPACT_SCHEMA_FILE_REQUIRED')
    if (g4s.get('status')!='FROZEN_BEFORE_GENOTYPE_QC' or g4s.get('master_reference_rows')!=plan['expected_master_reference_rows'] or
        sb.get('schema_sha256')!=SCHEMA_SHA256 or sb.get('master_rows')!=plan['expected_master_reference_rows'] or
        sb.get('master_sha256')!=mapping[str(root/'qc/master_static.bin')] or sb.get('SNP_pool_sha256')!=mapping[str(root/'qc/master_SNP_utf8.bin')] or
        sb.get('static_rows')!=g4s.get('static_rows') or sb.get('logical_sha256',{}).get('static_selection')!=g4s.get('static_ordered_sha256') or
        sb.get('logical_sha256',{}).get('static_exclusions')!=g4s.get('static_exclusions_sha256')):raise RuntimeError('SOURCE_STATIC_BEFORE_QC_FROZEN_RECEIPT_CONTRADICTION')
    if qb.get('schema_sha256')!=SCHEMA_SHA256 or qb.get('genotype_QC_sha256')!=mapping[str(root/'qc/genotype_QC.bin')] or qb.get('QC_rows')!=sb.get('static_rows') or qb.get('J_rows')!=g4j.get('J_count'):raise RuntimeError('SOURCE_COMPLETE_COMPACT_QC_BINDING_REQUIRED')
    # Reconstruct all original logical streams, including exclusions, with exact
    # integer/float checks; a physical hash alone is not a logical equality proof.
    from collections import Counter
    actual_static_exclusions=Counter(r['exclusion'] for r in static_rows(root/'qc',sb,resource_guard,selected=False))
    actual_genotype_exclusions=Counter(r['exclusion'] for r in genotype_rows(root/'qc',sb,qb,resource_guard,eligible=False))
    actual_counts=Counter(str(r['CHR']) for r in genotype_rows(root/'qc',sb,qb,resource_guard))
    if dict(actual_static_exclusions)!=g4s.get('exclusion_counts') or dict(actual_genotype_exclusions)!=g4j.get('genotype_exclusion_counts'):raise RuntimeError('SOURCE_COMPLETE_EXCLUSION_RECONSTRUCTION_DIFFERS')
    counts=result.get('chromosome_J_counts',{})
    if dict(actual_counts)!=counts:raise RuntimeError('SOURCE_RECONSTRUCTED_EXACT_J_COUNTS_DIFFER')
    if set(counts)!={str(c) for c in range(1,23)} or any(type(n) is not int or n<192 for n in counts.values()):raise RuntimeError('SOURCE_EXACT22_NONEMPTY_G5_UNIVERSES_REQUIRED')
    if g4j.get('status')!='REFERENCE_ONLY_J_FROZEN' or g4j.get('J_ordered_sha256')!=result.get('J_ordered_sha256') or result.get('J_ordered_sha256')!=qb.get('logical_sha256',{}).get('J_ordered') or g4j.get('genotype_exclusions_sha256')!=qb.get('logical_sha256',{}).get('genotype_exclusions') or g4j.get('chromosome_J_counts')!=counts or g4j.get('J_count')!=sum(counts.values()) or g4j.get('no_MAF_threshold') is not True or g4j.get('no_GWAS_outcomes') is not True or g4j.get('no_native_mask_selection') is not True:raise RuntimeError('SOURCE_REFERENCE_ONLY_J_RECEIPT_CONTRADICTION')
    if g5.get('status')!='ACTUAL_AUTHENTICATED_SOURCE_CONTROLS_PASS' or [r.get('CHR') for r in g5.get('chromosomes',[])]!=list(range(1,23)) or g5.get('fixed_selection_sha256')!=mapping[str(root/'qc/controls/fixed_G5_selections.json')] or g5.get('maintained_B_padding_independent_rejection') is not False or g5.get('zero_variance_fallback_not_credited') is not True or g5.get('finite_reference_503_conditional_only') is not True or g5.get('population_LD_or_actual_GWAS_mask_adequacy_established') is not False:raise RuntimeError('SOURCE_FINITE_REFERENCE_G5_RECEIPT_CONTRADICTION')
    for r in g5['chromosomes']:
        c=r['CHR'];pref=root/'qc/controls'/f'chr{c}'
        if r.get('J_count')!=counts[str(c)] or not 192<=r.get('independent_rows',0)<=288 or r.get('source_window_sha256')!=mapping[str(pref)+'.bed'] or r.get('B_array_sha256')!=mapping[str(pref)+'.B.npz'] or r.get('B_receipt_sha256')!=mapping[str(pref)+'.B.receipt.json'] or r.get('controls')!=['first64','middle64','last64','early32_scatter','middle32_scatter','late32_scatter','combined_distant_and_windows'] or r.get('full_chromosome_energy_symmetry') is not True or r.get('between_storage_window_correlations_retained') is not True:raise RuntimeError('SOURCE_G5_NAMED_CONTROL_IDENTITY_OR_SCOPE_CONTRADICTION')
    for value in [g5,result]:
        if any(value.get(k)!=0 for k in ['S_simulations_run','C_calibration_credit','estimator_calls']) or value.get('biological41Cov_release_admitted') is not False:raise RuntimeError('SOURCE_G5_PROHIBITED_SCOPE_CREDIT')
    verify_outputs(mapping)
    if regular_sha(path)!=digest:raise RuntimeError('SOURCE_RESULT_CHANGED_AFTER_PHASE_PARSE')
    mapping=dict(mapping);mapping.update({str(path):digest,str(monitor_path):monitor_sha,str(root/'worker_stdout.log'):regular_sha(root/'worker_stdout.log')})
    targets=[*[root/name for name in ['', 'archive', 'extracted', 'extracted/members', 'qc', 'qc/controls', 'tmp']],*[Path(p) for p in mapping]]
    return mapping,targets,result


def execute(plan_path,plan_sha,admission_path,admission_sha):
    started=time.monotonic();plan=load_frozen(plan_path,plan_sha);admission=load_frozen(admission_path,admission_sha)
    root=Path(plan['output_root']);ownership=[];source_guard=Guard(plan,started)
    def identity():identity_gate(plan,plan_path,admission_path,plan_sha,admission_sha)
    def owned_cleanup():cleanup(plan,ownership)
    if root.exists() or root.is_symlink() or any(p.is_symlink() for p in root.parents):raise RuntimeError('SOURCE_FRESH_PRIVATE_NAMESPACE_REQUIRED_NO_RETRY')
    with deferred_termination_signals():
        # Both locks precede any source operation and remain held through audit,
        # persistence and unconditional owned-group cleanup. Busy means no launch.
        if not Path(plan['resource_contract']['shared_heavy_mutex']).is_file():raise RuntimeError('SOURCE_EXISTING_SHARED_HEAVY_LOCK_REQUIRED')
        with exclusive_heavy_lock(before_release=owned_cleanup) as heavy_fd:
            with shared_raw_family_lock(plan['raw_transfer_family_mutex'],owned_cleanup) as raw_fd:
                identity();runtime_profile_gate(plan,admission);resource_snapshot(plan,started,launch=True)
                root.mkdir()
                for name in ['archive','extracted','qc','tmp']:(root/name).mkdir()
                pending=root/'PENDING.json';seal=root/'terminal_seal.json';primary=root/'source_execution_receipt.json'
                binding=dict(plan_sha256=plan_sha,admission_sha256=admission_sha,scope='SOURCE_AUTHENTICATION_REFERENCE_QC_OPERATOR_ONLY')
                terminal=TerminalCommit(pending,seal,binding);frozen={str(pending):terminal.pending_sha};targets=[root,root/'archive',root/'extracted',root/'qc',root/'tmp',pending]
                try:
                    command=[plan['source_python'],'-B',plan['source_worker'],'--plan',str(plan_path),'--plan-sha256',plan_sha,'--started-monotonic',str(started),
                      '--admission',str(admission_path),'--admission-sha256',admission_sha,'--heavy-lock-fd',str(heavy_fd),'--raw-lock-fd',str(raw_fd)]
                    command=tuple(command)
                    monitor,msha=monitor_worker(plan,list(command),(heavy_fd,raw_fd),started,ownership,dict(identity=identity,plan_sha=plan_sha))
                    worker_map,more,result=consume_worker_result(plan,monitor,msha,plan_sha,command,resource_guard=source_guard);frozen.update(worker_map);targets+=more
                    freeze_registered(root,targets,frozen)
                    output=dict(status='SOURCE_ONLY_COMPLETE_NO_CALIBRATION_CREDIT',binding=binding,worker_result=result,monitor_sha256=msha,
                      precommit_ordinary_and_transport_output_sha256=dict(frozen),private_PENDING_removed_only_by_Terminal2=True,owned_cleanup_verified=True,source_G1_not_run=True,
                      full_source_build_qualified_unresolved=True,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False)
                    psha=save(primary,output);frozen[str(primary)]=psha;targets+=[primary,seal]
                    expected_seal=terminal_intended_sha(terminal,{str(primary):psha})
                    def final_identity():
                        identity();runtime_profile_gate(plan,admission);owned_cleanup()
                        if seal.exists():
                            if regular_sha(seal)!=expected_seal:raise RuntimeError('SOURCE_FIXED_INTENDED_TERMINAL_SEAL_CHANGED')
                            frozen[str(seal)]=expected_seal
                        freeze_registered(root,targets,frozen)
                        consume_worker_result(plan,monitor,msha,plan_sha,command,frozen[str(root/'qc/source_worker_result.json')],source_guard)
                    if not terminal.commit({str(primary):psha},identity_gate=final_identity,
                        resource_gate=lambda:resource_snapshot(plan,started),termination_gate=assert_no_termination):raise RuntimeError('SOURCE_TERMINAL_COMMIT_FAILED_PENDING_REMAINS')
                    # Deliberately no fallible audit/write after Terminal2 removes
                    # the private immutable-owned marker (reviewed contract).
                except BaseException as e:
                    owned_cleanup()
                    try:save(str(primary)+'.failure.json',dict(status='FAILED_PRESERVED_NO_SOURCE_CREDIT',error=repr(e),binding=binding,PENDING_authoritative=True))
                    except BaseException as second:safe_diagnostic(lambda:'SOURCE_FAILURE_ADDENDUM_WRITE_FAILED: '+repr(second))
                    raise


if __name__=='__main__':
    p=argparse.ArgumentParser(allow_abbrev=False);p.add_argument('--plan',default=str(PLAN));p.add_argument('--plan-sha256',required=True);p.add_argument('--admission',default=str(ADMISSION));p.add_argument('--admission-sha256',required=True);a=p.parse_args()
    execute(a.plan,a.plan_sha256,a.admission,a.admission_sha256)
