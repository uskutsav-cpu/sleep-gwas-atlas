#!/usr/bin/env python3
"""Corrected-schema private metadata controls; no real worker/fit/lock."""
import contextlib
import hashlib
import importlib.util
import io
import json
import math
from pathlib import Path
import resource
import sys
import time

R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_finn_old_fixture',R/'independent_finngen_h2_candidate_controls_v4_3_3.py')
F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)
F.OUT=R/'independent_finngen_h2_corrected_candidate_controls_v4_3_4'
sys.path.insert(0,str(F.S))
from terminal_commit_common_v2 import require_committed


def fixture(label):
    f=F.fixture(label);d=f['dir'];plan=dict(f['plan'])
    prefix=str(d/'corrected_h2');capture=Path(prefix+'.full_precision.json')
    stage=d/'corrected_stage.json';seal=d/'corrected_seal.json';pending=d/'corrected_PENDING.json'
    log=Path(prefix+'.log');delete=Path(prefix+'.delete');part=Path(prefix+'.part_delete')
    worker=d/'corrected_worker.json';journal=d/'corrected_worker.failure_journal.jsonl';stdout=d/'corrected_stdout.log'
    log.write_text('PRIVATE FIXTURE ONLY; NO ESTIMATOR EXECUTED\n')
    journal.write_text('{"fixture_only":true}\n');stdout.write_text('PRIVATE FIXTURE ONLY\n')
    r=json.loads(f['capture'].read_text());x=r['estimates'][0]
    x.update(n_annot=1,part_delete_values=x['tot_delete_values'])
    for p,field in [(delete,'tot_delete_values'),(part,'part_delete_values')]:
        p.write_text(''.join(format(a[0],'.17g')+'\n' for a in x[field]))
    argv=['--h2',plan['derivative'],'--ref-ld-chr','PRIVATE_FIXTURE_REFERENCE',
          '--w-ld-chr','PRIVATE_FIXTURE_REFERENCE','--n-blocks','200','--print-delete-vals','--out',prefix]
    job=dict(job_id='finngen_insomnia_observed_h2',kind='h2',inputs=[plan['derivative']],estimates=1,
        out_prefix=prefix,output_prefix=prefix,result_receipt=str(capture),worker_receipt=str(worker),ldsc_args=argv,
        command_template=['PRIVATE_FIXTURE_PYTHON','-u','PRIVATE_FIXTURE_CAPTURE','--plan',str(d/'corrected_plan.json'),
                          '--plan-sha','{PLAN_SHA256}','--job-id','finngen_insomnia_observed_h2']+argv)
    plan.update(jobs=[job],stage_receipt=str(stage),terminal_pending=str(pending),terminal_seal=str(seal),
                executor_path=str(F.S/'58_run_finngen_feasibility_stage_v3_4.py'))
    if label=='two_jobs':plan['jobs'].append(dict(job,job_id='UNADMITTED_EXTRA_JOB'))
    if label=='liability_args':job['ldsc_args']+=['--pop-prev','0.1']
    path=d/'corrected_plan.json';F.write(path,plan);ph=F.sha(path)
    r.update(plan_sha256=ph,job_id=job['job_id'],arguments=dict(h2=plan['derivative'],rg=None,n_blocks=200,out=prefix,
             samp_prev=None,pop_prev=None,print_delete_vals=True),stock_log_sha256=F.sha(log),
             stock_delete_array_sha256={str(delete):F.sha(delete),str(part):F.sha(part)},
             final_intersections=[dict(input=plan['derivative'],final_ordered_SNP_count=12345,n_blocks=200,
                 final_ordered_SNP_sha256='1'*64,final_ordered_Z_float64_big_endian_sha256='2'*64,
                 final_ordered_N_float64_big_endian_sha256='3'*64)])
    if label=='wrong_stock_hash':r['stock_log_sha256']='0'*64
    if label=='stock_vector_disagreement':
        delete.write_text('0.9\n'*200);r['stock_delete_array_sha256'][str(delete)]=F.sha(delete)
    if label=='bad_intersection':r['final_intersections'][0]['final_ordered_SNP_sha256']='invalid'
    if label=='capture_liability':r['arguments']['pop_prev']=0.1
    F.write(capture,r)
    outputs={str(p):F.sha(p) for p in [capture,log,delete,part,journal,stdout]}
    if label=='worker_map_disagreement':outputs[str(capture)]='0'*64
    F.write(worker,dict(status='WORKER_COMPLETE_VERIFIED',plan_sha256=ph,returncode=0,stop_reason=None,
         owned_cleanup_verified=True,process_group_teardown=dict(remaining_group_members=[]),
         command=[ph if a=='{PLAN_SHA256}' else a for a in job['command_template']],output_sha256=outputs,
         failure_journal=str(journal)))
    frozen={str(p):F.sha(p) for p in [capture,log,delete,part,worker,journal,stdout]}
    if label=='omit_capture':del frozen[str(capture)]
    if label=='omit_stock':del frozen[str(log)];del frozen[str(delete)]
    if label=='missing_partition':del frozen[str(part)]
    F.write(stage,dict(status='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED',plan_sha256=ph,
         admission_sha256='a'*64,owned_cleanup_verified=True,result_sha256=frozen))
    binding=dict(plan_sha256=ph,admission_sha256='a'*64,executor_sha256=F.sha(plan['executor_path']))
    F.commit(pending,seal,binding,stage)
    f.update(plan=plan,plan_path=path,hash=ph,stage=stage,capture=capture,pending=pending,seal=seal,
             log=log,deletes=delete,part=part,worker=worker,binding=binding)
    return f


def replace_symlink(path):
    target=path.with_name(path.name+'.same_bytes');target.write_bytes(path.read_bytes())
    path.unlink();path.symlink_to(target)


def run():
    started=time.monotonic();F.OUT.mkdir(exist_ok=False)
    inputs={str(F.S/n):F.sha(F.S/n) for n in ['58_run_finngen_feasibility_stage_v3_4.py',
        '59_prepare_finngen_feasibility_stage_v3_4.py','60_verify_finngen_h2_diagnostic_v3.py']}
    cases=[]
    labels=['healthy','omit_capture','omit_stock','wrong_stock_hash','stock_vector_disagreement','two_jobs',
        'output_post_arithmetic_symlink','prior_pending','h2_pending','missing_partition','bad_intersection',
        'capture_liability','liability_args','worker_map_disagreement','postpersist_input_symlink',
        'postpersist_output_mutation','postpersist_second_write_failure']
    for label in labels:
        f=fixture(label);v=F.module('_corrected_finn_'+label,F.S/'60_verify_finngen_h2_diagnostic_v3.py')
        out=f['dir']/'diagnostic.json';table=out.with_suffix('.tsv')
        if label=='prior_pending':f['prior_pending'].write_text('VETO\n')
        if label=='h2_pending':f['pending'].write_text('VETO\n')
        if label=='output_post_arithmetic_symlink':
            old=v.centered_se;calls=[0]
            def wrap(xs):
                value=old(xs);calls[0]+=1
                if calls[0]==2:replace_symlink(f['log'])
                return value
            v.centered_se=wrap
        old_regular_hashes=v.regular_hashes
        if label in ['postpersist_input_symlink','postpersist_output_mutation']:
            def gate(mapping):
                if str(out) in mapping:
                    if label=='postpersist_input_symlink':replace_symlink(f['log'])
                    else:out.write_text(out.read_text()+' ')
                return old_regular_hashes(mapping)
            v.regular_hashes=gate
        if label=='postpersist_second_write_failure':table.write_text('PRESERVED EXISTING REPORT COLLISION\n')
        accepted=False;error=None
        try:
            with contextlib.redirect_stdout(io.StringIO()):v.verify(f['plan_path'],f['hash'],out)
            accepted=True
        except BaseException as e:error=type(e).__name__+': '+str(e)
        written=out.exists();arithmetic_true=False
        if written:
            try:arithmetic_true=json.loads(out.read_text())['arithmetic_pass'] is True
            except (ValueError,KeyError):pass
        if label=='healthy':
            report=json.loads(out.read_text());assert report['diagnostic']['native_SNP_count']==12345
            assert report['intercept_delete_export_scope'].startswith('Intercept arithmetic uses captured native array')
        cases.append(dict(label=label,accepted=accepted,expected_acceptance=label=='healthy',
            gate_control_pass=accepted==(label=='healthy'),error=error,diagnostic_written=written,
            persisted_arithmetic_pass_true=arithmetic_true,
            diagnostic_has_durable_terminal_or_failure_veto=False,private_fixture_only=True,
            plan_sha256=f['hash']))
    after={p:F.sha(p) for p in inputs}
    F.write(R/'independent_finngen_h2_corrected_candidate_controls_receipt_v4_3_4.json',
        dict(schema='independent_finngen_corrected_candidate_metadata_controls_v4_3_4',controls=cases,
            all_gate_controls_pass=all(c['gate_control_pass'] for c in cases),controls_count=len(cases),
            provisional_success_report_survives_failed_verifier=[c['label'] for c in cases if not c['accepted'] and c['persisted_arithmetic_pass_true']],
            reviewed_input_sha256=inputs,reviewed_inputs_unchanged=inputs==after,
            private_regular_file_sha256={str(p):F.sha(p) for p in sorted(F.OUT.rglob('*')) if p.is_file() and not p.is_symlink()},
            real_body_reads=False,workers=0,fits=0,network=0,real_heavy_mutex=False,
            elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    print(json.dumps(dict(controls=len(cases),all_gate_controls_pass=all(c['gate_control_pass'] for c in cases),
        leftover_provisional_success=sum(not c['accepted'] and c['persisted_arithmetic_pass_true'] for c in cases))))


if __name__=='__main__':run()
