#!/usr/bin/env python3
"""Private metadata-only prerequisite controls; never reads real GWAS assets.

The fabricated capture numbers are arithmetic fixtures, not scientific results.
Only pure gates/TerminalCommit and the scalar verifier are called: no executor,
preparer, workers, subprocesses, fits, network, or actual heavy-lock calls.
"""
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

P = Path(__file__).resolve().parents[1]
S = P / 'scripts'
OUT = P / 'reviews/independent_finngen_h2_candidate_controls_v4_3_3'
sys.path.insert(0, str(S))
from terminal_commit_common_v2 import TerminalCommit


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def write(p, x):
    Path(p).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')


def module(name, p):
    spec = importlib.util.spec_from_file_location(name, p)
    x = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(x)
    return x


def se(xs):
    m = math.fsum(xs) / len(xs)
    return math.sqrt(199 / 200 * math.fsum((x-m)**2 for x in xs))


def commit(pending, seal, binding, receipt):
    t = TerminalCommit(pending, seal, binding)
    assert t.commit({str(receipt):sha(receipt)}, identity_gate=lambda:None,
                    resource_gate=lambda:None, termination_gate=lambda:None)


def fixture(label, changes=None):
    """Genuine small Terminal2 seals around entirely invented metadata files."""
    d = OUT / label
    d.mkdir(parents=True, exist_ok=False)
    derivative = d/'fixture_derivative.txt'
    derivative.write_text('PRIVATE ARITHMETIC FIXTURE; NO GWAS CONTENT\n')
    prior_executor = d/'prior_executor.txt'
    prior_executor.write_text('PRIVATE FIXTURE EXECUTOR IDENTITY\n')
    prior_admit = d/'prior_admission.json'
    write(prior_admit, {'fixture_only':True})
    pp = d/'prior_plan.json'
    ps = d/'prior_stage.json'
    pr = d/'prior_result.json'
    pending = d/'prior_PENDING.json'
    seal = d/'prior_seal.json'
    write(pp, dict(terminal_pending=str(pending), terminal_seal=str(seal),
                   executor_path=str(prior_executor), admission=str(prior_admit)))
    ph = sha(pp)
    write(pr, dict(status='QUALIFIED_PREPROCESSING_DIAGNOSTIC_PASS',plan_sha256=ph,
                   derivative=str(derivative),derivative_sha256=sha(derivative),
                   full_source_gzip_CRC_and_EOF_verified=True,
                   derivative_full_gzip_CRC_and_EOF_verified=True,
                   scientific_source_admitted=False,independent_replication_established=False,
                   verified_per_variant_N=False,verified_per_variant_INFO=False))
    result_map = {str(pr):sha(pr),str(derivative):sha(derivative)}
    write(ps, dict(status='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED',plan_sha256=ph,
                   admission_sha256=sha(prior_admit),result_sha256=result_map))
    commit(pending, seal, dict(plan_sha256=ph,admission_sha256=sha(prior_admit),
                              executor_sha256=sha(prior_executor)),ps)
    metadata = [pp,ps,pr,seal,prior_admit,prior_executor]
    capture = d/'fixture_h2.full_precision.json'
    stage = d/'h2_stage.json'
    plan_path = d/'h2_plan.json'
    h_pending = d/'h2_PENDING.json'
    h_seal = d/'h2_seal.json'
    log = d/'fixture_h2.log'
    deletes = d/'fixture_h2.delete'
    log.write_text('PRIVATE ARITHMETIC FIXTURE STOCK LOG\n')
    totals = [.12 + (.001 if i%2 else -.001) for i in range(200)]
    intercepts = [1.05 + (.0002 if i%2 else -.0002) for i in range(200)]
    deletes.write_text(''.join(format(x,'.17g')+'\n' for x in totals))
    plan = dict(stage='observed_h2',scope='FINNGEN_INSOMNIA_SOURCE_FEASIBILITY_ONLY',
        new_rg_commands=0,assumed_effective_N=185146.70377332723,derivative=str(derivative),
        input_sha256={str(derivative):sha(derivative)},stage_receipt=str(stage),
        executor_path=str(S/'58_run_finngen_feasibility_stage_v3_3.py'),
        terminal_pending=str(h_pending),terminal_seal=str(h_seal),
        jobs=[dict(job_id='finngen_insomnia_observed_h2',result_receipt=str(capture),
                   ldsc_args=['--h2',str(derivative),'--n-blocks','200','--print-delete-vals'])],
        prior_preprocessing_terminal=dict(plan_path=str(pp),plan_sha256=ph,
            stage_receipt=str(ps),stage_receipt_sha256=sha(ps),preprocessing_receipt=str(pr),
            executor_sha256=sha(prior_executor),terminal_seal_sha256=sha(seal),
            result_sha256=result_map,metadata_sha256={str(x):sha(x) for x in metadata}))
    if changes == 'two_jobs':
        plan['jobs'].append(dict(plan['jobs'][0],job_id='UNADMITTED_EXTRA_FIXTURE_JOB'))
    write(plan_path,plan)
    hh=sha(plan_path)
    row=dict(input=str(derivative),tot=.12,tot_se=se(totals),intercept=1.05,
             intercept_se=se(intercepts),mean_chisq=1.2,lambda_gc=1.1,n_blocks=200,
             constrain_intercept=False,n_snp=None,tot_delete_values=[[x] for x in totals],
             intercept_delete_values=[[x] for x in intercepts])
    r=dict(plan_sha256=hh,job_id=plan['jobs'][0]['job_id'],estimates=[row],
           final_intersections=[dict(input=str(derivative),final_ordered_SNP_count=12345)],
           input_sha256_before=plan['input_sha256'],input_sha256_after=plan['input_sha256'],
           warnings_and_errors=[],stock_log_sha256=sha(log),stock_delete_array_sha256={str(deletes):sha(deletes)},
           arguments=dict(h2=str(derivative),rg=None,n_blocks=200,samp_prev=None,pop_prev=None))
    if changes == 'wrong_stock_hash':
        r['stock_log_sha256']='0'*64
        r['stock_delete_array_sha256'][str(deletes)]='0'*64
    if changes == 'stock_vector_disagreement':
        deletes.write_text(''.join('0.9\n' for _ in range(200)))
        r['stock_delete_array_sha256'][str(deletes)]=sha(deletes)
    write(capture,r)
    consumed = {str(capture):sha(capture),str(log):sha(log),str(deletes):sha(deletes)}
    if changes == 'omit_capture':
        del consumed[str(capture)]
    if changes == 'omit_stock':
        del consumed[str(log)]
        del consumed[str(deletes)]
    write(stage,dict(status='QUALIFIED_FEASIBILITY_STAGE_COMPLETE_VERIFIED',plan_sha256=hh,
                     admission_sha256='a'*64,owned_cleanup_verified=True,result_sha256=consumed))
    commit(h_pending,h_seal,dict(plan_sha256=hh,admission_sha256='a'*64,
                                 executor_sha256=sha(plan['executor_path'])),stage)
    return dict(dir=d,plan=plan,plan_path=plan_path,hash=hh,stage=stage,capture=capture,
                prior_pending=pending,prior_seal=seal,prior_plan=pp,derivative=derivative,
                pending=h_pending,seal=h_seal,log=log,deletes=deletes)


def run():
    started=time.monotonic()
    before={str(S/n):sha(S/n) for n in [
        '58_run_finngen_feasibility_stage_v3_3.py','59_prepare_finngen_feasibility_stage_v3_3.py',
        '60_verify_finngen_h2_diagnostic_v2.py','terminal_commit_common_v2.py',
        'canonical_calibration_common_v4_5.py','sensitivity_capture_common.py',
        '38_sensitivity_ldsc_capture.py']}
    OUT.mkdir(exist_ok=False)
    cases=[]
    labels=['healthy','prior_pending','prior_failure','prior_metadata_drift','prior_derivative_symlink',
        'h2_pending','h2_failure','capture_drift','output_initial_symlink','omit_capture','omit_stock',
        'wrong_stock_hash','stock_vector_disagreement','two_jobs','output_post_arithmetic_symlink']
    inherited=['omit_capture','omit_stock','wrong_stock_hash','stock_vector_disagreement','two_jobs']
    for label in labels:
        f=fixture(label,label if label in inherited else None)
        if label=='prior_pending':f['prior_pending'].write_text('VETO\n')
        if label=='prior_failure':Path(str(f['prior_seal'])+'.failure.json').write_text('VETO\n')
        if label=='prior_metadata_drift':f['prior_plan'].write_text(f['prior_plan'].read_text()+' ')
        if label=='prior_derivative_symlink':
            target=f['dir']/'same_bytes_derivative.txt';target.write_bytes(f['derivative'].read_bytes());f['derivative'].unlink();f['derivative'].symlink_to(target)
        if label=='h2_pending':f['pending'].write_text('VETO\n')
        if label=='h2_failure':Path(str(f['seal'])+'.failure.json').write_text('VETO\n')
        if label=='capture_drift':f['capture'].write_text(f['capture'].read_text()+' ')
        if label=='output_initial_symlink':
            target=f['dir']/'same_bytes_log.txt';target.write_bytes(f['log'].read_bytes());f['log'].unlink();f['log'].symlink_to(target)
        verifier=module('_finn_scalar_'+label,S/'60_verify_finngen_h2_diagnostic_v2.py')
        if label=='output_post_arithmetic_symlink':
            original=verifier.centered_se;calls=[0]
            def wrap(xs):
                val=original(xs);calls[0]+=1
                if calls[0]==2:
                    target=f['dir']/'same_bytes_log.txt';target.write_bytes(f['log'].read_bytes());f['log'].unlink();f['log'].symlink_to(target)
                return val
            verifier.centered_se=wrap
        accepted=False;error=None
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                verifier.verify(f['plan_path'],f['hash'],f['dir']/'diagnostic.json')
            accepted=True
        except BaseException as e:error=type(e).__name__+': '+str(e)
        expected=label=='healthy'
        cases.append(dict(label=label,accepted=accepted,expected_acceptance=expected,
                          expected_control_pass=accepted==expected,error=error,
                          fixture_plan_sha256=f['hash'],fixture_only=True,
                          diagnostic_written=(f['dir']/'diagnostic.json').exists()))
    artifacts={str(p):sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and not p.is_symlink()}
    after={p:sha(p) for p in before}
    receipt=dict(schema='independent_finngen_candidate_metadata_controls_v4_3_3',
        reviewed_input_sha256=before,reviewed_inputs_unchanged=before==after,
        controls=cases,controls_count=len(cases),unexpected_acceptances=sum(c['accepted'] and not c['expected_acceptance'] for c in cases),
        private_fixture_regular_file_sha256=artifacts,
        scope=dict(real_GWAS_or_reference_bytes_read=False,workers_launched=0,fits_launched=0,
                   real_heavy_mutex_acquired=False,private_Terminal2_calls_only=True),
        elapsed_seconds=time.monotonic()-started,max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    write(P/'reviews/independent_finngen_h2_candidate_controls_receipt_v4_3_3.json',receipt)
    print(json.dumps(dict(controls=len(cases),unexpected_acceptances=receipt['unexpected_acceptances'],inputs_unchanged=before==after)))


if __name__=='__main__':run()
