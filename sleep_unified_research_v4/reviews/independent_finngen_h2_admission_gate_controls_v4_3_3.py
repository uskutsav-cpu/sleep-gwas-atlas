#!/usr/bin/env python3
"""Pure candidate58 gates on private metadata fixtures; no execute/prepare."""
import importlib.util
import json
from pathlib import Path
import sys

R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('_fixture_builder',R/'independent_finngen_h2_candidate_controls_v4_3_3.py')
F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)
sys.path.insert(0,str(F.S))
C=F.module('_candidate_admission',F.S/'58_run_finngen_feasibility_stage_v3_3.py')
original_lock=C.common.SHARED_LOCK
labels=['healthy','dependency_drift','review_drift','prior_pending','prior_failure',
        'admission_executor_drift','frozen_plan_drift','scope_drift','guard_drift','two_jobs',
        'capture_input_drift']
cases=[]
try:
    for label in labels:
        f=F.fixture('admission_'+label)
        d=f['dir'];campaign=d/'private_campaign';campaign.mkdir()
        C.common.SHARED_LOCK=campaign/'native_heavy_worker.lock'
        approved=campaign/'new_source_feasibility/finngen_R13_F5_INSOMNIA'
        namespace=approved/'observed_h2_v3_3';namespace.mkdir(parents=True)
        dep=d/'bound_dependency.txt';dep.write_text('PRIVATE METADATA DEPENDENCY\n')
        baseline=d/'baseline.json';F.write(baseline,{'fixture_only':True})
        admit=d/'admission.json'
        reviews=[]
        for i in range(7):
            p=d/('review_'+str(i)+'.txt');p.write_text('PRIVATE REVIEW IDENTITY FIXTURE\n');reviews.append(str(p))
        plan=f['plan']
        plan.update(namespace=str(namespace),source_namespace=str(approved),shared_heavy_worker_lock=str(C.common.SHARED_LOCK),
            guard=dict(worker_count=1,BLAS_threads=1,internal_floor_bytes=3<<30,SSD_floor_bytes=5<<30,
                observed_aggregate_worker_RSS_limit_bytes=2<<30,new_output_limit_bytes=2<<30,deadline_seconds=7200),
            admission=str(admit),required_independent_review_paths=reviews,
            dependencies_sha256={str(dep):F.sha(dep)},baseline_gate_plan=str(baseline),baseline_gate_plan_sha256=F.sha(baseline))
        if label=='scope_drift':plan['new_rg_commands']=1
        if label=='guard_drift':plan['guard']['internal_floor_bytes']-=1
        if label=='two_jobs':plan['jobs'].append(dict(plan['jobs'][0],job_id='UNADMITTED_EXTRA'))
        F.write(f['plan_path'],plan);expected=F.sha(f['plan_path'])
        a=dict(execution_admitted=True,plan_sha256=expected,scope=plan['scope'],independent_binding_review_pass=True,
               resource_plan_review_pass=True,independent_review_sha256={p:F.sha(p) for p in reviews},
               executor_sha256=F.sha(F.S/'58_run_finngen_feasibility_stage_v3_3.py'))
        if label=='admission_executor_drift':a['executor_sha256']='0'*64
        F.write(admit,a)
        capture=json.loads(f['capture'].read_text());capture['plan_sha256']=expected
        if label=='capture_input_drift':capture['input_sha256_after']={}
        F.write(f['capture'],capture)
        if label=='dependency_drift':dep.write_text('DRIFT\n')
        if label=='review_drift':Path(reviews[0]).write_text('DRIFT\n')
        if label=='prior_pending':f['prior_pending'].write_text('VETO\n')
        if label=='prior_failure':Path(str(f['prior_seal'])+'.failure.json').write_text('VETO\n')
        if label=='frozen_plan_drift':f['plan_path'].write_text(f['plan_path'].read_text()+' ')
        accepted=False;error=None
        try:
            C.require_admission(f['plan_path'],expected)
            C.result_gate(plan,expected)
            accepted=True
        except BaseException as e:error=type(e).__name__+': '+str(e)
        cases.append(dict(label=label,accepted=accepted,expected_acceptance=label=='healthy',
                          passed=accepted==(label=='healthy'),error=error,
                          plan_sha256=expected,private_fixture_only=True))
finally:C.common.SHARED_LOCK=original_lock
receipt=dict(schema='independent_finngen_pure_admission_gates_v4_3_3',controls=cases,all_expected_pass=all(c['passed'] for c in cases),
    candidate_sha256=F.sha(F.S/'58_run_finngen_feasibility_stage_v3_3.py'),
    fixture_builder_sha256=F.sha(R/'independent_finngen_h2_candidate_controls_v4_3_3.py'),
    fixture_regular_sha256={str(p):F.sha(p) for p in sorted(F.OUT.glob('admission_*/*')) if p.is_file() and not p.is_symlink()},
    real_GWAS_reference_reads=False,worker_calls=0,fits=0,real_heavy_lock_opened=False,
    note='Only common.SHARED_LOCK path identity was redirected; no lock API or execution API was called. Final actual-plan review still pending.')
F.write(R/'independent_finngen_h2_admission_gate_controls_receipt_v4_3_3.json',receipt)
print(json.dumps(dict(controls=len(cases),all_expected_pass=receipt['all_expected_pass'])))
