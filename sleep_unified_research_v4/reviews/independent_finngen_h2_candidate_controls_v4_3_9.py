"""Pure candidate admission/result helpers, private metadata only.

No executor/preparer, worker, native estimator, raw/reference reader, network,
actual lock or real mount/resource query is called. Old fixture numbers are
invented arithmetic metadata, never biological results or calibration evidence.
"""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time
from types import SimpleNamespace

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
S = P / 'scripts'
sys.path.insert(0, str(S))

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def run():
    started = time.monotonic()
    F = load('_qualified_finn_tiny_builder', R / 'independent_finngen_h2_candidate_controls_v4_3_3.py')
    F.OUT = R / 'independent_finngen_h2_candidate_controls_v4_3_9'
    F.OUT.mkdir(exist_ok=False)
    candidate = S / '58_run_finngen_feasibility_stage_v3_9.py'
    reviewed = {str(path): sha(path) for path in [candidate, S / '59_prepare_finngen_feasibility_stage_v3_9.py',
        S / '60_verify_finngen_h2_diagnostic_v5.py', S / '58_run_finngen_feasibility_stage_v3_8.py']}
    labels = ['healthy', 'prior_pending', 'prior_failure', 'prior_failed_status', 'prior_metadata_drift',
        'dependency_drift', 'review_drift', 'executor_admission_drift', 'plan_drift', 'two_jobs',
        'scope_rg', 'guard_drift', 'two_estimates', 'two_intersections', 'input_before_drift',
        'input_after_drift', 'args_rg', 'args_population_prev', 'args_sample_prev',
        'capture_mutation_after_parse', 'prior_derivative_drift_after_capture_parse',
        'prior_metadata_drift_after_capture_parse']
    cases = []
    for label in labels:
        C = load('_finn_h2_candidate_3_9_' + label, candidate)
        previous_lock = C.common.SHARED_LOCK
        f = F.fixture(label)
        d, plan = f['dir'], f['plan']
        campaign = d / 'private_campaign'
        approved = campaign / 'new_source_feasibility/finngen_R13_F5_INSOMNIA'
        namespace = approved / 'observed_h2_v3_9'
        namespace.mkdir(parents=True)
        C.common.SHARED_LOCK = campaign / 'native_heavy_worker.lock'
        dep = d / 'dependency.txt'
        dep.write_text('PRIVATE BOUND METADATA\n')
        baseline = d / 'baseline.json'
        write(baseline, {'private_metadata_fixture': True})
        reviews = []
        for i in range(7):
            path = d / ('review_' + str(i) + '.txt')
            path.write_text('PRIVATE EXACT REVIEW BINDING\n')
            reviews.append(str(path))
        admission = d / 'admission.json'
        plan.update(namespace=str(namespace), source_namespace=str(approved), shared_heavy_worker_lock=str(C.common.SHARED_LOCK),
            executor_path=str(candidate), admission=str(admission), required_independent_review_paths=reviews,
            guard=dict(worker_count=1, BLAS_threads=1, internal_floor_bytes=3 << 30, SSD_floor_bytes=5 << 30,
                observed_aggregate_worker_RSS_limit_bytes=2 << 30, new_output_limit_bytes=2 << 30, deadline_seconds=7200),
            dependencies_sha256={str(dep): sha(dep)}, baseline_gate_plan=str(baseline), baseline_gate_plan_sha256=sha(baseline))
        job = plan['jobs'][0]
        job.update(kind='h2', inputs=[str(f['derivative'])], estimates=1,
            output_prefix=str(d / 'fixture_h2'), out_prefix=str(d / 'fixture_h2'), worker_receipt=str(d / 'worker.json'))
        job['ldsc_args'] = ['--h2', str(f['derivative']), '--ref-ld-chr', 'PRIVATE_REF_', '--w-ld-chr', 'PRIVATE_REF_',
            '--n-blocks', '200', '--print-delete-vals', '--out', job['out_prefix']]
        if label == 'args_rg':
            job['ldsc_args'] += ['--rg', 'UNADMITTED']
        if label == 'args_population_prev':
            job['ldsc_args'] += ['--pop-prev', '0.1']
        if label == 'args_sample_prev':
            job['ldsc_args'] += ['--samp-prev', '0.1']
        if label == 'two_jobs':
            plan['jobs'].append(dict(job, job_id='UNADMITTED_EXTRA'))
        if label == 'scope_rg':
            plan['new_rg_commands'] = 1
        if label == 'guard_drift':
            plan['guard']['internal_floor_bytes'] -= 1
        if label == 'prior_failed_status':
            proof = plan['prior_preprocessing_terminal']
            stage = Path(proof['stage_receipt'])
            old = json.loads(stage.read_text())
            old['status'] = 'FAILED_PRESERVED'
            write(stage, old)
            seal = f['prior_seal']
            seal_data = json.loads(seal.read_text())
            seal_data['result_receipt_sha256'] = {str(stage): sha(stage)}
            write(seal, seal_data)
            proof['stage_receipt_sha256'] = sha(stage)
            proof['terminal_seal_sha256'] = sha(seal)
            proof['metadata_sha256'][str(stage)] = sha(stage)
            proof['metadata_sha256'][str(seal)] = sha(seal)
        write(f['plan_path'], plan)
        ph = sha(f['plan_path'])
        admit = dict(execution_admitted=True, plan_sha256=ph, scope=plan['scope'], independent_binding_review_pass=True,
            resource_plan_review_pass=True, executor_sha256=sha(candidate), independent_review_sha256={path: sha(path) for path in reviews})
        if label == 'executor_admission_drift':
            admit['executor_sha256'] = '0' * 64
        write(admission, admit)
        capture = json.loads(f['capture'].read_text())
        capture['plan_sha256'] = ph
        if label == 'two_estimates':
            capture['estimates'] *= 2
        if label == 'two_intersections':
            capture['final_intersections'] *= 2
        if label == 'input_before_drift':
            capture['input_sha256_before'] = {}
        if label == 'input_after_drift':
            capture['input_sha256_after'] = {}
        write(f['capture'], capture)
        if label == 'prior_pending':
            f['prior_pending'].write_text('PRIVATE PENDING VETO\n')
        if label == 'prior_failure':
            Path(str(f['prior_seal']) + '.failure.json').write_text('PRIVATE FAILURE VETO\n')
        if label == 'prior_metadata_drift':
            f['prior_plan'].write_text(f['prior_plan'].read_text() + ' ')
        if label == 'dependency_drift':
            dep.write_text('PRIVATE DRIFT\n')
        if label == 'review_drift':
            Path(reviews[0]).write_text('PRIVATE DRIFT\n')
        if label == 'plan_drift':
            f['plan_path'].write_text(f['plan_path'].read_text() + ' ')
        if label in ['capture_mutation_after_parse', 'prior_derivative_drift_after_capture_parse', 'prior_metadata_drift_after_capture_parse']:
            done = [False]
            def parsed(payload):
                value = json.loads(payload)
                if isinstance(value, dict) and 'estimates' in value and not done[0]:
                    done[0] = True
                    if label == 'capture_mutation_after_parse':
                        f['capture'].write_text(f['capture'].read_text() + ' ')
                    if label == 'prior_derivative_drift_after_capture_parse':
                        f['derivative'].write_text('PRIVATE INPUT DRIFT\n')
                    if label == 'prior_metadata_drift_after_capture_parse':
                        f['prior_plan'].write_text(f['prior_plan'].read_text() + ' ')
                return value
            C.json = SimpleNamespace(loads=parsed, dumps=json.dumps)
        accepted, error, result = False, None, None
        try:
            C.require_admission(f['plan_path'], ph)
            result = C.result_gate(plan, ph)
            accepted = True
        except BaseException as exc:
            error = type(exc).__name__ + ': ' + str(exc)
        finally:
            C.common.SHARED_LOCK = previous_lock
        cases.append(dict(label=label, accepted=accepted, expected_acceptance=label == 'healthy',
            passed=accepted == (label == 'healthy'), error=error, result_identity=result,
            fixture_plan_sha256=ph, only_pure_helpers=True))
    receipt = dict(schema='independent_finngen_h2_prospective_candidate_helpers_v4_3_9',
        controls=cases, all22_expected_results=all(case['passed'] for case in cases),
        reviewed_identity_sha256=reviewed, reviewed_inputs_unchanged=all(sha(path) == digest for path, digest in reviewed.items()),
        fixture_regular_sha256={str(path): sha(path) for path in sorted(F.OUT.rglob('*')) if path.is_file() and not path.is_symlink()},
        actual_scientific_workers=0, native_fits=0, arithmetic_verifier_runs=0, actual_mutex_calls=0,
        raw_reference_body_reads=0, network_calls=0, actual_plan_prepared=False,
        qualifications='Private miniature producer identities and invented capture values; only admission/result helper behavior. Stock/intersection/delete arithmetic is inherited from the unchanged60v5 consumer review, not recomputed here.',
        elapsed_seconds=time.monotonic() - started, max_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    with (R / 'independent_finngen_h2_candidate_controls_receipt_v4_3_9.json').open('x') as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({'controls': len(cases), 'all_pass': receipt['all22_expected_results']}))

if __name__ == '__main__':
    run()
