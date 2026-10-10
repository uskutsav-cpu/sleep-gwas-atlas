"""Independent v1 review: aggregate metadata and entirely invented controls.

No source body, runtime census, network URL, production worker or mutex runs.
The candidate is imported read-only; only test-process objects are mocked.
"""
import ast
from contextlib import contextmanager
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
from unittest.mock import patch

P = Path(__file__).resolve().parents[1]
R = P / 'reviews'
S = P / 'scripts'
FIX = P.parents[1] / 'independent_signed_reference_source_fixtures_v1_2'
PLAN = P / 'manifests/signed_reference_source_execution_plan_v1.json'
PLAN_SHA = '6d73f45ee3f73f7db7c0cc727d9b2d3e1af30d3ee27391cf1154f1f8ab8482b4'
AUTHOR = R / 'signed_reference_source_author_preparation_seal_v1.json'
AUTHOR_SHA = '798be8797263ab9d8d4f4f05cf2c0966ed0d93c09dd40508fca0878f82c47311'
START = time.monotonic()
CASES = []
ASSERTIONS = 0
RESOURCE = json.loads((R / 'independent_signed_reference_source_prelaunch_resource_plan_v1.json').read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(65536), b''):
            h.update(b)
    return h.hexdigest()


def check(value, message):
    global ASSERTIONS
    ASSERTIONS += 1
    if not value:
        raise AssertionError(message)


def guard():
    check(time.monotonic()-START < RESOURCE['own_deadline_seconds'], 'own deadline')
    check(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss < RESOURCE['own_RSS_stop_bytes'], 'own RSS')
    check(sum(p.stat().st_size for p in FIX.rglob('*') if p.is_file()) < RESOURCE['own_output_cap_bytes'], 'own fixture cap')


def group(name, fn):
    before = ASSERTIONS
    detail = fn()
    guard()
    CASES.append(dict(case=name, assertions=ASSERTIONS-before, status='PASS', detail=detail))


def put(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')
    return sha(path)


def reject(fn):
    try:
        fn()
    except Exception as e:
        return type(e).__name__+': '+str(e)
    raise AssertionError('expected refusal')


def bindings():
    check(sha(PLAN) == PLAN_SHA, 'actual frozen plan')
    check(sha(AUTHOR) == AUTHOR_SHA, 'actual author seal')
    plan = json.loads(PLAN.read_text())
    author = json.loads(AUTHOR.read_text())
    check(author['actual_plan_sha256'] == PLAN_SHA, 'author actual plan')
    check(len(author['file_sha256']) == author['artifact_count'] == 53, '53 author artifacts')
    for path, digest in author['file_sha256'].items():
        check(sha(path) == digest, 'author sealed file '+path)
    deferred = plan['preparation_deferred_dependency_sha256']
    fresh = {}
    for path, digest in plan['dependency_sha256'].items():
        if path not in deferred:
            check(sha(path) == digest, 'nonbody plan dependency '+path)
            fresh[path] = digest
    check(len(fresh) == 28 and len(deferred) == 26, 'fresh/deferred cardinality')
    for path, digest in deferred.items():
        check(plan['dependency_sha256'][path] == digest, 'deferred source map binding')
    # Reuse the already sealed complete B identity review and root consumption.
    bseal = R / 'independent_decoderB_runtime_identity_review_v1_seal.json'
    check(sha(bseal) == '7d7f388d4380632682f3dde920b32f93e2b7e454d90b2a72dac1afde3c7a0a51', 'prior independent B seal')
    bs = json.loads(bseal.read_text())
    bm = bs.get('file_sha256', bs.get('review_artifact_sha256'))
    check(bm is not None and len(bm) == 6, 'B six review artifacts')
    for path, digest in bm.items():
        check(sha(path) == digest, 'prior B review artifact')
    rootproof = P / 'logs/signed_reference_source_author_seal_root_consumed_v1.json'
    check(sha(rootproof) == 'dd3f668873987c0512740eb2118dc34dbda845cb682da87d0442289a5bcd978e', 'root author consumption')
    # Read only aggregate metadata about the author's 4,770-file census.
    inventory_candidates = [p for p in author['file_sha256'] if 'inventory' in Path(p).name]
    check(author['invented_fixture_regular_file_count'] == 4770, 'author fixture census stated')
    check(sha(PLAN) == PLAN_SHA and sha(AUTHOR) == AUTHOR_SHA, 'post-read plan/author')
    return dict(fresh_nonbody_plan_dependencies=fresh, inherited_deferred_dependencies=deferred,
                author_artifact_count=53, author_control_groups=42,
                author_invented_fixture_census_not_repeated=4770,
                author_inventory_metadata_paths=inventory_candidates,
                completed_B_review_reused_no_runtime_census=True)


sys.path.insert(0, str(S))
import signed_reference_source_common_v1 as common
import signed_reference_source_controller_v1 as controller
import signed_reference_source_worker_v1 as worker
import signed_reference_decode_v1 as decode
import terminal_commit_common_v2 as terminal


def ambient_curl():
    captured = []
    class Stop(Exception):
        pass
    def capture(command, **kwargs):
        captured.append((command, kwargs))
        raise Stop('capture before subprocess')
    with patch.object(worker.subprocess, 'run', capture):
        try:
            worker.curl_file(dict(curl='/usr/bin/curl', candidate=dict(bytes=288277344)),
                             'https://example.invalid/body', FIX/'never_written')
        except Stop:
            pass
    command, kw = captured[0]
    check(command[1] == '--fail' and '-q' not in command and '--disable' not in command, 'v1 argv lacks first -q')
    check('env' not in kw, 'curl inherits ambient environment')
    ownhome = FIX / 'invented_curl_home'
    ownhome.mkdir()
    config = ownhome / '.curlrc'
    # The only executed curl operation is --version: no URL or transfer.
    config.write_text('insecure\nlocation\nretry = 7\nthis_is_a_deliberately_invalid_review_option\n')
    env = {**os.environ, 'CURL_HOME': str(ownhome)}
    unguarded = subprocess.run(['/usr/bin/curl', '--version'], env=env, capture_output=True, text=True)
    guarded = subprocess.run(['/usr/bin/curl', '-q', '--version'], env=env, capture_output=True, text=True)
    check('review_option' in unguarded.stderr, 'actual curl consumed own config')
    check(guarded.returncode == 0 and guarded.stdout.startswith('curl ') and not guarded.stderr, 'first -q disables own config')
    for name, result in [('config_consumed', unguarded), ('config_disabled', guarded)]:
        put(FIX/(name+'.json'), dict(argv_no_URL=True, returncode=result.returncode, stdout=result.stdout, stderr=result.stderr))
    return dict(measured_finding='V1_AMBIENT_CURL_CONFIG_NOT_DISABLED', captured_argv=command,
                config_sha256=sha(config), unguarded_returncode=unguarded.returncode,
                first_q_returncode=guarded.returncode, network_requests=0,
                qualifier='No real config was read or transfer run; own invalid config proves ordinary ambient config consumption.')


def representation():
    p = R/'signed_reference_source_storage_exposure_v1.json'
    value = json.loads(p.read_text())
    static, final = [Path(k) for k in value['fixture_file_sha256'] if Path(k).name == 'static_selection.jsonl'][0], [Path(k) for k in value['fixture_file_sha256'] if Path(k).name == 'J_ordered.jsonl'][0]
    for k, v in value['fixture_file_sha256'].items():
        check(sha(k) == v, 'invented storage source SHA')
    n = sum(1 for _ in static.open())
    check(n == sum(1 for _ in final.open()) == 17, '17 invented retained rows')
    total = static.stat().st_size + final.stat().st_size
    per = Fraction(total, n)
    cap = 512 * (1 << 20)
    maximum = cap * per.denominator // per.numerator
    hypothetical = math.ceil(1290028*per)
    check(total == 10966 and maximum == 832282 and hypothetical == 832143944, 'independent exact rational capacity')
    check(hypothetical > cap, 'hypothetical two streams exceed subcap')
    # No real count is inferred. A prospective compact schema must retain the
    # same source/master ranks, all QC values and exclusions; not remove rows.
    return dict(invented_rows=n, combined_bytes=total, exact_combined_bytes_per_row=str(per),
                strict_QC_cap=cap, max_similar_eligible_rows_before_any_other_QC=maximum,
                hypothetical_all_master_eligible_two_stream_bytes=hypothetical,
                actual_eligible_rows_unknown=True, actual_source_failure_observed=False,
                compact_equivalent_prospective_successor_recommended=True)


def independent_arithmetic():
    import numpy as np
    # Independently construct 503 doses and scalar sample moments without
    # using the candidate's decoded values to define the expected answer.
    dose = np.array(([0, 1, 2, 2, 0, 1]*84)[:503], dtype=np.int8)
    dose[:25] = -1
    obs = [int(v) for v in dose if v != -1]
    n = len(obs); total = sum(obs); sq = sum(v*v for v in obs)
    mean = total/n
    expected = np.array([(float(v)-mean) if v != -1 else 0. for v in dose])
    expected /= math.sqrt(sum(float(v)*float(v) for v in expected)/503)
    for sign in [-1, 1]:
        check(np.allclose(decode.standardized(dose, sign, np), sign*expected, rtol=1e-13, atol=1e-14), 'independent normalized signed dose')
    q = decode.qc_column(dose, np)
    check((q['observed_count'], q['missing_count'], q['sum_dose'], q['sum_dose2'], q['integer_variance_numerator']) == (n,25,total,sq,n*sq-total*total), 'integer sufficient moments exact')
    # A1=-centered A2, with the same mean-imputed 503 denominator.
    a2 = np.array([2-int(v) if v >= 0 else -1 for v in dose])
    check(np.allclose(decode.standardized(a2, 1, np), -expected, rtol=1e-13, atol=1e-14), 'A1/A2 sign conversion')
    return dict(rows=1, individuals=503, imputed=25, discrete_moments_exact=True,
                denominator=503, A1_to_A2_standardized_sign=-1,
                qualified_source_level_decoderB_controls_inherited_not_repeated=True)


def metadata(root, ps, command, contradiction=None):
    mapping = {}
    for p in sorted(controller.expected_worker_outputs(root)):
        p = Path(p);p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(b'INVENTED_ONLY\n');mapping[str(p)] = sha(p)
    candidate=dict(bytes=14, upstream_md5='md5:INVENTED', file_name='INVENTED')
    archive=root/'archive/1000G_Phase3_plinkfiles.tgz'
    ai=dict(bytes=14,md5='INVENTED',sha256=mapping[str(archive)])
    counts={str(c):211 for c in range(1,23)}
    def phase(name, value):
        p=root/'qc'/name;p.write_text(json.dumps(value,indent=2)+'\n');mapping[str(p)]=sha(p)
    phase('G2_authentication.json',dict(status='AUTHENTICATED_NEW_CANDIDATE_NOT_OLD_COMPACT_IDENTITY',candidate=candidate,archive_identity=ai))
    phase('G3_format.json',dict(status='EXACT_SOURCE_FORMAT_PASS',samples=503,all22_fam_byte_identical=True,padding_all_rows_checked=True,chromosome_variant_counts=counts,source_index_sha256=mapping[str(root/'extracted/source_index.sqlite')]))
    phase('G4_static_seal.json',dict(status='FROZEN_BEFORE_GENOTYPE_QC',master_reference_rows=1290028,static_ordered_sha256=mapping[str(root/'qc/static_selection.jsonl')],static_exclusions_sha256=mapping[str(root/'qc/static_exclusions.jsonl')]))
    phase('G4_J_seal.json',dict(status='REFERENCE_ONLY_J_FROZEN',J_ordered_sha256=mapping[str(root/'qc/J_ordered.jsonl')],chromosome_J_counts=counts,J_count=sum(counts.values()),no_MAF_threshold=True,no_GWAS_outcomes=True,no_native_mask_selection=True))
    rows=[]
    for c in range(1,23):
        pref=root/'qc/controls'/f'chr{c}'
        rows.append(dict(CHR=c,J_count=211,independent_rows=211,source_window_sha256=mapping[str(pref)+'.bed'],B_array_sha256=mapping[str(pref)+'.B.npz'],B_receipt_sha256=mapping[str(pref)+'.B.receipt.json'],controls=['first64','middle64','last64','early32_scatter','middle32_scatter','late32_scatter','combined_distant_and_windows'],full_chromosome_energy_symmetry=True,between_storage_window_correlations_retained=True))
    phase('G5_operator_controls.json',dict(status='ACTUAL_AUTHENTICATED_SOURCE_CONTROLS_PASS',fixed_selection_sha256=mapping[str(root/'qc/controls/fixed_G5_selections.json')],chromosomes=rows,maintained_B_padding_independent_rejection=False,zero_variance_fallback_not_credited=True,finite_reference_503_conditional_only=True,population_LD_or_actual_GWAS_mask_adequacy_established=False,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False))
    result=dict(status='SOURCE_G2_G5_PASS_PRIVATE_PROVISIONAL',plan_sha256=ps,ordinary_output_sha256=mapping,archive_identity=ai,source_index_sha256=mapping[str(root/'extracted/source_index.sqlite')],J_ordered_sha256=mapping[str(root/'qc/J_ordered.jsonl')],chromosome_J_counts=counts,S_simulations_run=0,C_calibration_credit=0,estimator_calls=0,biological41Cov_release_admitted=False)
    common.save(root/'qc/source_worker_result.json',result)
    monitor=dict(status='COMPLETE_REAPED_SOURCE_ONLY',command=command,returncode=0,stop_reason=None,reaped_immediate_child=True,signal_zero_group_absent=True,plan_sha256=ps,cleanup_errors=[],run_error=None,process_group_teardown=dict(remaining_group_members=[]),first_teardown_attempt={})
    if contradiction == 'wrong_command':monitor['command']=['INVENTED_WRONG_COMMAND']
    if contradiction == 'cleanup_error':monitor['cleanup_errors']=['INVENTED_BAD_CLEANUP']
    if contradiction == 'extra_directory':(root/'qc/UNREGISTERED_EMPTY_DIRECTORY').mkdir()
    return monitor


@contextmanager
def fake_lock(*args, **kwargs):
    yield 91
    cb=kwargs.get('before_release')
    if cb:cb()


def terminal_case(label, fault=None):
    root=FIX/label
    lock=FIX/(label+'.fake_lock');lock.write_text('INVENTED_NO_FLOCK\n')
    plan=dict(output_root=str(root),candidate=dict(bytes=14,upstream_md5='md5:INVENTED',file_name='INVENTED'),expected_master_reference_rows=1290028,
              source_python='INVENTED_PYTHON',source_worker='INVENTED_WORKER',raw_transfer_family_mutex=str(lock),resource_contract=dict(shared_heavy_mutex=str(lock)))
    pp=FIX/(label+'.plan.json');ap=FIX/(label+'.admission.json')
    ps=put(pp,plan);ads=put(ap,dict(invented=True));calls=[]
    def mock_monitor(p,command,fds,started,ownership,gates):
        (root/'worker_stdout.log').write_bytes(b'')
        m=metadata(root,ps,command,fault);mp=root/'monitor_receipt.json'
        return mp,common.save(mp,m)
    def resource_gate(*args,**kwargs):
        calls.append('resource')
        if fault == 'postseal_resource' and (root/'terminal_seal.json').exists():raise RuntimeError('INVENTED_POSTSEAL_RESOURCE')
        return dict(invented=True)
    def identity(*args,**kwargs):
        calls.append('identity')
        if fault == 'postseal_output_mutation' and (root/'terminal_seal.json').exists():
            (root/'qc/J_ordered.jsonl').write_bytes(b'INVENTED_LATE_MUTATION\n')
    with patch.object(controller,'exclusive_heavy_lock',fake_lock), patch.object(controller,'shared_raw_family_lock',fake_lock), patch.object(controller,'identity_gate',identity), patch.object(controller,'runtime_profile_gate',lambda *a:None), patch.object(controller,'resource_snapshot',resource_gate), patch.object(controller,'monitor_worker',mock_monitor):
        try:
            controller.execute(pp,ps,ap,ads);outcome='COMPLETE'
        except RuntimeError as e:
            outcome='REJECTED: '+str(e)
    committed=not (root/'PENDING.json').exists() and (root/'terminal_seal.json').is_file()
    if label == 'healthy':
        check(committed and outcome == 'COMPLETE','healthy terminal')
        primary=root/'source_execution_receipt.json'
        binding=dict(plan_sha256=ps,admission_sha256=ads,scope='SOURCE_AUTHENTICATION_REFERENCE_QC_OPERATOR_ONLY')
        terminal.require_committed(root/'PENDING.json',root/'terminal_seal.json',binding,{str(primary):sha(primary)})
    elif fault in ['wrong_command','extra_directory']:
        check(committed and outcome == 'COMPLETE','measured missing consumer check')
    else:
        check(not committed and (root/'PENDING.json').is_file() and outcome.startswith('REJECTED'), 'failure retains PENDING')
    return dict(outcome=outcome,terminal_committed=committed,PENDING_present=(root/'PENDING.json').exists(),
                ordinary_fixture_file_count=sum(1 for p in root.rglob('*') if p.is_file()),gate_calls=calls,
                actual_worker_network_body_mutex_operations=0,
                qualified_mock_only='Fake phase metadata is not actual G2-G5 proof; current producer fixes its command itself.')


def main():
    check(shutil.disk_usage('/System/Volumes/Data').free > RESOURCE['internal_launch_floor_bytes'], 'own launch floor')
    FIX.mkdir()
    group('frozen_actual_plan_author_closure_and_prior_B_review',bindings)
    group('actual_curl_readonly_owned_hostile_config_and_v1_argv',ambient_curl)
    group('independent_exact_rational_invented_storage_capacity',representation)
    group('independent_scalar_signed_mean_imputation_math',independent_arithmetic)
    for name, fault in [('healthy',None),('wrong_command','wrong_command'),('extra_directory','extra_directory'),('cleanup_error_reject','cleanup_error'),('postseal_resource_reject','postseal_resource'),('postseal_output_mutation_reject','postseal_output_mutation')]:
        group('actual_controller_private_metadata_'+name,lambda name=name,fault=fault:terminal_case(name,fault))
    files={str(p):sha(p) for p in FIX.rglob('*') if p.is_file()}
    receipt=dict(schema='INDEPENDENT_SIGNED_REFERENCE_SOURCE_PRELAUNCH_CONTROLS_V1',
                 actual_plan_sha256=PLAN_SHA,author_seal_sha256=AUTHOR_SHA,
                 cases=CASES,case_count=len(CASES),assertions=ASSERTIONS,
                 fixture_file_sha256=files,fixture_bytes=sum(Path(p).stat().st_size for p in files),
                 elapsed_seconds=time.monotonic()-START,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 actual_body_reads_or_network_or_workers_or_mutex_operations=0,
                 runtime_recensus=0,old_author_suites_rerun=0,estimator_calls=0)
    put(R/'independent_signed_reference_source_prelaunch_controls_v1_2.json',receipt)
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['cases','fixture_file_sha256']}))


if __name__=='__main__':main()
