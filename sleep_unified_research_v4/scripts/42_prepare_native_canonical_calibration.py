#!/usr/bin/env python3
"""Prepare (default) or later execute only 16 independently admitted controls.

Preparation imports no numerical runtime, opens no GWAS outcome columns, and
launches no estimator. An explicit reviewed admission bound to the plan and all
190 historical completion receipts is required for every future worker.
"""
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from canonical_calibration_common import ADMISSION, PACKAGE, PLAN, ROOT, SHARED_LOCK, admit, check_hashes, exclusive_heavy_lock, sha, utc, write_new

HISTORICAL = PACKAGE / 'manifests/ssd_native_execution_plan_v4_3.json'
BLOCKS = PACKAGE / 'manifests/canonical_200_interval_preparation_v4.json'
CONTROLS = PACKAGE / 'statistical_validation/genomicsem_native_calibration_controls_v4.tsv'
SSD = Path('/Volumes/Extreme SSD/sleep-unified-research-v1/research-completion-2026-10-09/sleep_unified_research_v4/statistical_validation/native_canonical_calibration_v4')
ALLOWED = ['CORE_INSOMNIA_BMI_DUPLICATE', 'FIRST_FROZEN_VALIDATION_PAIR_DUPLICATE']


def prepare():
    if PLAN.exists() or ADMISSION.exists():
        raise RuntimeError('EXISTING_FROZEN_PREPARATION_REQUIRES_REVIEW_NO_OVERWRITE')
    native = json.loads(HISTORICAL.read_text())
    blocks = json.loads(BLOCKS.read_text())
    if blocks['native_plan_sha256'] != sha(HISTORICAL) or blocks['build'] != 'GRCh37/hg19' or blocks['block_count'] != 200 or not blocks['all_source_reference_hashes_unchanged_after']:
        raise RuntimeError('CANONICAL_PARTITION_NOT_ADMITTED_FOR_PREPARATION')
    if len(blocks['reference_sha256']) != 22 or any(native['dependencies_sha256'].get(p) != h for p, h in blocks['reference_sha256'].items()):
        raise RuntimeError('CONSUMED_REFERENCE_IDENTITIES_DO_NOT_MATCH')
    if sha(blocks['interval_tsv_path']) != blocks['interval_tsv_sha256'] or sha(blocks['coordinate_map_path']) != blocks['coordinate_map_sha256']:
        raise RuntimeError('FROZEN_COORDINATES_OR_INTERVALS_CHANGED')
    verified = {r['path']: r['actual_sha256'] for r in native['inputs_verified'] if r['match']}
    with CONTROLS.open(newline='') as f:
        controls = list(csv.DictReader(f, delimiter='\t'))
    if [r['control_id'] for r in controls] != ALLOWED:
        raise RuntimeError('PRESPECIFIED_CONTROL_MEMBERSHIP_CHANGED')
    dependencies = dict(native['dependencies_sha256'])
    small_paths = [Path(__file__), PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py',
                   PACKAGE / 'scripts/canonical_calibration_common.py', PACKAGE / 'scripts/canonical_native_adapter.py',
                   PACKAGE / 'scripts/canonical_calibration_capture.py', PACKAGE / 'scripts/canonical_calibration_arithmetic.py',
                   PACKAGE / 'scripts/canonical_calibration_fixture.py',
                   PACKAGE / 'scripts/41_prepare_canonical_blocks.py', HISTORICAL, BLOCKS, CONTROLS,
                   PACKAGE / 'statistical_validation/common_boundary_adapter_prototype_v4.py',
                   PACKAGE / 'statistical_validation/canonical_native_calibration_fixture_receipt_v4_2.json',
                   Path(blocks['interval_tsv_path']), Path(blocks['coordinate_map_path'])]
    dependencies.update({str(p): sha(p) for p in small_paths})
    python = next(p for p in native['dependencies_sha256'] if p.endswith('/.ldsc-env/bin/python'))
    ldsc_dir = str(ROOT.parent / 'ldsc-code')
    refs = str(Path(next(iter(blocks['reference_sha256']))).parent) + '/'
    jobs = []
    source_hashes = {}
    for control in controls:
        inputs = [control['input1'], control['input2']]
        identities = {control['input1']: control['input1_sha256'], control['input2']: control['input2_sha256']}
        if any(verified.get(p) != h for p, h in identities.items()):
            raise RuntimeError('CONTROL_SOURCES_NOT_EXACTLY_HISTORICALLY_ADMITTED')
        source_hashes.update(identities)
        for method in ('one_step', 'two_step30'):
            for partition in ('stock_default200', 'canonical200'):
                for duplicate in ('A', 'B'):
                    job_id = control['control_id'] + '__' + method + '__' + partition + '__' + duplicate
                    out = SSD / job_id
                    ldsc_args = ['--rg', ','.join(inputs), '--ref-ld-chr', refs, '--w-ld-chr', refs,
                                 '--n-blocks', '200', '--print-delete-vals', '--out', str(out / 'stock_ldsc')]
                    if method == 'two_step30':
                        ldsc_args += ['--two-step', '30']
                    jobs.append({'job_id': job_id, 'control_id': control['control_id'], 'method': method, 'partition': partition,
                                 'duplicate': duplicate, 'inputs': inputs, 'input_sha256': identities, 'output_dir': str(out), 'ldsc_args': ldsc_args})
    # Compressed input identity reads only; no GWAS rows/phenotype results read.
    current_inputs = check_hashes(source_hashes)
    if len(jobs) != 16 or len(current_inputs) != 3:
        raise RuntimeError('EXPECTED_16_CONTROLS_AND_THREE_UNIQUE_SOURCES')
    record = {'schema': 'native_canonical_calibration_operational_plan_v1', 'prepared_utc': utc(), 'self_path': str(PLAN),
              'scope': 'METHOD_CALIBRATION_ONLY', 'allowed_control_ids': ALLOWED, 'jobs': jobs,
              'historical_plan_path': str(HISTORICAL), 'historical_plan_sha256': sha(HISTORICAL),
              'block_receipt_path': str(BLOCKS), 'block_receipt_sha256': sha(BLOCKS),
              'build': blocks['build'], 'interval_tsv_sha256': blocks['interval_tsv_sha256'],
              'ordered_all_reference_records_sha256': blocks['ordered_all_reference_records_sha256'],
              'ordered_construction_records_sha256': blocks['ordered_construction_records_sha256'],
              'coordinate_map_path': blocks['coordinate_map_path'], 'coordinate_map_sha256': blocks['coordinate_map_sha256'],
              'partition_amendment': blocks['partition_amendment'], 'static_construction_exclusions': blocks['static_construction_exclusions'],
              'reference_sha256': blocks['reference_sha256'], 'dependencies_sha256': dependencies,
              'inputs_current_sha256': current_inputs, 'python': python, 'ldsc_dir': ldsc_dir, 'ssd_output_root': str(SSD),
              'environment': native['environment'], 'worker_count': 1, 'blas_threads': 1,
              'shared_lock_path': str(SHARED_LOCK), 'shared_lock_policy': 'LOCK_EX|LOCK_NB for entire executor lifetime, inherited descriptor required by every native/arithmetic worker.',
              'internal_floor_bytes': 3 * 1024**3, 'ssd_floor_bytes': 5 * 1024**3, 'worker_rss_limit_bytes': 2 * 1024**3,
              'output_limit_bytes': 8 * 1024**3, 'per_worker_seconds_limit': 2 * 3600, 'runtime_poll_seconds': 2,
              'resource_storage_assumption': '16 controls retain exact weighted arrays; 8 GiB total cap, plus 5 GiB SSD free floor; admission review must assess capacity before launch.',
              'operational_resource_amendment': 'Capture output cap raised from prototype 4 GiB to 8 GiB to retain all 16 prespecified control designs/deletes; this separate operational amendment requires root review before execution. Scientific thresholds/inputs/partition remain unchanged.',
              'arithmetic': {'point_delete_atol': 1e-10, 'point_delete_rtol': 1e-8, 'covariance_atol': 1e-14, 'covariance_rtol': 1e-6,
                             'duplicate_zero_delta_se_atol': 1e-12, 'independent_row_solves': 'replicate A of all 8 control/arm/partition combinations',
                             'full_fit_weights_and_Nbar_frozen_for_deletion': True, 'no_calibrated_difference_p_values': True},
              'native_jobs_launched': 0, 'GWAS_outcome_columns_read': False, 'allow_41_covariance_outcomes': False,
              'allow_calibrated_biological_p_values': False, 'realistic_LD_sampling_calibration_pass': False,
              'remaining_method_gate': 'Verified source-matched signed LD/genotype asset and frozen realistic-LD simulations; synthetic/duplicate arithmetic cannot establish real genomic correlation length, selection, or reference validity.'}
    write_new(PLAN, record)
    write_new(ADMISSION, {'schema': 'native_canonical_calibration_execution_admission_v1', 'prepared_utc': utc(),
                         'plan_sha256': sha(PLAN), 'scope': 'METHOD_CALIBRATION_ONLY', 'execution_admitted': False,
                         'independent_binding_review_pass': False, 'resource_plan_review_pass': False,
                         'independent_review_path': None, 'independent_review_sha256': None,
                         'historical_campaign_receipt_sha256': {}, 'allow_41_covariance_outcomes': False,
                         'allow_calibrated_biological_p_values': False,
                         'blocked_reason': 'Historical 190-job campaign active; independent executor binding/resource review and exact completion receipts required before any fit.'})
    print(json.dumps({'status': 'PREPARED_EXECUTION_BLOCKED', 'jobs': len(jobs), 'plan_sha256': sha(PLAN), 'native_jobs_launched': 0}))


def monitor_helpers():
    path = PACKAGE / 'scripts/30_prepare_and_run_ssd_native_campaign.py'
    spec = importlib.util.spec_from_file_location('canonical_calibration_owned_process_helpers', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def snapshot(output_root):
    return {'recorded_utc': utc(), 'internal_free_bytes': shutil.disk_usage('/System/Volumes/Data').free,
            'ssd_free_bytes': shutil.disk_usage(output_root.parent).free}


def output_bytes(output_root):
    return sum(p.stat().st_size for p in output_root.rglob('*') if p.is_file())


def monitored_worker(plan, admission_sha, job_id, command, output_dir, lock_fd):
    helper = monitor_helpers()
    root = Path(plan['ssd_output_root'])
    state = snapshot(root)
    if state['internal_free_bytes'] < plan['internal_floor_bytes'] or state['ssd_free_bytes'] < plan['ssd_floor_bytes']:
        raise RuntimeError('FULL_NATIVE_RESOURCE_FLOOR_NO_WORKER_LAUNCHED')
    env = {**os.environ, 'TMPDIR': str(root / 'tmp'), 'XDG_CACHE_HOME': str(root / 'cache'),
           'MPLCONFIGDIR': str(root / 'cache/matplotlib'), 'PYTHONDONTWRITEBYTECODE': '1',
           'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'}
    for path in ('tmp', 'cache/matplotlib'):
        (root / path).mkdir(parents=True, exist_ok=True)
    logfile = output_dir / 'worker_stdout.log'
    proc = None
    started = time.monotonic()
    peak = 0
    reason = None
    teardown = None
    returncode = None
    samples = []
    try:
        with logfile.open('xb') as stream:
            proc = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True, pass_fds=(lock_fd,))
            while True:
                returncode = proc.poll()
                rss, pids = helper.owned_rss(proc.pid)
                peak = max(peak, rss)
                state = snapshot(root)
                size = output_bytes(root)
                elapsed = time.monotonic() - started
                if state['internal_free_bytes'] < plan['internal_floor_bytes']: reason = 'INTERNAL_RESOURCE_FLOOR'
                elif state['ssd_free_bytes'] < plan['ssd_floor_bytes']: reason = 'SSD_RESOURCE_FLOOR'
                elif rss > plan['worker_rss_limit_bytes']: reason = 'OWNED_WORKER_RSS_LIMIT'
                elif size > plan['output_limit_bytes']: reason = 'OUTPUT_LIMIT'
                elif elapsed > plan['per_worker_seconds_limit']: reason = 'WORKER_DEADLINE'
                elif sha(ADMISSION) != admission_sha: reason = 'ADMISSION_CHANGED_DURING_RUN'
                if len(samples) < 100 or reason or returncode is not None:
                    samples.append({**state, 'elapsed_seconds': elapsed, 'rss_bytes': rss, 'owned_pids': pids, 'output_bytes': size})
                if reason or returncode is not None: break
                time.sleep(plan['runtime_poll_seconds'])
    finally:
        if proc is not None:
            teardown = helper.terminate_owned(proc)
            returncode = proc.returncode
        receipt = {'job_id': job_id, 'command': command, 'completed_utc': utc(), 'plan_sha256': sha(PLAN),
                   'admission_sha256': admission_sha, 'returncode': returncode, 'stop_reason': reason,
                   'elapsed_seconds': time.monotonic() - started, 'peak_owned_rss_bytes': peak,
                   'resource_samples': samples, 'process_group_teardown': teardown,
                   'all_output_sha256': {str(p): sha(p) for p in output_dir.rglob('*') if p.is_file()}}
        write_new(output_dir / 'monitor_receipt.json', receipt)
    if returncode != 0 or reason:
        raise RuntimeError('CALIBRATION_WORKER_FAILED_PRESERVED_PARTIAL_OUTPUT: ' + job_id)


def execute(lock_fd):
    plan, _ = admit()
    admission_sha = sha(ADMISSION)
    root = Path(plan['ssd_output_root'])
    if root.exists():
        raise RuntimeError('EXISTING_CALIBRATION_NAMESPACE_REQUIRES_REVIEW_NO_OVERWRITE')
    root.mkdir(parents=True, exist_ok=False)
    for job in plan['jobs']:
        # Recheck admission/history before EACH worker, never overlap native work.
        admit()
        out = Path(job['output_dir'])
        out.mkdir(exist_ok=False)
        command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_capture.py'),
                   '--plan', str(PLAN), '--admission', str(ADMISSION), '--plan-sha', sha(PLAN), '--job-id', job['job_id'], '--shared-lock-fd', str(lock_fd)]
        monitored_worker(plan, admission_sha, job['job_id'], command, out, lock_fd)
        receipt = json.loads((out / 'native_capture_receipt.json').read_text())
        if receipt['job'] != job or receipt['plan_sha256'] != sha(PLAN):
            raise RuntimeError('NATIVE_CAPTURE_RECEIPT_NOT_BOUND')
        print(json.dumps({'completed_method_control': job['job_id'], 'biological_outcomes_computed': 0}), flush=True)
    verify(lock_fd)


def verify(lock_fd):
    plan, _ = admit()
    for job in plan['jobs']:
        out = Path(job['output_dir'])
        r = json.loads((out / 'monitor_receipt.json').read_text())
        if r['returncode'] != 0 or r['stop_reason'] is not None or r['process_group_teardown']['remaining_group_members'] != []:
            raise RuntimeError('ALL_16_NATIVE_CONTROLS_MUST_COMPLETE_BEFORE_ARITHMETIC')
        check_hashes(r['all_output_sha256'])
        if not (out / 'native_capture_receipt.json').exists():
            raise RuntimeError('MISSING_NATIVE_CAPTURE')
    out = Path(plan['ssd_output_root']) / 'independent_arithmetic_worker'
    out.mkdir(exist_ok=False)
    command = [plan['python'], '-B', str(PACKAGE / 'scripts/canonical_calibration_arithmetic.py'), '--plan', str(PLAN), '--admission', str(ADMISSION), '--shared-lock-fd', str(lock_fd)]
    monitored_worker(plan, sha(ADMISSION), 'INDEPENDENT_CAPTURED_ROW_ARITHMETIC', command, out, lock_fd)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=('prepare', 'execute', 'verify'), nargs='?', default='prepare')
    args = p.parse_args()
    if args.mode == 'prepare':
        prepare()
    else:
        with exclusive_heavy_lock() as lock_fd:
            {'execute': execute, 'verify': verify}[args.mode](lock_fd)


if __name__ == '__main__':
    main()
