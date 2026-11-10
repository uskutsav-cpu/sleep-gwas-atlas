"""Actual matched-trainer execution qualification; default is read-only readiness.

Six unmodified trainer processes exercise four updates versus two plus two for
both arms. Two separately attested resource processes exercise 1000 updates and
complete development inference. Every weight remains engineering-only.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

from permitted_AVI_experiment import ROOT,OUT,AVI,get_AVI_plan,readiness,sha


def live_fitting_processes():
    text = subprocess.run(['ps','-axo','pid,command'],check=True,capture_output=True,text=True).stdout
    names = {'run_resumable_comparison.py','readonly_native_runtime.py','readonly_native_runtime_768.py',
             'train_resumable_development.py','train_permitted_AVI_development.py','run_permitted_AVI_comparison.py',
             'audit_permitted_AVI_training.py','bench_permitted_AVI_resources.py'}
    result = []
    for row in text.splitlines()[1:]:
        fields = row.strip().split(None,1)
        if len(fields) != 2 or int(fields[0]) == os.getpid():
            continue
        parts = shlex.split(fields[1])
        if parts and Path(parts[0]).name.startswith('python') and any(Path(p).name in names for p in parts[1:]):
            result.append({'PID':int(fields[0]),'command':fields[1]})
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--execute',action='store_true')
    args = parser.parse_args()
    live = live_fitting_processes()
    if not args.execute:
        status = readiness()
        blockers = [s for s in status['blockers'] if not s.startswith('The new matched trainer still needs')]
        if live:
            blockers.append('A real fitting process is live; do not overlap this audit.')
        print(json.dumps({'execution_audit_ready':not blockers,'blockers':blockers,
                          'actual_live_fitting_processes':live,'actual_execution_audit_done':False}),flush=True)
        return
    assert not live, 'No concurrent fitting during the matched-trainer execution audit.'
    plan = get_AVI_plan(for_execution_audit=True)
    output = AVI/'matched_AVI_trainer_execution_audit.json'
    assert not output.exists(), 'Do not repeat an already completed execution audit.'
    runs = OUT/'AVI_recovery_audit_runs'
    assert not list(runs.glob('AVI_execution_*')) and not list(runs.glob('AVI_resources_*')), 'Inspect preserved prior audit attempts before any new execution.'
    logs = AVI/'execution_audit_logs';logs.mkdir(exist_ok=True)
    recipe = plan['recipes'][0]
    env = {**os.environ,'PYTHONDONTWRITEBYTECODE':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2',
           'VECLIB_MAXIMUM_THREADS':'2','PYTORCH_MPS_HIGH_WATERMARK_RATIO':'0.4','PYTORCH_MPS_LOW_WATERMARK_RATIO':'0.3'}
    work = Path(__file__).parent

    def run(arm,tag,updates,resume=False,resource_probe=False):
        if resource_probe:
            command = [sys.executable,'-u',str(work/'bench_permitted_AVI_resources.py'),'--arm',arm,'--audit-tag',tag]
        else:
            command = [sys.executable,'-u',str(work/'train_permitted_AVI_development.py'),
                       '--recipe',recipe,'--held-study','none','--arm',arm,
                       '--audit-tag',tag,'--audit-stop-after-updates',str(updates)] + (['--resume'] if resume else [])
        path = logs/(tag+('_resumed' if resume else '')+'.log')
        assert not path.exists()
        with path.open('w') as handle:
            subprocess.run(command,check=True,stdout=handle,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env)
        print('COMPLETED AVI ENGINEERING PROCESS',arm,tag,updates,'resume',resume,flush=True)

    def checkpoint(tag,arm):
        return runs/tag/'development'/f'{recipe}_none_{arm}_{plan["seed"]}_recovery.pt'

    started_metadata = {}
    for arm in plan['arms']:
        direct_tag = f'AVI_execution_{arm}_direct4'
        continued_tag = f'AVI_execution_{arm}_continued4'
        run(arm,direct_tag,4)
        run(arm,continued_tag,2)
        prefix = checkpoint(continued_tag,arm)
        preserved = runs/continued_tag/'committed_two_update_prefix.pt'
        os.link(prefix,preserved)  # The continuation uses atomic replacement.
        run(arm,continued_tag,4,resume=True)
        started_path = checkpoint(direct_tag,arm).with_name(checkpoint(direct_tag,arm).name.replace('_recovery.pt','_started.json'))
        started_metadata[arm] = json.loads(started_path.read_text())
    resources = {}
    for arm in plan['arms']:
        tag = f'AVI_resources_{arm}_v1'
        run(arm,tag,1000,resource_probe=True)
        path = runs/tag/'resource_probe_result.json'
        proof = json.loads(path.read_text());assert proof['passed'] and proof['actual_training_updates'] == 1000
        assert proof['maximum_RSS_bytes_on_macOS'] < 3*1024**3
        assert proof['minimum_free_disk_at_recorded_checkpoints'] > 300*1024**2
        profile_path = runs/tag/'resource_profile.json'
        profile = json.loads(profile_path.read_text())
        assert sha(profile_path) == proof['resource_profile_sha256']
        assert [row['updates'] for row in profile] == list(range(0,1001,100))
        late_growth = profile[-1]['maximum_RSS_bytes_on_macOS']-profile[5]['maximum_RSS_bytes_on_macOS']
        assert late_growth < 256*1024**2, 'Late resource growth is too large for scientific launch.'
        validation = proof['full_development_validation']
        assert validation['rows'] == started_metadata[arm]['prediction_rows']
        assert validation['prediction_indices_sha256'] == started_metadata[arm]['prediction_index_sha256']
        assert validation['native_lengths'] == [150,199,200,270]
        assert validation['parameters_and_buffers_unchanged'] and validation['CPU_and_MPS_RNG_unchanged']
        resources[arm] = {'result':proof,'result_sha256':sha(path),
                          'late_500_update_and_full_validation_peak_RSS_growth_bytes':late_growth}
    # Tensor verification runs only after every fitting process is terminal.
    verification_log = logs/'CPU_state_verification.log'
    assert not verification_log.exists()
    with verification_log.open('w') as handle:
        subprocess.run([sys.executable,'-u',str(work/'verify_permitted_AVI_execution_states.py')],
                       check=True,stdout=handle,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env)
    verification_path = AVI/'matched_AVI_execution_state_verification.json'
    state_verification = json.loads(verification_path.read_text())
    assert state_verification['passed']
    assert state_verification['code_sha256']['verify_permitted_AVI_execution_states.py'] == sha(work/'verify_permitted_AVI_execution_states.py')
    continuations = state_verification['continuation_results']
    result = {'passed':True,'time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'both_arms_executed':True,'common_initialization_and_batch_schedule_verified':True,
              'actual_saved_state_continuation_verified':True,'fitting_device_resource_preflight_passed':True,
              'score_features_reached_actual_training_gradients':True,'same_availability_shape_and_initial_parameters':True,
              'continuation_numeric_tolerance':{'absolute':1e-9,'relative':1e-6,'scope':'Four updates only; initial restore and all RNG/schedule states exact.'},
              'continuation_results':continuations,'resource_results':resources,
              'actual_unmodified_trainer_processes':6,'separately_attested_resource_processes':2,
              'separate_CPU_state_verification_processes':1,'state_verification_sha256':sha(verification_path),
              'coordinator_has_no_Torch_import':True,
              'actual_updates_audit_only':2016,'engineering_weights_used_for_scientific_initialization':False,
              'new_external_outcomes_read':False,'predictive_gain_established':False,
              'scope':'Actual runtime, paired initialization/data, recovery, gradients reaching permitted inputs, and full-development resource behavior. No biological success or outside-study evidence.',
              'code_sha256':{p.name:sha(p) for p in [Path(__file__),work/'bench_permitted_AVI_resources.py',
                                work/'verify_permitted_AVI_execution_states.py',
                                work/'train_permitted_AVI_development.py',work/'permitted_AVI_experiment.py',
                                work/'permitted_AVI_network.py',work/'permitted_AVI_head.py',work/'permitted_AVI_features.py',
                                work/'readonly_native_inputs.py',work/'resumption.py',work/'audit_readonly_runtime_recovery.py']}}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print('ACTUAL MATCHED AVI TRAINER EXECUTION AUDIT PASS',flush=True)


if __name__ == '__main__':
    main()
